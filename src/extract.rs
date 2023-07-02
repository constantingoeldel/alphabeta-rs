use std::{fs, sync::Mutex};

use anyhow::bail;
use indicatif::MultiProgress;
use itertools::Itertools;

use crate::{
    files::{file_name, lines_from_file, load_methylome, open_file},
    genes::{Gene, Genome},
    methylation_site::Chromosome,
    setup::setup_output_dir,
    windows::Windows,
    *,
};
use rayon::prelude::*;

pub fn extract(args: arguments::Windows) -> Result<(u32, Vec<i32>)> {
    let start = std::time::Instant::now();
    let mut args = args;

    if args.alphabeta && (args.nodes.is_none() || args.edges.is_none()) {
        panic!("You need to specify the nodes and edges file when using the alphabeta method");
    }

    // Adj ust window_step to default value
    if args.window_step == 0 {
        args.window_step = args.window_size;
    }

    let methylome_files = load_methylome(&args.methylome)?;
    let annotation_lines = lines_from_file(&args.genome)?;
    let mut genes: Vec<Gene> = Vec::new();

    // Parse annotation file to extract genes
    for line in annotation_lines {
        let line = line?;
        let gene = Gene::from_annotation_file_line(&line, args.invert);
        if let Some(gene) = gene {
            genes.push(gene)
        }
    }

    if genes.is_empty() {
        bail!(Error::Simple("Could not parse a single annotation from the annotation file. Please check your input or add a parser implemenation for your data format."))
    }

    // Extract all the different chromosomes, which are used to construct the genome HashMap
    let chromosomes: Vec<&Chromosome> = genes.iter().map(|g| &g.chromosome).unique().collect();
    // Determine the maximum gene length by iterating over all genes
    let mut max_gene_length: u32 = 100; // if not using absolute window sizes, the maximum gene length will be 100%
    if args.absolute {
        max_gene_length = genes
            .iter()
            .map(|g| g.end - g.start)
            .max()
            .expect("There are no genes in the annotation file");
        println!("The maximum gene length is {max_gene_length} bp");
    }
    setup_output_dir(args.clone(), max_gene_length)?;

    // Structure genes first by chromosome
    let mut genome = Genome::new(chromosomes);
    for gene in genes {
        genome.insert_gene(gene.chromosome.clone(), gene);
    }
    genome.sort();

    let distributions = Mutex::new(Vec::new());
    let steady_state_methylations = Mutex::new(Vec::new());

    let bars = MultiProgress::new();

    methylome_files
        .par_iter()
        .try_for_each_with(genome, |genome, path| -> Result<()> {
            let file = open_file(path)?;
            let mut windows = Windows::extract(
                file,
                genome.to_owned(),
                max_gene_length,
                args.clone(),
                file_name(path),
                &bars,
            )?;
            if args.invert {
                windows = windows.inverse();
            }
            windows.save(args.clone(), file_name(path))?;
            let distribution = windows.distribution();

            distributions.lock().unwrap().push(distribution);

            let methylation = windows.steady_state_methylation();

            steady_state_methylations.lock().unwrap().push(methylation);

            Ok(())
        })?;

    // Removing the mutexes as the paralell part is over
    let steady_state_meth = steady_state_methylations.into_inner().unwrap();
    let distributions = distributions.into_inner().unwrap();
    let sample_size = steady_state_meth.len();
    let mut average_methylation = vec![0.0; steady_state_meth[0].len()];
    for source in steady_state_meth.iter() {
        for (i, window) in source.iter().enumerate() {
            average_methylation[i] += window / sample_size as f64
        }
    }

    let methylation_file = format!(
        "{}/steady_state_methylation.txt",
        &args.output_dir.display()
    );
    let all_methylations_file = format!(
        "{}/all_steady_state_methylation.txt",
        &args.output_dir.display()
    );
    let all_distributions_file = format!("{}/distributions.txt", &args.output_dir.display());

    let names = methylome_files.iter().map(|f| file_name(f)).collect();

    for (distribution, file) in distributions.iter().zip(&methylome_files) {
        fs::write(
            format!(
                "{}/distribution_{}",
                &args.output_dir.display(),
                file_name(file)
            ),
            Windows::print_distribution(distribution),
        )
        .unwrap();
    }

    fs::write(
        methylation_file,
        Windows::print_steady_state_methylation(&average_methylation),
    )?;
    fs::write(
        all_methylations_file,
        Windows::print_all_steady_state_methylations(names, steady_state_meth),
    )?;

    fs::write(
        all_distributions_file,
        Windows::print_all_distributions(
            methylome_files.iter().map(|f| file_name(f)).collect(),
            &distributions,
        ),
    )?;

    println!("Done in: {:?}", start.elapsed());
    Ok((max_gene_length, distributions[0].clone()))
}
