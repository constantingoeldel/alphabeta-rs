use std::{fs, sync::Mutex};

use crate::{
    config::set,
    files::{file_name, lines_from_file, load_methylome, open_file},
    genes::{Gene, Genome, Region},
    setup::setup_output_dir,
    windows::Windows,
    *,
};
use alphabeta::Analysis;
use itertools::Itertools;
// use methylome::genes::{Gene, Genome};
// use methylome::methylation_site::Chromosome;

use methylome::{steady_state, Chromosome};
use progress_bars::MultiProgress;
use rayon::prelude::*;

pub fn extract(args: config::Metaprofile) -> Return<()> {
    set(args.clone());
    let start = std::time::Instant::now();

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
        return Err(Error::NoGenesFound);
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
        .try_for_each_with(genome, |genome, path| -> Return<()> {
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

    // Removing the mutexes as the parallel part is over
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

    println!("Created metaprofile in: {:?}", start.elapsed());

    // If the alphabeta subcommand is set, run the alphabeta method
    if args.command.is_some() {
        alphabeta_multiple(args, max_gene_length, distributions[0].clone())
    }

    Ok(())
}

fn alphabeta_multiple(args: Metaprofile, max_gene_length: u32, distribution: Vec<i32>) {
    let regions = vec![
        (Region::Upstream, args.cutoff),
        (Region::Gene, max_gene_length),
        (Region::Downstream, args.cutoff),
    ];
    let mut results = Vec::new();
    let total_steps = if args.absolute {
        (max_gene_length + 2 * args.cutoff) / args.window_step
    } else {
        (3 * 100) / args.window_step
    };

    let (multi, pb) = progress_bars::multi(total_steps as u64);

    // Create an empty 3D array to store the raw results
    for region in regions {
        let max = if args.absolute { region.1 } else { 100 };

        for window in (0..max).step_by(args.window_step as usize) {
            pb.inc(1);

            let args = alphabeta::AlphaBeta::default(
                args.output_dir
                    .join(region.0.to_string())
                    .join(window.to_string()),
                args.iterations,
            );

            let alphabeta_result = alphabeta::run(args, &multi);
            match alphabeta_result {
                Err(e) => println!("Error: {e}"),
                Ok((model, analysis, _, obs_meth_lvl)) => {
                    results.push((model, analysis, region.0.clone(), obs_meth_lvl));
                }
            }
        }
    }
    pb.finish();
    let mut print = String::from("run;window;cg_count;region;alpha;beta;1/2*(alpha+beta);pred_steady_state;obs_steady_state;sd_alpha;sd_beta;ci_alpha_0.025;ci_alpha_0.975;ci_beta_0.025;ci_beta_0.975\n");
    for (i, ((model, analysis, region, obs_meth_lvl), d)) in
        results.iter().zip(distribution.iter()).enumerate()
    {
        print += &format!(
            "{};{};{};{};{};{};{};{};{};{};{};{};{};{};{}\n",
            args.name,
            i,
            d,
            region,
            model.alpha,
            model.beta,
            0.5 * (model.alpha + model.beta),
            steady_state(model.alpha, model.beta),
            obs_meth_lvl,
            analysis.sd_alpha,
            analysis.sd_beta,
            analysis.ci_alpha.0,
            analysis.ci_alpha.1,
            analysis.ci_beta.0,
            analysis.ci_beta.1
        )
    }

    println!("{print}");

    fs::write(args.output_dir.join("results.txt"), print).expect("Could not save results to file.");
    // let db = db::connect()
    //     .await
    //     .expect("Could not connect to database: Did you provide a connection string?");
    // import_results(&db, args.name, results).await.expect("Could not save results to a database. Your data is stored in files in each directory");
    let analyses = results
        .iter()
        .map(|r| r.1.clone())
        .collect::<Vec<Analysis>>();

    plot::metaplot(&analyses, &args).expect("Could not plot results");
}
