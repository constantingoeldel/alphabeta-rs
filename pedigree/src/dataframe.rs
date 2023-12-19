use std::{
    fs::File,
    io::{BufRead as _, BufReader},
    path::Path,
};

use polars::prelude::*;

pub fn load<P>(path: P, custom_column_names: Option<&[&'static str]>) -> DataFrame
where
    P: AsRef<Path>,
{
    let tic = std::time::Instant::now();

    // Get first line of file
    let mut reader = BufReader::new(File::open(&path).unwrap());
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();

    // Check if file has header
    let has_header = line
        .split('\t')
        .next()
        .map(|s| s.parse::<u32>().is_err())
        .unwrap_or(false);
    let separator = if line.contains('\t') {
        "\t"
    } else if line.contains(',') {
        ","
    } else if line.contains(';') {
        ";"
    } else {
        " "
    };

    let mut df = CsvReader::from_path(path.as_ref())
        .unwrap()
        .has_header(has_header)
        .with_separator(separator.as_bytes()[0])
        .with_ignore_errors(true)
        .finish()
        .unwrap();

    if let Some(names) = custom_column_names {
        df.set_column_names(names).unwrap();
    } else if !has_header {
        let default_names = [
            "chromosome",
            "location",
            "strand",
            "context",
            "counts.methylated",
            "counts.total",
            "posteriorMax",
            "status",
            "rc.meth.lvl",
            "context.trinucleotide",
        ];
        let current_column_names: Vec<String> = df
            .get_column_names()
            .into_iter()
            .map(|s| s.into())
            .collect();
        current_column_names.iter().enumerate().for_each(|(i, s)| {
            df.rename(s, default_names.get(i).unwrap_or(&s.as_str()))
                .unwrap();
        });
    } else if df.get_column_index("seqnames").is_some() {
        df.rename("seqnames", "chromosome").unwrap();
    }

    let df = if df.get_column_index("chromosome").is_some() {
        df.filter(&df.column("chromosome").unwrap().is_not_null())
            .unwrap()
    } else {
        df
    };

    let colums_to_drop: Vec<String> = df
        .get_column_names()
        .iter()
        .filter(|c| c.starts_with('_'))
        .into_vec();

    let df = df.drop_many(&colums_to_drop);

    // if df.get_column_index("strand").is_some() {
    //     df.with_column(
    //         df.column("strand")
    //             .unwrap()
    //             .cast(&DataType::Categorical(None))
    //             .unwrap(),
    //     )
    //     .unwrap();
    // }

    // if df.get_column_index("status").is_some() {
    //     df.with_column(
    //         df.column("status")
    //             .unwrap()
    //             .cast(&DataType::Categorical(None))
    //             .unwrap(),
    //     )
    //     .unwrap();
    // }

    println!("Elapsed time: {:?}", tic.elapsed());
    df
}

fn sites_within(df: DataFrame, chromosome: u32, start: u32, end: u32, strand: &str) -> DataFrame {
    let cx = if df.get_column_index("location").is_some() {
        "location"
    } else {
        "start"
    };

    assert!(df.get_column_index(cx).is_some());
    assert!(df.get_column_index("strand").is_some());
    assert!(df.get_column_index("chromosome").is_some());

    df.lazy()
        .filter(
            col(cx)
                .gt(start)
                .and(col("chromosome").eq(lit(chromosome)))
                .and(col(cx).lt(end))
                .and(col("strand").eq(lit(strand))),
        )
        .collect()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn load_small() {
        println!("Hello, world!");

        let path = PathBuf::from("../data/methylome/G0.txt");

        let df = load(path, None);

        assert_eq!(df.shape(), (2000, 10));
    }

    #[test]
    fn sites_within_gene() {
        let path = PathBuf::from("../data/methylome/G0.txt");

        let df = load(path, None);

        let df = sites_within(df, 1, 100, 200, "+");
        assert_eq!(df.shape(), (12, 10));
    }
    #[test]
    #[allow(non_snake_case)]
    fn sites_within_gbM() {
        let path =
            PathBuf::from("../data/methylome_full/within_gbM_genes/methylome_Col0_G0_All.txt");

        let df = load(path, None);

        assert_eq!(df.shape(), (950044, 9)); // context.trinucleotide is missing in this dataset
    }
    #[test]
    fn sites_huge() {
        let path = PathBuf::from(
            "../../MA3_new_total_original_methylome/methylome/methylome_Col0_G0_All.txt",
        );

        let df = load(path, None);

        assert_eq!(df.shape(), (42859516, 10));
    }
    #[test]
    fn sites_histone_mods() {
        let path = PathBuf::from("../../modifications/final/h2az.bed");

        let df = load(path, Some(&["chromosome", "start", "end", "_id", "_score"]));

        assert_eq!(df.shape(), (33799373, 3));
    }

    #[test]
    fn avg_methylation_lvl_by_chromosome() {
        let path = PathBuf::from("../data/methylome/G0.txt");

        let df = load(path, None);
        let df = df
            .lazy()
            .group_by(["chromosome"])
            .agg([col("rc.meth.lvl").mean().alias("avg_methylation_lvl")])
            .collect()
            .unwrap();
        println!("{}", df);
        assert_eq!(df.shape(), (5, 2));
    }

    // #[test]
    // fn convert_status_to_num() {
    //     let path = PathBuf::from("../data/methylome/G0.txt");

    //     let mut df = load(path, None);
    //     println!("{}", df);
    //     let df = df
    //         .apply("status", |s| {
    //             s.categorical()
    //                 .unwrap()
    //                 .iter_str()
    //                 .map(|l| match l {
    //                     Some("U") => Some(0.0),
    //                     Some("M") => Some(1.0),
    //                     Some("I") => Some(0.5),
    //                     _ => None,
    //                 })
    //                 .collect::<Float32Chunked>()
    //                 .into_series()
    //         })
    //         .unwrap();
    //     println!("{}", df);
    //     assert_eq!(df.shape(), (2000, 10));
    // }
}
