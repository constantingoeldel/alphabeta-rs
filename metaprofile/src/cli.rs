use metaprofile::{extract, Metaprofile as Args, Subcommands};

use clap::Parser;
use methylome::steady_state;
use ndarray::{Array, Axis};
use ndarray_npy::write_npy;
use std::fs;

fn main() {
    let mut args = Args::parse();

    println!("Starting run {}", args.name);

    // TODO validate you have to pass nodes and edges with alphabeta subcommand
    // match args.command {
    //     Some(Subcommands::AlphaBeta(s)) => {
    //         if s.nodes.is_none() || s.edges.is_none() {
    //             panic!(
    //                 "You need to specify the nodes and edges file when using the alphabeta method"
    //             );
    //         }
    //     }
    // }

    if args.window_step == 0 {
        args.window_step = args.window_size;
    }
    // move to lib.rs, rename to metaprofile
    let result = extract(args.clone());

    match (result) {
        Err(e) => println!("Error: {e}"),
        _ => println!("Done"),
    }
}

// #[cfg(test)]
// mod test {
//     use std::path::PathBuf;

//     use alphabeta::arguments::Windows;

//     #[test]
//     fn end_to_end() {
//         let args = Windows {
//             name: String::from("test"),
//             genome: PathBuf::from("data/annotation.bed"),
//             methylome: PathBuf::from("data/methylome"),
//             output_dir: String::from("data/output_metaplot"),
//             alphabeta: false,
//         };
//         }
//     }
// }
