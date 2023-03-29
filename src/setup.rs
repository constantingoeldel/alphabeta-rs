use std::fs;

use crate::*;

pub fn set_up_output_dir(args: arguments::Windows, max_gene_length: u32) -> Result<()> {
    fs::read_dir(&args.output_dir).or(Err(Error::File(args.output_dir.clone())))?; // Throw error if base output dir does not exist

    // // Replace existing content of output dir
    // fs::remove_dir_all(&args.output_dir).unwrap();
    // fs::create_dir(&args.output_dir).unwrap();

    if args.alphabeta {
        // If alphabeta is set, a nodelist and edgelist are required, this is checked for in lib.rs
        let edgelist = fs::read_to_string(
            args.edges
                .expect("Edgelist must be provided when alphabeta runs"),
        )?;
        let nodes = fs::read_to_string(
            args.nodes
                .expect("Nodelist must be provided when alphabeta runs"),
        )?;

        let sides = vec![
            ("upstream", args.cutoff),
            ("gene", max_gene_length),
            ("downstream", args.cutoff),
        ];

        for side in sides {
            let max = if args.absolute { side.1 } else { 100 };
            let side = side.0;

            for window in (0..max).step_by(args.window_step as usize) {
                let mut nodelist = String::new();
                let lines = nodes.split('\n');
                for line in lines {
                    if line.starts_with('/') {
                        let old_file = line.split('\t').next().unwrap();
                        let filename = old_file.split('/').last().unwrap();
                        let file = format!(
                            "{}/{}/{}/{}",
                            &args.output_dir.to_string_lossy(),
                            side,
                            window,
                            filename
                        );
                        nodelist += &line.replace(old_file, &file);
                    } else {
                        nodelist += line;
                    }
                    nodelist += "\n";
                }

                let path = format!("{}/{}/{}", &args.output_dir.to_string_lossy(), side, window);
                let window_dir = fs::read_dir(&path);

                match window_dir {
                    Ok(_) => {
                        print_dev!("Directory already exists, skipping: {}", path);
                    }
                    Err(_) => {
                        fs::create_dir_all(&path).expect("Could not create output directory");
                        fs::write(path.to_owned() + "/nodelist.txt", nodelist)
                            .expect("Nodelist not writable at ");
                        fs::write(path.to_owned() + "/edgelist.txt", &edgelist).expect("msg");
                    }
                };
            }
        }
    }
    Ok(())
}
