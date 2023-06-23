# Alphabeta (but it's fast)

This repostiory contains a fast implementation of the AlphaBeta algorithm first proposed by [Yadollah Shahryary, Frank Johannes and Rashmi Hazarika](https://doi.org/10.5281/zenodo.3992612). The original R implemtation is accessible on Github [here](https://github.com/jlab-code/AlphaBeta/). I matched the original parameters needed for the program, to the old [documentation](https://github.com/jlab-code/AlphaBeta/blob/master/vignettes/AlphaBeta.pdf) is still useful.

Additionally, this repository also contains a program for creating metaprofiles of (epi)genetic data. They are connected, so output from the metaprofile program can automcatically be fed into AlphaBeta.

## How to use

You'll need to install [Rust](https://www.rust-lang.org/tools/install), the programming language used for this project. Then enter the following commands into your terminal:

```bash
cargo install alphabeta
```

This will install the programs on your system. You can then ensure everything works by running:

```bash
alphabeta --help
# and
metaprofile --help
```

### Windows

Windows is currently not supported, but using [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) on Windows works.

Open powershell and run:

```powershell
wsl --install -d ubuntu
```

Configure your user, open Ubuntu and then install the dependencies and the program:

```bash
sudo apt update
sudo apt install build-essential pkg-config libssl-dev
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install alphabeta
```


<details>
<summary>I get an error!</summary>

If you received an error message about libopenblas, you will need to run the code from inside the repository with cargo (I don't really understand this issue)

```bash
cargo run --release --bin alphabeta
# or
cargo run --release --bin metaprofile
```

</details>

### Updating

If you want to use a new version of the program, just run the following command again:

```bash
cargo install alphabeta
```

## Parameters

### Alphabeta

```bash
Usage: alphabeta [OPTIONS] --edges <EDGES> --nodes <NODES> --output <OUTPUT>

Options:
  -i, --iterations <ITERATIONS>
          Number of iterations to run for Nelder-Mead optimization, even 100 is enough [default: 1000]
  -e, --edges <EDGES>
          Relative or absolute path to an edgelist, see /data for an example
  -n, --nodes <NODES>
          Relative or absolute path to a nodelist, see /data for an example
  -p, --posterior-max-filter <POSTERIOR_MAX_FILTER>
          Minimum posterior probability for a singe basepair read to be included in the estimation [default: 0.99]
  -o, --output <OUTPUT>
          Relative or absolute path to an output directory, must exist, EXISTING FILES WILL BE OVERWRITTEN
  -h, --help
          Print help
  -V, --version
          Print version
```

### Metaprofile

```bash
Usage: metaprofile [OPTIONS] --methylome <METHYLOME> --genome <GENOME> --output-dir <OUTPUT_DIR>

Options:
  -m, --methylome <METHYLOME>      Path of directory containing the methlyome files from which to extract the CG-sites
  -g, --genome <GENOME>            Path of the annotation file containing information about beginning and end of gbM-genes
  -w, --window-size <WINDOW_SIZE>  Size of the window in percent of the gbM-gene length or in basepair number if --absolute is supplied [default: 5]
  -s, --window-step <WINDOW_STEP>  Size of the step between the start of each window. Default value is window-size, so no overlapp happens [default: 0]
  -o, --output-dir <OUTPUT_DIR>    Path of the directory where extracted segments shall be stored
  -a, --absolute                   Use absolute length in base-pairs for window size instead of percentage of gene length
  -c, --cutoff <CUTOFF>            Number of basepairs to include upstream and downstream of gene [default: 2048]
  -i, --invert                     Invert strands, to switch from 5' to 3' and vice versa
      --db                         Use a Postgres database to do everything
  -e, --edges <EDGES>              Provide an edgefile
  -n, --nodes <NODES>              Provide a nodefile - paths will be updated to match the output directory
      --alphabeta                  Also run AlphaBeta on every window after extraction, results will be stored in the same directory as the segments
      --name <NAME>                Name of the run to be used when storing the result in Postgres [default: "Anonymous Run 1684308724"]
  -f, --force                      Overwrite existing content in output directory? If false (default) it will reuse existing windows
      --cutoff-gene-length         Let the cutoff be the gene length instead of a fixed number. So if the gene is 1000 bp long, the cutoff will be 1000 bp instead of 2048 bp (the default). This option takes preference over the cutoff option
  -h, --help                       Print help
  -V, --version                    Print version
```

About window step and size:

Size determines the "length" of each window, for example for `--window-size 5`, each window will span 5% of the length of the gene it is in. If you supply `--absolute`, the size will be interpreted as the number of basepairs instead of a percentage, so 5 bp.

Step determines the distance between the start of each window. If you supply `--window-step 1 and --window-size 5`, the first window will go from 0% to 5% and the second from 1% to 6% and so on. If you supply `--window-step 5 and --window-size 5`, the first window will go from 0% to 5% and the second from 5% to 10% and so on. In the latter case, you can also omit the `step` paramter, as it will default to the same value as `size`.

## Examples

### Run alphabeta

```bash
alphabeta \
--edges ./data/edgelist.txt \
--nodes ./data/nodelist.txt \
--output ./data/output
```

### Create a metaprofile and feed it into AlphaBeta

```bash
metaprofile \
--methylome ../methylome/within_gbM_genes/ \
--genome ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
--output-dir /mnt/extStorage/workingDir/./windows/wt \
--edges ../methylome/edgelist.txt \
--nodes ../methylome/nodelist.txt \
--alphabeta \
--name wildtype \
--window-step 1 --window-size 5
```
