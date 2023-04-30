# Alphabeta (but it's fast)

This repostiory contains a fast implementation of the AlphaBeta algorithm first proposed by [Yadollah Shahryary, Frank Johannes and Rashmi Hazarika](https://doi.org/10.5281/zenodo.3992612). The original R implemtation is accessible on Github [here](https://github.com/jlab-code/AlphaBeta/). I matched the original parameters needed for the program, to the old [documentation](https://github.com/jlab-code/AlphaBeta/blob/master/vignettes/AlphaBeta.pdf) is still useful.

Additionally, this repository also contains a program for creating metaprofiles of (epi)genetic data. They are connected, so output from the metaprofile program can automcatically be fed into AlphaBeta.

## How to use

### Dependencies

The program depends on [OpenBLAS](https://www.openblas.net/) for fast matrix calculations. On Linux, you can install them with your package manager, on MacOS you can use [Homebrew](https://brew.sh/). On Windows, you can download prebuilt binaries from the [OpenBLAS website](https://www.openblas.net/).

Original Installation instructions [here](https://github.com/xianyi/OpenBLAS/wiki/Precompiled-installation-packages).

Debian/Ubuntu:

```bash
sudo apt update
sudo apt install libopenblas-dev
```

If you don't have sudo rights, ask your system administrator.

MacOS:

```bash
brew install openblas
```

Windows:

Download the prebuilt binaries from the [OpenBLAS website](https://www.openblas.net/) (Big .zip button). For building from source, see [this](https://github.com/blas-lapack-rs/openblas-src#windows-and-vcpkg).

<!-- ### Prebuilt binaries

In the releases tab, you can find prebuilt binaries for Windows, Linux and MacOS. These should work out of the box, but if they don't, please [open an issue](https://github.com/constantingoeldel/alphabeta-rs/issues/new). You'll need to download the program for your platform, unzip it and then run it from the command line. -->

### Building from source

If you are somewhat familiar with coding and git, I'd recommend this approach: You'll need to install [Rust](https://www.rust-lang.org/tools/install) and [git](https://git-scm.com/downloads) first. Then, you can clone this repository and build the program yourself:

```bash
git clone https://github.com/constantingoeldel/alphabeta-rs.git
cd alphabeta-rs
cargo install --path .
```

This will install the programs on your system. You can then ensure everything works by running:

```bash
alphabeta --help
# and
metaprofile --help
```

If you received an error message about libopenblas, you will need to run the code with cargo (I don't really understand this issue)

```bash
cargo run --release --bin alphabeta
# or
cargo run --release --bin metaprofile
```

### Updating

If you want to use a new version of the program, either download the new binaries from the same source or run:

```bash
cd alphabeta-rs
git pull
cargo install --path .
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
  -m, --methylome <METHYLOME>      Path to directory containing the methlyome files from which to extract the CG-sites
  -g, --genome <GENOME>            Path of the annotation file containing information about beginning and end of gbM-genes
  -w, --window-size <WINDOW_SIZE>  Size of the window in percent of the gbM-gene length or in basepair number if --absolute is supplied [default: 5]
  -s, --window-step <WINDOW_STEP>  Size of the step between the start of each window. Default value is window-size, so no overlapp happens
  -o, --output-dir <OUTPUT_DIR>    Path of the directory where extracted segments shall be stored
  -a, --absolute                   Use absolute length in base-pairs for window size instead of percentage of gene length
  -c, --cutoff <CUTOFF>            Number of basepairs to include upstream and downstream of gene [default: 2048]
  -i, --invert                     Invert strands, to switch from 5' to 3' and vice versa
      --db                         Use a Postgres database to do everything
  -e, --edges <EDGES>              Provide an edgefile
  -n, --nodes <NODES>              Provide a nodefile - paths will be updated to match the output directory
      --alphabeta                  Also run AlphaBeta on every window after extraction, results will be stored in the same directory as the segments
      --name <NAME>                Name of the run to be used when storing the result in Postgres [default: "Instant { tv_sec: 36502, tv_nsec: 792133216 }"]
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
