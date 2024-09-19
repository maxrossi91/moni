[![Release](https://img.shields.io/github/release/maxrossi91/moni.svg)](https://github.com/maxrossi91/moni/releases)
[![Downloads](https://img.shields.io/github/downloads/maxrossi91/moni/total?logo=github)](https://github.com/maxrossi91/moni/releases/download/v0.2.2/moni_0.2.2_amd64.deb)

# MONI
```console
                           __  __  ____  _   _ _____
                          |  \/  |/ __ \| \ | |_   _|
                          | \  / | |  | |  \| | | |
                          | |\/| | |  | | . ` | | |
                          | |  | | |__| | |\  |_| |_
                          |_|  |_|\____/|_| \_|_____|
                                            ver 0.2.2
```
A Pangenomics Index for Finding MEMs.

MONI index uses the prefix-free parsing of the text [2][3] to build the Burrows-Wheeler Transform (BWT) of the reference genomes, the suffix array (SA) samples at the beginning and at the end of each run of the BWT, and the threshold positions of [1]. 


### Construction of the index:
```
usage: moni build [-h] -r REFERENCE [-w WSIZE] [-p MOD] [-t THREADS] [-k] [-v]
                  [-f] [--moni-ms] [--spumoni]
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        reference file name (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: same as reference)
  -w WSIZE, --wsize WSIZE
                        sliding window size (default: 10)
  -p MOD, --mod MOD     hash modulus (default: 100)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 0)
  -k                    keep temporary files (default: False)
  -v                    verbose (default: False)
  -f                    read fasta (default: False)
  -g GRAMMAR, --grammar GRAMMAR
                        select the grammar [plain, shaped] (default: plain)

```


### Computing the matching statistics with MONI:
```
usage: moni ms [-h] -i INDEX -p PATTERN [-o OUTPUT] [-t THREADS]
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index base name (default: None)
  -p PATTERN, --pattern PATTERN
                        the input query (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: .)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 1)
  -g GRAMMAR, --grammar GRAMMAR
                        select the grammar [plain, shaped] (default: plain)
```

### Computing the matching statistics with MONI:
```
usage: moni mems [-h] -i INDEX -p PATTERN [-o OUTPUT] [-e] [-s] [-t THREADS]
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index base name (default: None)
  -p PATTERN, --pattern PATTERN
                        the input query (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: .)
  -e, --extended-output
                        output MEM occurrence in the reference (default: False)
  -s, --sam-output
                        output MEM in a SAM formatted file. (default: False)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 1)
  -g GRAMMAR, --grammar GRAMMAR
                        select the grammar [plain, shaped] (default: plain)
```

### Computing the MEM extension with MONI and ksw2:
```
usage: moni extend [-h] -i INDEX -p PATTERN [-o OUTPUT] [-t THREADS] [-b BATCH] [-g GRAMMAR] [-L EXTL] [-A SMATCH] [-B SMISMATCH] [-O GAPO] [-E GAPE]

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index folder (default: None)
  -p PATTERN, --pattern PATTERN
                        the input query (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: .)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 1)
  -b BATCH, --batch BATCH
                        number of reads per thread batch (default: 100)
  -g GRAMMAR, --grammar GRAMMAR
                        select the grammar [plain, shaped] (default: plain)
  -L EXTL, --extl EXTL  length of reference substring for extension (default: 100)
  -A SMATCH, --smatch SMATCH
                        match score value (default: 2)
  -B SMISMATCH, --smismatch SMISMATCH
                        mismatch penalty value (default: 4)
  -O GAPO, --gapo GAPO  coma separated gap open penalty values (default: 4,13)
  -E GAPE, --gape GAPE  coma separated gap extension penalty values (default: 2,1)
```

# Example
### Install prerequisite packages

```console
apt-get update
apt-get install -y build-essential cmake git python3 zlib1g-dev
```

### Download

```console
git clone https://github.com/maxrossi91/moni
```

### Compile

```console
mkdir build
cd build
cmake ..
make
```
### Install

```console
make install
```

This command will install the binaries to the default install location (e.g., `/usr/local/bin` for Ubuntu users). If the user wants the binary in some other custom location, this can be done using `cmake -DCMAKE_INSTALL_PERFIX=<dest> ..` instead of `cmake ..` in the compile sequence of commands, where `<dest>` is the preferred destination directory.

### Run

##### Build the index for `SARS-CoV2.1k.fa.gz` in the `data/SARS-CoV2` folder
```console
moni build -r data/SARS-CoV2/SARS-CoV2.1k.fa.gz -o sars-cov2 -f
```
It produces three files `sars-cov2.plain.slp`, `sars-cov2.thrbv.ms`, and `sars-cov2.idx` in the current folder which contain the grammar, the rlbwt and the thresholds, and the starting position and name of each fasta sequence in the reference file respectively.

##### Compute the matching statistics of `reads.fastq.gz ` against `SARS-CoV2.1k.fa.gz` in the `data/SARS-CoV2` folder
```console
moni ms -i sars-cov2 -p data/SARS-CoV2/reads.fastq.gz -o reads
```
It produces two output files `reads.lengths` and `reads.pointers` in the current folder which store the lengths and the positions of the matching statistics of the reads against the reference in a fasta-like format.  

##### Compute the MEMs of `reads.fastq.gz ` against `SARS-CoV2.1k.fa.gz` in the `data/SARS-CoV2` folder
```console
moni mems -i sars-cov2 -p data/SARS-CoV2/reads.fastq.gz -o reads
```
It produces one output file `reads.mems` in the current folder which store the MEMs reposted as pairs of position and lengths in a fasta-like format.  

##### Compute the MEM extension of `reads.fastq.gz ` against `SARS-CoV2.1k.fa.gz` in the `data/SARS-CoV2` folder
```console
moni extend -i sars-cov2 -p data/SARS-CoV2/reads.fastq.gz -o reads
```
It produces one output file `reads.sam` in the current folder which stores the information of the MEM extensions in SAM format.  
# External resources

* [Big-BWT](https://github.com/alshai/Big-BWT.git)
    * [gSACA-K](https://github.com/felipelouza/gsa-is.git)
    * [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [Divsufsort](https://github.com/simongog/libdivsufsort.git)
* [klib](https://github.com/attractivechaos/klib)
* [ksw2](https://github.com/lh3/ksw2)
* [r-index](https://github.com/maxrossi91/r-index.git)
* [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds.git)
* [bigrepair](https://gitlab.com/manzai/bigrepair.git)
* [shaped_slp](https://github.com/koeppl/ShapedSlp.git)
<!-- * [Google Benchmark](https://github.com/google/benchmark.git)
    * [Google Test](https://github.com/google/googletest) -->

# Citation 

Please, if you use this tool in an academic setting cite the following papers:

    @article{RossiOLGB21,
    author      = { Massimiliano Rossi and 
                    Marco Oliva and
                    Ben Langmead and
                    Travis Gagie and
                    Christina Boucher},
    title       = {MONI: A Pangenomics Index for Finding Maximal Exact Matches},
    booktitle   = {Research in Computational Molecular Biology - 25th Annual 
                    International Conference, {RECOMB} 2021, Padova, Italy},
    journal     = {Journal of Computational Biology},
    volume      = {29},
    number      = {2},
    pages       = {169--187},
    year        = {2022},
    publisher   = {Mary Ann Liebert, Inc., publishers 140 Huguenot Street, 3rd Floor New~…}
    }


# Authors

### Theoretical results:

* Christina Boucher
* Travis Gagie
* Ben Langmead
* Massimiliano Rossi

### Implementation:

* [Massimiliano Rossi](https://github.com/maxrossi91)

### Experiments

* [Marco Oliva](https://github.com/marco-oliva)
* [Massimiliano Rossi](https://github.com/maxrossi91)

# Why "MONI"?

**Moni** is the Finnish word for *multi*.

# References

[1] Hideo Bannai, Travis Gagie, and Tomohiro I, *"Refining ther-index"*, Theoretical Computer Science, 812 (2020), pp. 96–108

[2] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini, *"Prefix-Free Parsing for Building Big BWTs"*, In Proc. of the 18th International Workshop on Algorithms in Bioinformatics (WABI 2018).

[3] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.