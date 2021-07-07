# MONI
```console
                           __  __  ____  _   _ _____
                          |  \/  |/ __ \| \ | |_   _|
                          | \  / | |  | |  \| | | |
                          | |\/| | |  | | . ` | | |
                          | |  | | |__| | |\  |_| |_
                          |_|  |_|\____/|_| \_|_____|
                                            ver 0.1.0
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

### Computing the MEM extension with MONI and ksw2:
```
usage: moni extend [-h] -i INDEX -p PATTERN [-o OUTPUT] [-t THREADS]
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index base name (default: None)
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
```

# Example
### Download

```console
git clone https://github.com/maxrossi91/moni
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

### Run

##### Build the index for `yeast.fasta` in the `data` folder
```console
./moni build -r ../data/yeast.fasta -f -t 4
```
It produces two files `yeast.fasta.plain.slp` and `yeast.fasta.thrbv.ms` in the `data` folder which contain the grammar and the rlbwt and the thresholds respectively.

##### Compute the matching statistics of `reads.fastq` against `yeast.fasta` in the `data` folder
```console
./moni ms -i ../data/yeast.fasta -p ../data/reads.fastq -t 4
```
It produces two output files `reads.fastq_yeast.fasta.lengths` and `reads.fastq_yeast.fasta.pointers` in the `data` folder which store the lengths and the positions of the matching statistics of the reads against the reference in a fasta-like format.  

##### Compute the MEM extension of `reads.fastq` against `yeast.fasta` in the `data` folder
```console
./moni extend -i ../data/yeast.fasta -p ../data/reads.fastq -t 4
```
It produces one output file `reads.fastq_yeast.fasta_25.sam` in the `data` folder which stores the information of the MEM extensions in sam format.  
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

    @inproceedings{RossiOLGB21,
    author      = { Massimiliano Rossi and 
                    Marco Oliva and
                    Ben Langmead and
                    Travis Gagie and
                    Christina Boucher},
    title       = {MONI: A Pangenomics Index for Finding MEMs},
    booktitle   = {Research in Computational Molecular Biology - 25th Annual 
                    International Conference, {RECOMB} 2021, Padova, Italy},
    volume      = {},
    series      = {Lecture Notes in Computer Science},
    pages       = {},
    publisher   = {Springer},
    year        = {2021}
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