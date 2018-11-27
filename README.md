# epic2

epic2 is an ultraperformant reimplementation of SICER. It focuses on speed, low memory overhead and ease of use.

#### Under active development

We are extremely responsive to bugs, issues and installation problems. We are
proud to say that epic was 20% more downloaded on PyPI than MACS2 even, and we
believe that was (among other things) because it was easy to install and use. We
wish to make epic2 similarly easy to install and use, so please report any
issues you have.

#### Features

* easy to install and use
* reads sam, single-end bam, bed and bedpe (.gz)
* extremely fast
* very low memory requirements
* works both with and without input
* metadata for ~80 UCSC genomes built in
* easily use custom genomes and assemblies with --chromsizes and --effective-genome-fraction args

#### Quick Start

```
git clone git@github.com:endrebak/epic2.git
cd epic2
python setup.py install
epic2 -t examples/test.bed -c examples/control.bed > results.txt
```

#### Install

The installation currently requires bioconda htslib to be installed for setup.py
to find the appropriate headers. I will update the install script with more ways
to include the headers.

```
pip install epic2
```

#### CLI

```
usage: epic2 [-h] --treatment TREATMENT [TREATMENT ...]
             [--control CONTROL [CONTROL ...]] [--genome GENOME]
             [--keep-duplicates] [--bin-size BIN_SIZE]
             [--gaps-allowed GAPS_ALLOWED] [--fragment-size FRAGMENT_SIZE]
             [--false-discovery-rate-cutoff FALSE_DISCOVERY_RATE_CUTOFF]
             [--effective-genome-fraction EFFECTIVE_GENOME_FRACTION]
             [--chromsizes CHROMSIZES] [--verbose]

epic2. (Visit github.com/endrebak/epic2 for examples and help.)

optional arguments:
  -h, --help            show this help message and exit
  --treatment TREATMENT [TREATMENT ...], -t TREATMENT [TREATMENT ...]
                        Treatment (pull-down) file(s) in one of these formats:
                        bed, bedpe, bed.gz, bedpe.gz or (single-end) bam, sam.
                        Mixing file formats is allowed.
  --control CONTROL [CONTROL ...], -c CONTROL [CONTROL ...]
                        Control (input) file(s) in one of these formats: bed,
                        bedpe, bed.gz, bedpe.gz or (single-end) bam, sam.
                        Mixing file formats is allowed.
  --genome GENOME, -gn GENOME
                        Which genome to analyze. Default: hg19. If
                        --chromsizes and --egf flag is given, --genome is not
                        required.
  --keep-duplicates, -kd
                        Keep reads mapping to the same position on the same
                        strand within a library. Default: False.
  --bin-size BIN_SIZE, -bin BIN_SIZE
                        Size of the windows to scan the genome. BIN-SIZE is
                        the smallest possible island. Default 200.
  --gaps-allowed GAPS_ALLOWED, -g GAPS_ALLOWED
                        This number is multiplied by the window size to
                        determine the number of gaps (ineligible windows)
                        allowed between two eligible windows. Must be an
                        integer. Default: 3.
  --fragment-size FRAGMENT_SIZE, -fs FRAGMENT_SIZE
                        (Single end reads only) Size of the sequenced
                        fragment. Each read is extended half the fragment size
                        from the 5' end. Default 150 (i.e. extend by 75).
  --false-discovery-rate-cutoff FALSE_DISCOVERY_RATE_CUTOFF, -fdr FALSE_DISCOVERY_RATE_CUTOFF
                        Remove all islands with an FDR below cutoff. Default
                        0.05.
  --effective-genome-fraction EFFECTIVE_GENOME_FRACTION, -egf EFFECTIVE_GENOME_FRACTION
                        Use a different effective genome fraction than the one
                        included in epic2. The default value depends on the
                        genome and readlength, but is a number between 0 and
                        1.
  --chromsizes CHROMSIZES, -cs CHROMSIZES
                        Set the chromosome lengths yourself in a file with two
                        columns: chromosome names and sizes. Useful to analyze
                        custom genomes, assemblies or simulated data. Only
                        chromosomes included in the file will be analyzed.
  --quiet, -q           Do not write output messages to stderr.
```

#### Performance

<img src="graphs/speed_epic2_vs_SICER_no_bigwig.png" />
