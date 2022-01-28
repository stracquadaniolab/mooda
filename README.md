# MOODA: Multi-Objective Optimization for DNA design and assembly

Current version: 0.10.1

![build](https://github.com/stracquadaniolab/mooda/workflows/release/badge.svg)
![platform](https://anaconda.org/stracquadaniolab/mooda/badges/platforms.svg)
![anaconda](https://anaconda.org/stracquadaniolab/mooda/badges/version.svg)

MOODA is a multi-objective optimisation algorithm for DNA sequence design and assembly.

It takes in input an annotated sequence in GenBank format, and optimize it with
respect to user-defined objectives.

Currently, some of the most common common operations in synthetic biology are
built-in, including:

- The `GCOptimizationOperator` introduces silent mutation in coding regions to
  obtain DNA constructs with a user-defined GC content.

- The `CodonUsageOperator` probabilistically recodes coding regions by
  probabilistically selecting the most frequent codon for an aminoacid in a host
  organism.

- The `BlockJoin` and `BlockSplit` operators allow the division of a sequence
  into fragments (or blocks). After the optimisation, each block is then adapted
  to the assembly method selected by the user. Currently, only the Gibson
  assembly is supported.

New operators, objective functions or assembly method can be added by extending
the `Operator`, `ObjectiveFunction` and `Assembly` classes.

## Installation

The easiest and fastest method to use `mooda` is using Docker:

```
    docker pull ghcr.io/stracquadaniolab/mooda
```

You can also install `mooda` using `conda`:

```
    $ conda install -c stracquadaniolab -c bioconda -c conda-forge mooda
```

or using `pip`:

```
    $ pip install mooda-dna
```

Please note, that `pip` will not install non Python requirements.

## Getting started

A typical `mooda` analysis consists of 3 steps:

1. Select a DNA sequence in Genbank format.

2. Write a MOODA configuration file. A `.yaml` file defining operators,
   objective functions, assemblies strategy and their parameters.

3. Run MOODA.

### Example: optimizing GC content, E. coli codon usage, number of fragments and the variance of their length

Create a test directory as follows:

```
    $ mkdir example-run
```

Move to your test directory as follows:

```
    $ cd example-run
```

Download test data from Github as follows:

```
    $ curl -LO https://github.com/stracquadaniolab/mooda/raw/master/examples/mooda-example1.tar.gz
```

Extract test data as follows:

```
    $ tar xvzf mooda-example1.tar.gz
```

Run `mooda` as follows:

```
    $ docker run -it --rm -v $PWD:$PWD -w $PWD ghcr.io/stracquadaniolab/mooda -i seq_5_5.gb  -c config.yaml -p 10 -it 20 -a 100 -mns 200 -mxs 2000 -bss 50 -js 40 -dir example-opt -gf True
```

Results will be available in the `example-opt` directory, where you will find:

- `Genbank` files of the Pareto optimal sequence.
- `FASTA` files with the fragments for Gibson assembly for each Pareto optimal
  sequence.
- `_logfile.yaml` file with information about the analysis.
- `_mooda_result.csv` file with objective function value information for each
  sequence.

#### Command line options

- **-i**: Input DNA sequence to process.

- **-c**: Configuration file to set operators, objective functions and their
  parameters.

- **-p**: Pool size. Number of candidate solutions sampled at each iteration.
  The pool size should increase with the length and complexity of the input
  sequence.

- **-it**: Number of iterations.  The number of iterations should increase with
  the length and complexity of the input sequence, although it will take longer
  to run.

- **-a**: Archive size. The number of non-dominated solutions to store at each
  iteration, which allows to use smaller pools for improved efficiency.

- **-mns**: Block minimum size.

- **-mxs**: Block maximum size.

- **-bss**: Sequence block step size, define the minimum variance between block
  size. Default: 50.

- **-js**: Sequence block assembly overlap size, define the amount of overlap
  between blocks. Default: 40.

- **-dir**: Output directory for MOODA results.

- **-gf**: Allow the writing of FASTA and GenBank files. Default: False.

## Authors

- Angelo Gaeta, a.gaeta@sms.ed.ac.uk
- Giovanni Stracquadanio, giovanni.stracquadanio@ed.ac.uk

## Citation

[Design and assembly of DNA molecules using multi-objective optimization](https://academic.oup.com/synbio/article-abstract/6/1/ysab026/6387748).
A Gaeta, V Zulkower, G Stracquadanio - Synthetic Biology, 2021

```
@article{10.1093/synbio/ysab026,
    author = {Gaeta, Angelo and Zulkower, Valentin and Stracquadanio, Giovanni},
    title = "{Design and assembly of DNA molecules using multi-objective optimization}",
    journal = {Synthetic Biology},
    volume = {6},
    number = {1},
    year = {2021},
    month = {10},
    issn = {2397-7000},
    doi = {10.1093/synbio/ysab026},
    url = {https://doi.org/10.1093/synbio/ysab026},
    note = {ysab026},
    eprint = {https://academic.oup.com/synbio/article-pdf/6/1/ysab026/40977182/ysab026.pdf},
}
```

## Issues

Please post an issue to report a bug or request new features.
