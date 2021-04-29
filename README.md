# MOODA: Multi-Objective Optimization for DNA design and assembly

Current version: 0.9.5

![build](https://github.com/stracquadaniolab/mooda/workflows/release/badge.svg)
![platform](https://anaconda.org/stracquadaniolab/mooda/badges/platforms.svg)
![anaconda](https://anaconda.org/stracquadaniolab/mooda/badges/version.svg)

MOODA is a multi-objective optimisation algorithm for DNA sequence design and assembly.

It takes in input an annotated sequence in GenBank format, and optimize it with
respect to user-specified objectives.

Currently, some of the most common common operations in synthetic biology are implemented:

- The `GC content operator` reduces the difference between the GC content of a
  sequence and the GC content set as the target. It introduces silent mutation
  inside CDSs, to increase or decrease the GC content.

- The `Codon usage` operator allows the recoding of CDSs according to the
  specified codon distribution. At each iteration, a specified number of codons
  is replaced by synonymous

- The `BlockJoin` and `BlockSplit` operators allow the division of the
  sequence into blocks, given a minimum and maximum size. After the
  optimisation, each block is then adapted to the selected assembly method.
  Currently, only Gibson assembly is supported.

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
    $ docker run -it --rm -v $PWD:$PWD -w $PWD ghcr.io/stracquadaniolab/mooda -ag mo -i seq_5_5.gb  -c config.yaml -p 10 -it 20 -a 100 -mns 200 -mxs 2000 -bss 50 -js 40 -dir example-opt -gf True
```

Results will be available in the `example-opt` directory, where you will find:

- `Genbank` files of the Pareto optimal sequence.
- `FASTA` files with the fragments for Gibson assembly for each Pareto optimal
  sequence.
- `_logfile.yaml` file with information about the analysis.
- `_mooda_result.csv` file with objective function value information for each
  sequence.

#### Command line options

- **-ag**: Algorithm to run can be either mo for Multi-Objective, either mc for Monte Carlo, mo is suggested for long sequences, Monte Carlo for small sequences and codon usage optimization. Default: mo.

- **-i**: Input DNA sequence to process.

- **-c**: Configuration file to set MOODA operators, objective functions and
  their parameters.

- **-p**: Pool size. The -p parameter should increase with the sequence size. It
  improves solution quality, however the computing time increase as well.

- **-it**: Number of iterations. The -it parameter should increase with the
  sequence size. It improves solution quality more than -p parameter, however
  the computing time increase as well.

- **-a**: Archive size, amount of non-dominated solutions to store at each
  algorithm iteration, allow to use smaller values for the pool size.

- **-mns**: Sequence block minimum size.

- **-mxs**: Sequence block maximum size.

- **-bss**: Sequence block step size, define the minimum variance between block
  lengths. Default: 50.

- **-js**: Sequence block assembly overlap size, define the amount of overlap
  between sequence blocks. Default: 40.

- **-dir**: Output directory for MOODA results.

- **-gf**: Allow the writing of FASTA and GenBank files, related to MOODA
  solution if set as True. Default: False.

## Authors

- Angelo Gaeta, a.gaeta@sms.ed.ac.uk
- Giovanni Stracquadanio, giovanni.stracquadanio@ed.ac.uk

## Citation

Design and assembly of DNA molecules using multi-objective optimisation.
Angelo Gaeta, Valentin Zulkower and Giovanni Stracquadanio.
bioRxiv. https://www.biorxiv.org/content/10.1101/761320v1

## Issues

Please post an issue to report a bug or request new features.
