# MOODA: Multi-Objective Optimization for DNA sequence Design and Assembly

Current version: 0.7.3-dev

![build](http://drone.stracquadaniolab.org/api/badges/stracquadaniolab/mooda/status.svg)
![platform](https://anaconda.org/stracquadaniolab/mooda/badges/platforms.svg)
![anaconda](https://anaconda.org/stracquadaniolab/mooda/badges/version.svg)


MOODA is a multi-objective optimisation algorithm for sequence Design and Assembly.

It takes as input an annotated sequence in GenBank format, and optimize it with respect to user-specified objectives.

Currently, some of the most common common operations in synthetic biology are implemented:

- The **GC content operator** reduces the difference between the GC content of a sequence and the GC content set as the target. It introduces silent mutation inside CDSs, to increase or decrease the GC content.

- The **Codon usage** operator allows the recoding of CDSs according to the specified codon distribution. At each iteration, a specified number of codons is replaced by synonymous

- The **Block Join** and **Block split** operators allow the division of the sequence into blocks, given a minimum and maximum size. After the optimisation, each block is then adapted to the selected assembly method. Currently, only Gibson assembly is supported.

New operators, objective functions or assembly method can be integrated into the algorithm as python sub-classes.



## Installation

The easiest and fastest way to install `mooda` using `conda`:

    $ conda install -c stracquadaniolab -c bioconda -c conda-forge mooda

Alternatively, you can install `mooda` through `pip`:

    $ pip install mooda-web

Please note, that pip will not install non Python requirements.

## Getting started

A typical `mooda` analysis consists of 3 steps:

1. Select a DNA sequence in Genbank format.

2. Write a MOODA configuration file. A .yaml file defining operators, objective functions, assemblies strategy and their parameters, this is how a MOODA configuration file looks like:

```
    Algorithm :

        operators :
            mooda.operator.SplitBlockOperator :
                min_block_size : 200
                max_block_size : 2000
                step_size : 50

            mooda.operator.JoinBlockOperator :
                min_block_size : 200
                max_block_size : 2000
                junction_size : 40
                step_size : 50

            mooda.operator.GCOptimizationOperator :
                codon_GC_table: "e_coli_codon_usage.yaml"
                target_gc : 50
                step_size : 0.05


            mooda.operator.CodonUsageOperator :
                step_size : 0.05
                codon_usage_table : "e_coli_codon_usage.yaml"

        objective_functions :

                mooda.objective_function.GCContentObjective :
                    target_gc : 50
                    junction_size : 40

                mooda.objective_function.BlockVarianceObjective:
                    junction_size : 40

                mooda.objective_function.BlockNumberObjective:

                mooda.objective_function.CodonUsageObjective :
                    codon_usage_table:"e_coli_codon_usage.yaml"

        assemblies :
                mooda.assembly.Gibson:
                    junction_size : 40
```

3. Run MOODA.


### Example
Test data are provided in `test/mooda_test.zip`.

You can run `mooda` on the test data as follows:


    $ mooda -ag mo -i seq_5_5.gb  -c gc_codonusage_blockvariance_blocknumber.yaml -p 10 -it 20 -a 100 -mns 200 -mxs 2000 -bss 50 -js 40 -dir mooda_results_dir -gf True

**-ag** Algorithm to run can be either mo for Multi-Objective, either mc for Monte Carlo, mo is suggested for long sequences, Monte Carlo for small sequences and codon usage optimization. Default=mo.

**-i** Input DNA sequence to process.

**-c** Configuration file to set MOODA operators, objective functions and their parameters.

**-p** Pool size. The -p parameter should increase with the sequence size. It improves solution quality, however the computing time increase as well.

**-it** Number of iterations. The -it parameter should increase with the sequence size. It improves solution quality more than -p parameter, however the computing time increase as well 

**-a** Archive size, amount of non-dominated solutions to store at each algorithm iteration, allow to use smaller values for the pool size.

**-mns** Sequence block minimum size.

**-mxs** Sequence block maximum size.

**-bss** Sequence block step size, define the minimum variance between block lengths. Default: 50.

**-js** Sequence block assembly overlap size, define the amount of overlap between sequence blocks. Default: 40.

**-dir** Output directory for MOODA results.

**-gf** Allow the writing of FASTA and GenBank files, related to MOODA solution if set as True. Default=False.

## Authors

- Angelo Gaeta, a.gaeta@sms.ed.ac.uk
- Giovanni Stracquadanio, giovanni.stracquadanio@ed.ac.uk

## Citation

Design and assembly of DNA molecules using multi-objective optimisation.
Angelo Gaeta, Valentin Zulkower and Giovanni Stracquadanio.
bioRxiv XX; doi: XX

## Issues
Please post an issue to report a bug or request new features.
