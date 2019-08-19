# MOODA: Multi-Objective Optimization for DNA sequence Design and Assembly

Current version: 0.4.1-dev

![build](https://circleci.com/gh/stracquadaniolab/baghera/tree/master.svg?style=svg)
![platform](https://anaconda.org/stracquadaniolab/baghera/badges/platforms.svg)
![anaconda](https://anaconda.org/stracquadaniolab/baghera/badges/version.svg)


MOODA is a multi-objective optimisation algorithm for sequence Design and assembly. 
It takes as input an annotated sequence in GENBANK format, and optimize it throughout to defined operators. They to perform some of the most common operations in synthetic biology :

The **GC content operator** reduces the difference between the GC content of a sequence and the GC content set as the target. It introduces silent mutation inside CDSs, to increase or decrease the GC content.

The **Codon usage** operator allows the recoding of CDSs according to the specified codon distribution. At each iteration, a specified number of codons is replaced by synonymous

The **Block Join** and **Block split** operators allow the division of the sequence into blocks, given a minimum and maximum size. After the optimisation, each block is then adapted to the selected assembly method.
New operators, objective functions or assembly method can be integrated into the algorithm as python sub-classes.



## Installation

The easiest and fastest way to install TOOL using `conda`:

    $ conda install -c stracquadaniolab -c bioconda -c conda-forge mooda

Alternatively, you can install `MOODA` through `pip`:

    $ pip install mooda

Please note, that pip will not install non Python requirements.

## Getting started

A typical `MOODA` analysis consists of 3 steps:

1. Select a DNA sequence in Genbank format, an example file seq_5_5.gb is provided in test/mooda_test.zip

2. Write MOODA configuration file, a .yaml file defining operators, objective functions, assemblies strategy and their parameters, this is how a MOODA configuration file looks like:


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

an example of a MOODA configuration file: gc_codonusage_blockvariance_blocknumber.yaml is already available in tests/mooda_test.zip. To run MOODA with this configuration, a .yaml specifying the codon to use for each amino acid and their frequency must be specified in the configuration file, an example is provided in mooda_test.zip, where E. coli codon distribution is taken as an example.
 3. Run MODOA according to your parameters.

### Example

Running `MOODA` on this input as follows:


    $ mooda -ag mo -i seq_5_5.gb  -c gc_codonusage_blockvariance_blocknumber.yaml -p 10 -it 20 -a 100 -mns 200 -mxs 2000 -bss 50 -js 40 -dir mooda_results_dir -gf True 

**-ag** Algorithm to run can be either mo for Multi-Objective, either mc for Monte Carlo, mo is suggested for long sequences, Monte Carlo for small sequences and codon usage optimization. Default=mo.

**-i** Input DNA sequence to optimise.

**-c** Configuration file to set MOODA operators, objective function and their parameters.

**-p** Population size, set the number of clones of the input sequence, on each clone the it number of operators will be applied. The -p parameter should increase with the sequence size. It improves solution quality, however the computing time increase as well.

**-it** Number of iterations, set the number of operators to apply to each clone defined by -p. The -it parameter should increase with the sequence size. It improves solution quality more than -p parameter, however the computing time increase as well -a Archive size, amount of non-dominated solutions to store at each algorithm iteration, allow to use smaller values for **-p**.

**-mns** Sequence block minimum size.

**-mxs** Sequence block maximum size.

**-bss** Sequence block step size, define the minimum variance between block lengths. Default=50.

**-js** Sequence block assembly overlap size, define the amount of overlap between sequence blocks. Default=40.
-dir Output directory for MODOA results

**-gf** Allow the writing of FASTA and GENBANK files, related to MOODA solution if set as True. Default=False.

## Authors

- Angelo Gaeta, a.gaeta@sms.ed.ac.uk
- Giovanni Stracquadanio, giovanni.stracquadanio@ed.ac.uk

## Citation

Design and assembly of long DNA molecules using multi-objective optimisation..
Angelo Gaeta and  Giovanni Stracquadanio.
bioRxiv XX; doi: XX

## Issues

Please post an issue to report a bug or request new features.
