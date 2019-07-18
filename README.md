# MOODA: Multi-Objective Optimization for DNA sequence Design and Assembly

Current version: 0.3.1-dev

![build](https://circleci.com/gh/stracquadaniolab/baghera/tree/master.svg?style=svg)
![platform](https://anaconda.org/stracquadaniolab/baghera/badges/platforms.svg)
![anaconda](https://anaconda.org/stracquadaniolab/baghera/badges/version.svg)

MOODA is a multi-objective optimisation algorithm for sequence Design and assembly. 
It takes as input an annoted sequence in GENBANK format, and optimize it throughout to defined operators. They to perform some of the most common operation in synthetic biology : 

The **GC content** operator reduces the difference between the GC content of a sequence and the GC content set as target. 
It introduce silent mutation inside CDSs, to incresease or decrease the GC content.

The **Codon usage** operator allow the recoding of CDSs according to the specified codon distribtuion. At each iteration a specificed number of codons is replaced by synonimous

The **Block Join** and **Block split** operators allow the divsion of the sequence into blocks, given a minimum and maximum size. After the optimisation each block is then adapted to the selected assembly method.

New operators, objective functions or assembly method can be integrated to the algortihm as python sub-classes.


## Installation

The easiest and fastest way to install TOOL using `conda`:

    $ conda install -c stracquadaniolab -c bioconda -c conda-forge mooda

Alternatively you can install `mooda` through `pip`:

    $ pip install mooda

Please note, that `pip` will not install non Python requirements.

## Getting started

A typical `mooda` analysis consists of 3 steps:

1. Select a DNA sequence in Genbank format.
2. Write MOODA configuration file, a .yaml file defining operators, objective functions and assemblies strategy e.g.:


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
                                                    codon_GC_table: "{your/path/codon_usage.yaml}"
                                                    target_gc : 50
                                                    step_size : 0.05


                                    mooda.operator.CodonUsageOperator :
                                                 step_size : 0.05
                                                 codon_usage_table : "{your/path/codon_usage.yaml}"

                            objective_functions :

                                    mooda.objective_function.GCContentObjective :
                                                                        target_gc : 50
                                                                        junction_size : 40
                                          
                                    mooda.objective_function.BlockVarianceObjective:
                                                                        junction_size : 40
                                                                        
                                    mooda.objective_function.BlockNumberObjective:

                                    mooda.objective_function.CodonUsageObjective :
                                                                        codon_usage_table: "{your/path/codon_usage.yaml}"
                            assemblies :
                                    mooda.assembly.Gibson:
                                                            junction_size : 40


3. Run mooda according to your parameters .

### Example

Running `mooda` on this input as follows:

    $ mooda mooda -ag mo -i your/path/genbank_sequence.gb -c your/path/mooda_config.yaml -p 10 -it 20 -a 100  -mns 200 -mxs        2000 -bss 50  -js 40 -dir your/path/mooda_results_dir -gf True 

**-ag** Algorithm to run, can be either mo for Multi-Objective, either mc for Monte Carlo, mo is suggested for long sequences,
Monte Carlo for small sequences and for codon usage optimization. Default=mo

**-i** Input DNA sequence to optimise

**-c** Configuration file to set MOODA operators, objective function and their parameters

**-p** Population size, set the number of clones of the input sequence, on each clone the **it** number of operators will be applied. The **-p** parameter should increase with the sequence size. It improves solution quality, however the computing time increase as well 

**-it** Number of iterations, set the number of operators to apply to each clone defined by **-p**.The **-it** parameter should increase with the sequence size. It improves solution quality more than **-p** parameter, however the computing time increase as well 
**-a** Archive size, amount of non-dominated solutions to store at each algorithm iteration, allow to use smaller values for **-p**

**-mns** Sequence block minimum size

**-mxs** Sequence block maximum size

**-bss** Sequence block step size, define the minimum variance between block lengths. Default=50

**-js** Sequence block assembly overlap size, define the amount of overlap between sequence blocks. Default=40

**-dir** Output directory for MODOA results

**-gf** Allow the writing of FASTA and GENBANK files,related to MOODA solution if set as True. Default=False


## Documentation
The official documentation for mooda can be found on [readthedocs](https://mooda.readthedocs.io/).

## Authors

- Angelo Gaeta, a.gaeta@sms.ed.ac.uk
- Giovanni Stracquadanio, giovanni.stracquadanio@ed.ac.uk

## Citation

Design and assembly of long DNA molecules using multi-objective optimisation..
Angelo Gaeta and  Giovanni Stracquadanio.
bioRxiv XX; doi: XX

## Issues

Please post an issue to report a bug or request new features.
