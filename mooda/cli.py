import logging
import argparse
import importlib
import sys

from mooda.objective_function import ObjectiveFunction
from mooda.operator import Operator
from mooda.assembly import Assembly
from mooda.io import MoodaFile
from mooda.evolutionary_algorithm import EvolutionaryAlgorithm
from mooda.montecarlo_algorithm import MonteCarloAlgorithm
from mooda.config import YamlConfig


def load_external_classes(config, section, class_type):
    cls_list = []

    for op, conf in config["Algorithm"][section].items():
        fields = op.split(".")

        # determining the name of the module and the class name
        # GS: You're trying to import the whole class name
        # GS: for example: to load splitblock operator, you have to pass something like
        # GS: mooda.Operator.SplitBlockOperator
        # GS: once the package is in your path, that's all you need.

        mod_name, cls_name = ".".join(fields[:-1]), fields[-1]

        logging.debug("Loading {} from {}".format(cls_name, mod_name))

        mod = importlib.import_module(mod_name)

        # tries to load and create an instance of the class
        try:
            cls = getattr(mod, cls_name)(config)
        except AttributeError:
            logging.error("Can not find class {}".format(cls))
            sys.exit(-1)

        # check if the class is of type Plugin
        if isinstance(cls, class_type):
            # exectute something
            cls.initialise()
            cls_list.append(cls)
        else:
            logging.error("{} is not a valid {}.".format(op, class_type))
            sys.exit(-1)
    return cls_list


##############################################################
###### CMD line options parser
##############################################################
def parse_command_line_options():
    parser = argparse.ArgumentParser(
        description="Multi Objective Optimisation for DNA design and Assembly."
    )
    parser.add_argument(
        "--input-file",
        "-i",
        type=str,
        required=True,
        help="A GenBank file containing the construct to redesign for assembly.",
    )

    parser.add_argument(
        "--algorithm",
        "-ag",
        type=str,
        choices=['mo', 'mc'],
        required=False,
        default='mo',
        help="Algorithm to use caa be either mo= Multi-Objective or mc= MonteCarlo.",
    )

    parser.add_argument(
        "--yaml-config",
        "-c",
        type=str,
        required=True,
        help="A .yaml configuration file describing all the parameters",
    )
    parser.add_argument(
        "--iterations",
        "-it",
        type=int,
        required=True,
        help="Algorithm number of iterations ",
    )
    parser.add_argument(
        "--population-size",
        "-p",
        type=int,
        required=False,
        help="number of individuals per iteration",
    )
    parser.add_argument(
        "--sequence_output",
        "-gf",
        type=bool,
        required=False,
        default=False,
        help="if true generates genbank and fasta files related to non"
             "dominated solutions",
    )
    parser.add_argument(
        "--archive-size",
        "-a",
        type=int,
        required=False,
        default=0,
        help="number of non dominated solutions stored at each iteration",
    )
    parser.add_argument(
        "--block-min-size", "-mns", type=int, required=True, help="block minimum size"
    )
    parser.add_argument(
        "--block-max-size", "-mxs", type=int, required=True, help="block maximum size"
    )
    parser.add_argument(
        "--block-step-size", "-bss", type=int, required=False, default=50, help="block maximum size"
    )
    parser.add_argument(
        "--block-junction-size", "-js", type=int, required=False, default=40, help="block junction size"
    )
    parser.add_argument(
        "--target-genetic-code", "-gc",
        choices=[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
        required=False, default=1, help=
        """
        codon usage table
        
         1 = Standard
         2 = Vertebrate Mitochondrial
         3 = Yeast Mitochondrial
         4 = Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma
         5 = Invertebrate Mitochondrial
         6 = Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
         9 = Echinoderm Mitochondrial; Flatworm Mitochondrial
         10 = Euplotid Nuclear
         11 = Bacterial, Archaeal and Plant Plastid
         12 = Alternative Yeast Nuclear
         13 = Ascidian Mitochondrial
         14 = Alternative Flatworm Mitochondrial
         15 =  Blepharisma Macronuclear 
         16 = Chlorophycean Mitochondrial
         21 = Trematode Mitochondrial
         22 = Scenedesmus obliquus Mitochondrial
         23 = Thraustochytrium Mitochondrial 
         24 = Pterobranchia Mitochondrial 
         25 = Candidate Division SR1 and Gracilibacteri
         26 = Pachysolen tannophilus Nuclear 
         27 = Karyorelict Nuclear
         28 = Condylostoma Nuclear
         29 = Mesodinium Nuclear 
         30 = Peritrich Nuclear
         31 = Blastocrithidia Nuclear
         32 = Balanophoraceae Plastid
        """
    )
    parser.add_argument(
        "--output-dir", "-dir", type=str, required=True, help="output directory"
    )

    options = parser.parse_args()
    return options


def main():
    # setting the logger for the whole section
    logging.basicConfig(level=logging.INFO)

    # parsing options
    options = parse_command_line_options()
    config = YamlConfig()
    config = config.load_config(options.yaml_config)

    # selecting the algorithm
    ag = None

    if options.algorithm == "mo":
        ag = EvolutionaryAlgorithm()
    elif options.algorithm == "mc":
        ag = MonteCarloAlgorithm()

    # loading settings from yaml and command line
    ag.load_setting(config, options)

    # setting operators
    ag.operators = load_external_classes(config, "operators", Operator)

    # setting objective functions
    ag.objective_functions = load_external_classes(
        config, "objective_functions", ObjectiveFunction
    )
    # setting assembly
    ag.assemblies = load_external_classes(
        config, "assemblies", Assembly
    )

    # running the algorithm
    ag.run()

    # saving files
    mf = MoodaFile(config)
    # GS: you don't need to pass config yaml again, it's already in the MoodaFile object
    mf.save(options, ag)


##############################################################
###### MAIN method
##############################################################
if __name__ == "__main__":
    main()
