"""
    Class EvolutionaryAlgorithm

    Our workhorse optimiser for the multi-objective DNA assembly problem.

"""
import logging
import numpy as np
import datetime
import random
from Bio import SeqIO
from mooda.individual import Individual
from mooda.archive import Archive
from mooda.operator import BlockInitialiserOperator


class MonteCarloAlgorithm:
    def __init__(self):
        self.gebank_record = None
        self.individual = None
        self.operators = []
        self.objective_functions = []
        self.assemblies = None
        self.current_iteration = 0
        self.max_iterations = None
        self.yaml = None
        self.dataframe = None
        self.archive = Archive()
        self.block_initialiser = BlockInitialiserOperator()
        self.time_start = None
        self.time_end = None
        self.outputfasta_gebank = None

    # load a genbank file assuming one record per file.

    def set_yaml_file(self, config):
        self.yaml = config

    def load_sequence(self, options):
        self.gebank_record = SeqIO.read(open(options.input_file, "r"), "genbank")

    def add_operator(self, op):
        self.operators.append(op)

    def add_objective_function(self, objf):
        self.objective_functions.append(objf)

    def set_max_iterations(self, options):
        self.max_iterations = options.iterations

    def set_genbank_fasta_output(self, options):
        self.outputfasta_gebank = options.sequence_output

    def set_archive_length(self, options):
        self.archive.length = options.archive_size

    def set_block_initialiser(self, options):
        self.block_initialiser.step_size = options.block_step_size
        self.block_initialiser.max_block_size = options.block_max_size - options.block_junction_size
        self.block_initialiser.min_block_size = options.block_min_size

    def load_setting(self, config, options):
        self.set_yaml_file(config)
        self.load_sequence(options)
        self.set_max_iterations(options)
        self.set_genbank_fasta_output(options)
        self.set_archive_length(options)
        self.set_block_initialiser(options)

    def __initialize(self):
        self.individual = Individual(self.gebank_record)
        self.individual.initialise()
        # evaluate objective functions for each individual
        self.block_initialiser.apply(self.individual)
        self.__eval_individual(self.individual)

    def run(self):
        self.time_start = datetime.datetime.now()
        self.__initialize()
        for it in range(self.max_iterations):
            self.current_iteration += 1
            logging.info("itr=" + str(self.current_iteration))
            offspring = self.individual.clone()
            op_curr = np.random.choice(self.operators)
            op_curr.apply(offspring)
            self.__eval_individual(offspring)
            self.__duel(self.individual, offspring)
        self.__termination()
        self.time_end = datetime.datetime.now()
        duration = "duration= " + str(self.time_end - self.time_start)
        logging.info(duration)

    def __duel(self, parent, offspring):
        ind_domination_counter = 0
        offspring_domination_counter = 0
        for objective in range(0, (len(parent.objectives))):
            if parent.objectives[objective] < offspring.objectives[objective]:
                ind_domination_counter += 1
        for objective in range(0, (len(parent.objectives))):
            if offspring.objectives[objective] < parent.objectives[objective]:
                offspring_domination_counter += 1
        if offspring_domination_counter > ind_domination_counter:
            self.individual = offspring
            self.archive.individuals.append(offspring)
        elif offspring_domination_counter == ind_domination_counter:
            self.archive.individuals.append(offspring)
            self.individual = random.choice([offspring, parent])

    def __termination(self):
        self.archive.compute_ranking()
        self.archive.initialise_fronts()
        self.archive.compute_crowding_distance(self.objective_functions)
        self.archive.select_only_rank1()
        self.archive.remove_duplicates()
        self.assemblies[0].apply(self.archive)



    def __eval_individual(self, ind):
        ind.objectives = []
        # evaluate objfun for each individual
        for obj in self.objective_functions:
            ind.objectives.append(obj.eval(ind))
