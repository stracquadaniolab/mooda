"""
    Class EvolutionaryAlgorithm

    Our workhorse optimiser for the multi-objective DNA assembly problem.

"""
import logging
import numpy as np
import datetime
from Bio import SeqIO
from mooda.individual import Individual
from mooda.population import Population
from mooda.archive import Archive
from mooda.operator import BlockInitialiserOperator
from mooda.utils import GenbankFile



class EvolutionaryAlgorithm:
    def __init__(self):
        self.gebank_record = None
        self.target_genetic_code = None
        self.population = Population()
        self.record = None
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

    def set_genetic_code(self, options):
        self.target_genetic_code = options.target_genetic_code

    def add_operator(self, op):
        self.operators.append(op)

    def add_objective_function(self, objf):
        self.objective_functions.append(objf)

    def add_assembly(self, assembly):
        self.assemblies = assembly

    def set_population_size(self, options):
        self.population.size = options.population_size

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
        self.set_genetic_code(options)
        self.record = GenbankFile(self.gebank_record)
        self.record.genbank_verification(options)
        self.set_population_size(options)
        self.set_max_iterations(options)
        self.set_genbank_fasta_output(options)
        self.set_archive_length(options)
        self.set_block_initialiser(options)

    def __initialize(self,options):
        for it in range(self.population.size):
            ind = Individual(self.gebank_record)
            ind.initialise(options)
            # evaluate objective functions for each individual
            self.block_initialiser.apply(ind)
            self.__eval_individual(ind)
            self.population.add_individual(ind)
        return self.population

    def run(self):
        self.time_start = datetime.datetime.now()
        self.__initialize(self.target_genetic_code)
        for it in range(self.max_iterations):
            self.current_iteration += 1
            #logging.info("itr=" + str(self.current_iteration))
            for i in range(len(self.population.individuals)):
                offspring = self.population.individuals[i].clone()
                op_curr = np.random.choice(self.operators)
                op_curr.apply(offspring)
                self.__eval_individual(offspring)
                self.population.individuals.append(offspring)
            self.population.compute_ranking()
            self.population.initialise_fronts()
            self.population.compute_crowding_distance(self.objective_functions)
            self.population.select_individuals()
            self.archive.add_to_archive(self.population)

            if len(self.archive.individuals) >= self.archive.length:
                self.archive.compute_ranking()
                self.archive.initialise_fronts()
                self.archive.select_only_rank1()
                self.archive.remove_duplicates()
                self.archive.random_deletion()

        self.assemblies[0].apply(self.population)
        self.__termination_process()
        self.time_end = datetime.datetime.now()
        duration = "duration= " + str(self.time_end - self.time_start)
        logging.info(duration)


    def __eval_individual(self, ind):
        ind.objectives = []
        # evaluate objfun for each individual
        for obj in self.objective_functions:
            ind.objectives.append(obj.eval(ind))


    def __termination_process(self):
        self.population.individuals = (
            self.population.individuals + self.archive.individuals
        )
        self.population.compute_ranking()
        self.population.initialise_fronts()
        self.population.compute_crowding_distance(self.objective_functions)
        self.population.select_only_rank1()
        self.population.remove_duplicates()
        self.population.individuals.sort(
            key=lambda i: i.crowding_distance, reverse=True
        )
