"""
    Class ObjectiveFunction

    desc
"""
import numpy as np
from Bio.SeqUtils import GC
from mooda.config import CodonTable
from mooda.config import RepetitionTable


class ObjectiveFunction:
    def __init__(self, yaml):
        self.yaml = yaml

    def eval(self, ind):
        pass


"""
    Class GCContentObj

    desc
"""


class GCContentObjective(ObjectiveFunction):
    def __set_target_gc(self):
        self.target_gc = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.GCContentObjective"
        ]["target_gc"]

    def initialise(self):
        self.__set_target_gc()

    def eval(self, ind):
        cum_sum = 0.0

        for pos in ind.blocks:
            gc_val = GC(ind.sequence[pos[0]: pos[1]])
            cum_sum += np.abs(gc_val - self.target_gc)
        return cum_sum / float(len(ind.blocks))


"""
    Class BlockVarianceObjective

    desc
"""


class BlockVarianceObjective(ObjectiveFunction):
    def initialise(self):
        pass

    def eval(self, ind):
        block_length_list = []
        for bb in ind.blocks:
            blocksize = bb[1] - bb[0]
            block_length_list.append(blocksize)
        block_variance = np.var(block_length_list, dtype=np.float64)
        return block_variance


class BlockNumberObjective(ObjectiveFunction):
    def initialise(self):
        pass

    def eval(self, ind):
        block_counter = len(ind.blocks)
        return block_counter


"""
    Class BasePairCostObjective

    desc
"""


class BasePairCostObjective(ObjectiveFunction):
    def set_basepair_cost(self):
        self.basepair_cost = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.BasePairCostObjective"
        ]["basepair_cost"]

    def set_block_cost(self):
        self.block_cost = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.BasePairCostObjective"
        ]["block_cost"]

    def initialise(self):
        self.set_basepair_cost()
        self.set_block_cost()

    def eval(self, ind):
        # turn blocks attribute into a list
        # for each block in the list
        cum_sum = 0.0

        for bb in ind.blocks:
            # evaluate the length of the
            blocksize = bb[1] - bb[0]
            cum_sum += self.block_cost + (self.basepair_cost * blocksize)
        return cum_sum


"""
    Class CodonUsage

    desc
"""


class CodonUsageObjective(ObjectiveFunction):
    def set_codon_usage_table(self):
        self.codon_usage_table = CodonTable()
        self.codon_usage_table.codons = self.codon_usage_table.load_config(
            self.yaml["Algorithm"]["objective_functions"][
                "mooda.objective_function.CodonUsageObjective"
            ]["codon_usage_table"]
        )

    # command to return the highest frequency from
    def __get_codon_highest_frequency(self, aa):

        codon_dict = self.codon_usage_table.codons[aa]
        codon_highest_frequency = max(codon_dict, key=codon_dict.get)
        highest_frequency = codon_dict[codon_highest_frequency]
        return highest_frequency

    def __get_codon_frequency(self, codon):
        aa = codon.translate()
        codon_freq = self.codon_usage_table.codons[aa][codon]
        return codon_freq

    def initialise(self):
        self.set_codon_usage_table()

    def eval(self, ind):
        cum_sum = 0
        for index in ind.cds_indexes_list:
            cds = ind.sequence[index.location.start:index.location.end]
            if index.location.strand != 1:
                cds = cds[::-1]
            cds_codons = [cds[i: i + 3] for i in range(0, len(cds), 3)]
            for codon in cds_codons:
                aa = codon.translate()
                codon_frequence = self.__get_codon_frequency(codon)
                codon_frequence_target = self.__get_codon_highest_frequency(aa)
                cum_sum += abs(codon_frequence - codon_frequence_target)
        return cum_sum


"""
    Class Repetition

    desc
"""


class RepetitionObjective(ObjectiveFunction):
    def initialise(self):
        self.repetition_table = RepetitionTable()
        self.repetition_table.motives_dict = self.repetition_table.load_config(
            self.yaml["Algorithm"]["objective_functions"][
                "mooda.objective_function.RepetitionObjective"
            ]["repetition_table"]
        )
        self.repetition_table.get_motives_list()

    def eval(self, ind):
        # objective f
        # unction value
        cum_sum = 0.0
        # for each CDS in
        for index in ind.cds_indexes_list:
            # sequence
            cds = ind.sequence[index.location.start:index.location.end]
            if index.location.strand != 1:
                cds = cds[::-1]
            for motive in self.repetition_table.motives:
                cum_sum += cds.count(motive)
        return cum_sum
