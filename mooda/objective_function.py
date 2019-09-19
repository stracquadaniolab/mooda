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

    def __set_junction_size(self):
        self.junction_size = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.GCContentObjective"
        ]["junction_size"]

    def initialise(self):
        self.__set_target_gc()
        self.__set_junction_size()

    def eval(self, ind):
        cum_sum = 0.0
        for pos in ind.blocks[:-1]:
            gc_val = GC(ind.sequence[pos[0]: (pos[1] + self.junction_size)])
            cum_sum += np.abs(gc_val - self.target_gc)

        pos = ind.blocks[-1]
        gc_val = GC(ind.sequence[pos[0]: pos[1]])
        cum_sum += np.abs(gc_val - self.target_gc)
        return cum_sum / float(len(ind.blocks))

    def __repr__(self):
        return 'GC content'


"""
    Class BlockVarianceObjective

    desc
"""


class BlockVarianceObjective(ObjectiveFunction):

    def __set_junction_size(self):
        self.junction_size = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.BlockVarianceObjective"
        ]["junction_size"]

    def initialise(self):
        self.__set_junction_size()

    def eval(self, ind):
        block_length_list = []

        for bb in ind.blocks[:-1]:
            blocksize = (bb[1] + self.junction_size) - bb[0]
            block_length_list.append(blocksize)

        bb = ind.blocks[-1]
        blocksize = bb[1] - bb[0]
        block_length_list.append(blocksize)

        block_variance = np.var(block_length_list, dtype=np.float64)
        return block_variance

    def __repr__(self):
        return 'Block variance'


class BlockNumberObjective(ObjectiveFunction):
    def initialise(self):
        pass

    def eval(self, ind):
        block_counter = len(ind.blocks)
        return block_counter

    def __repr__(self):
        return 'Block number'


"""
    Class BasePairCostObjective

    desc
"""


class BasePairCostObjective(ObjectiveFunction):

    def __set_junction_size(self):
        self.junction_size = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.BasePairCostObjective"
        ]["junction_size"]


    def __set_basepair_cost(self):
        self.basepair_cost = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.BasePairCostObjective"
        ]["basepair_cost"]

    def __set_block_cost(self):
        self.block_cost = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.BasePairCostObjective"
        ]["block_cost"]

    def initialise(self):
        self.__set_junction_size()
        self.__set_basepair_cost()
        self.__set_block_cost()

    def eval(self, ind):
        # turn blocks attribute into a list
        # for each block in the list
        cum_sum = 0.0

        for bb in ind.blocks[:-1]:
            blocksize = (bb[1] + self.junction_size) - bb[0]
            cum_sum += self.block_cost + blocksize * self.basepair_cost

        bb = ind.blocks[-1]
        blocksize = bb[1] - bb[0]
        cum_sum += self.block_cost + blocksize * self.basepair_cost
        return cum_sum

    def __repr__(self):
        return 'Cost'

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

    def __get_codon_frequency(self, codon, codon_table):
        aa = codon.translate(codon_table)
        codon_freq = self.codon_usage_table.codons[aa][codon]
        return codon_freq

    def initialise(self):
        self.set_codon_usage_table()

    def eval(self, ind):
        cum_sum = 0
        for cds in ind.cds_list:
            cds_codon_table =cds.translation_table_target
            cds_seq = ind.sequence[cds.pt.location.start:cds.pt.location.end]
            if cds.pt.location.strand != 1:
                cds_seq = cds_seq.complement()
                cds_seq = cds_seq[::-1]
            cds_codons = [cds_seq[i: i + 3] for i in range(0, len(cds_seq), 3)]
            for codon in cds_codons:
                aa = codon.translate()
                if len(codon) == 3:
                    codon_frequence = self.__get_codon_frequency(codon,cds_codon_table)
                    codon_frequence_target = self.__get_codon_highest_frequency(aa)
                    cum_sum += abs(codon_frequence - codon_frequence_target)
        return cum_sum

    def __repr__(self):
        return 'Codon usage'



"""
    Class Repetition

    desc
"""


class MotifObjective(ObjectiveFunction):
    def initialise(self):
        self.repetition_table = RepetitionTable()
        self.repetition_table.motives_dict = self.repetition_table.load_config(
            self.yaml["Algorithm"]["objective_functions"][
                "mooda.objective_function.MotifObjective"
            ]["motif_table"]
        )
        self.repetition_table.get_motives_list()
        self.junction_size = self.yaml["Algorithm"]["objective_functions"][
            "mooda.objective_function.MotifObjective"
        ]["junction_size"]

    def eval(self, ind):
        # objective Function value
        cum_sum = 0.0
        # for each CDS in
        for pos in ind.blocks[:-1]:
            block_sequence = ind.sequence[pos[0]: pos[1] + self.junction_size]
            for motive in self.repetition_table.motives:
                cum_sum += block_sequence.count(motive)

        pos = ind.blocks[-1]
        block_sequence = ind.sequence[pos[0]: pos[1]]
        for motive in self.repetition_table.motives:
            cum_sum += block_sequence.count(motive)
        return cum_sum

    def __repr__(self):
        return 'Motifs & Repeats'
