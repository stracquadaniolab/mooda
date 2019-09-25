import numpy as np
import random
import math
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from mooda.config import *

"""
    Class Operator

    desc
"""


class Operator:
    def __init__(self, yaml):
        self.yaml = yaml


    def apply(self, ind):
        pass


"""
Class Block initialisator
"""


class BlockInitialiserOperator:
    """
    This method is used during the initialisation to create block : min_size< size<max_size
    """

    def __init__(self):
        super().__init__()
        self.max_block_size = None
        self.min_block_size = None
        self.step_size = None
        self.junction_size = None

    def _get_block_maxsize(self, ind):
        max_length_block = 0
        for block in ind.blocks:
            blocksize = block[1] - block[0]
            if blocksize > max_length_block:
                max_length_block = blocksize
        return max_length_block

    def apply(self, ind):
        # getting the block with the highest length
        max_length_block = self._get_block_maxsize(ind)
        # while at least one block is >than the preset max size
        while max_length_block > self.max_block_size:
            # choose a random block index
            bkpt = random.randint(0, len(ind.blocks) - 1)
            # evalaute his blocksize
            blocksize = ind.blocks[bkpt][1] - ind.blocks[bkpt][0]
            # if that block > than the preset max size
            if blocksize > self.max_block_size:
                # split the block according old_block = [min- max]
                # newblock = [min- random(min-max,min_block_size)][random(min-max,min_block_size)-max]
                second_end = ind.blocks[bkpt][1]
                first_end = random.randrange(
                    ind.blocks[bkpt][0] + self.min_block_size,
                    second_end - self.min_block_size, self.step_size)
                second_start = first_end
                ind.blocks[bkpt][1] = first_end
                newblock = [second_start, second_end]
                ind.blocks.insert((bkpt + 1), newblock)
            max_length_block = self._get_block_maxsize(ind)


"""
    Class SplitBlock

    desc
"""


class SplitBlockOperator(Operator):

    def __set_max_block_size(self):
        self.max_block_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.SplitBlockOperator"
        ]["max_block_size"]

    def __set_min_block_size(self):
        self.min_block_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.SplitBlockOperator"
        ]["min_block_size"]

    def __set_step_size(self):
        self.step_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.SplitBlockOperator"
        ]["step_size"]

    def initialise(self):
        self.__set_max_block_size()
        self.__set_min_block_size()
        self.__set_step_size()

    # split a block in two non overlapping
    def apply(self, ind):

        # making sure that blocks are sorted
        bkpts = list(range(0, len(ind.blocks)))
        bkpt = random.choice(bkpts)
        blocksize = ind.blocks[bkpt][1] - ind.blocks[bkpt][0]

        if blocksize >= 2 * self.min_block_size:
            second_end = ind.blocks[bkpt][1]

            try:
                first_end = random.randrange(ind.blocks[bkpt][0] + self.min_block_size,
                                             second_end - self.min_block_size, self.step_size)
            except ValueError:
                first_end = (ind.blocks[bkpt][0] + self.min_block_size)

            second_start = first_end
            second_end = ind.blocks[bkpt][1]
            ind.blocks[bkpt][1] = first_end
            newblock = [second_start, second_end]
            ind.blocks.insert((bkpt + 1), newblock)




"""
    Class JoinBlockOperator

    desc
"""


class JoinBlockOperator(Operator):

    def __set_max_block_size(self):
        self.max_block_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.JoinBlockOperator"
        ]["max_block_size"] - self.junction_size

    def __set_min_block_size(self):
        self.min_block_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.JoinBlockOperator"
        ]["min_block_size"]

    def __set_step_size(self):
        self.step_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.JoinBlockOperator"
        ]["step_size"]

    def __set_junction_size(self):
        self.junction_size = self.yaml["Algorithm"]["operators"][
            "mooda.operator.JoinBlockOperator"
        ]["junction_size"]


    def join(self, ind, block_1_pt, direction):
        # If forward join a block with next block
        if direction == "forward":
            block_2_pt = block_1_pt + 1
            ind.blocks[block_1_pt][1] = ind.blocks[block_2_pt][1]
            ind.blocks.remove(ind.blocks[block_2_pt])
            first_start = ind.blocks[block_1_pt][0]
            first_end = ind.blocks[block_1_pt][1]
            new_block_pt = block_1_pt


        # If backward join a block with previous
        elif direction == "backward":
            block_2_pt = block_1_pt - 1
            ind.blocks[block_2_pt][1] = ind.blocks[block_1_pt][1]
            ind.blocks.remove(ind.blocks[block_1_pt])
            first_start = ind.blocks[block_2_pt][0]
            first_end = ind.blocks[block_2_pt][1]
            new_block_pt = block_2_pt

        # if blocksize > max size split until  min<size<max
        block_size = first_end - first_start
        if first_end - first_start > self.max_block_size:
            break_range = math.floor(self.max_block_size - (block_size/2))
            median_point = math.floor((first_start + first_end)/2)
            if self.max_block_size - block_size > self.min_block_size:
                first_end = random.randrange((median_point - break_range),(median_point + break_range), self.step_size)
            else :
                first_end = median_point
            second_start = first_end
            second_end = ind.blocks[new_block_pt][1]
            ind.blocks[new_block_pt][1] = first_end
            newblock = [second_start, second_end]
            ind.blocks.insert((new_block_pt + 1), newblock)


    def initialise(self):
        self.__set_junction_size()
        self.__set_step_size()
        self.__set_max_block_size()
        self.__set_min_block_size()


    def apply(self, ind):
        block_1_pt = random.randrange(0, len(ind.blocks) - 1)
        if block_1_pt == 0:
            direction = "forward"
            self.join(ind, block_1_pt, direction)
        elif 0 < block_1_pt < (len(ind.blocks) - 1):
            direction = random.choice(["forward", "backward"])
            self.join(ind, block_1_pt, direction)
        elif block_1_pt == (len(ind.blocks) - 1):
            direction = "forward"
            self.join(ind, block_1_pt, direction)


"""
    Class GCOptimizationOperator(Operator)
"""


class GCOptimizationOperator(Operator):
    def __set_codon_GC_table(self):
        self.codon_GC_table = CodonTable()
        self.codon_GC_table.codons = self.codon_GC_table.load_config(
            self.yaml["Algorithm"]["operators"]["mooda.operator.GCOptimizationOperator"][
                "codon_GC_table"
            ]
        )

    def initialise(self):
        self.__set_codon_GC_table()
        self.target_gc = self.yaml["Algorithm"]["operators"][
            "mooda.operator.GCOptimizationOperator"
        ]["target_gc"]
        self.step_size_gc = self.yaml["Algorithm"]["operators"][
            "mooda.operator.GCOptimizationOperator"
        ]["step_size"]

    def apply(self, ind):
        # reference sequence

        # picking a CDS at random
        cds_pt = random.randint(0, (len(ind.cds_list) - 1))
        # compute the GC of the CDS
        cds_gc = GC(ind.cds_list[cds_pt].seq)
        # get codons for the CDS
        # create and shuffle CDS indexes for hill climbing
        cds_codon_index = list(range(0, len(ind.cds_list[cds_pt].codons)))
        np.random.shuffle(cds_codon_index)

        # hill climbing GC mutator
        codon_it = 0
        accepted_step = False
        while not accepted_step and codon_it < len(ind.cds_list[cds_pt].codons):
            # picking a codon at random
            curr_codon = ind.cds_list[cds_pt].codons[cds_codon_index[codon_it]]
            #checking that the codon is a triplet
            if len(curr_codon) != 3 or codon_it in ind.cds_list[cds_pt].overlapping_codons_indexes:
            # translate codon to aminoacid
                break
            curr_aa = curr_codon.translate(table=ind.cds_list[cds_pt].translation_table_target)

        # get current codons and directionality
        # in a hill climbing fashion 333, 5, 4, 4,
            if cds_gc < self.target_gc:
                curr_aa_codons = self.codon_GC_table.get_codons(
                    str(curr_aa), order_asc_gc=False
                )
            elif cds_gc > self.target_gc:
                curr_aa_codons = self.codon_GC_table.get_codons(
                    str(curr_aa), order_asc_gc=False
                )
            else:
                curr_aa_codons = None

        # computing stop criterion
            if curr_aa_codons and (curr_aa_codons[0] != curr_codon):
                ind.cds_list[cds_pt].codons[cds_codon_index[codon_it]] = Seq(
                    curr_aa_codons[0], IUPAC.unambiguous_dna
                )
                # calculating GC just for the CDS.
                ind.cds_list[cds_pt].build_sequence_from_codons()

                accepted_step = (
                    np.abs(GC(ind.cds_list[cds_pt].seq) - cds_gc) >= self.step_size_gc
                )
            else:
                codon_it = codon_it + 1

            # updating individual sequence
        ind.cds_list[cds_pt].build_sequence_from_codons()
        edited_cds = ind.cds_list[cds_pt].seq
        if ind.cds_list[cds_pt].strand != 1:
            edited_cds = edited_cds.complement()
            edited_cds = edited_cds[::-1]

        ind.sequence = (
            ind.sequence[: ind.cds_list[cds_pt].pt.location.start]
            + edited_cds
            + ind.sequence[ind.cds_list[cds_pt].pt.location.end:]
        )


"""
    Class CodonUsageOperator
"""


class CodonUsageOperator(Operator):
    def __set_codon_usage_table(self):
        self.codon_usage_table = CodonTable()
        self.codon_usage_table.codons = self.codon_usage_table.load_config(
            self.yaml["Algorithm"]["operators"]["mooda.operator.CodonUsageOperator"][
                "codon_usage_table"
            ]
        )

    def __set_step_size_codon_usage(self):
        self.step_size_codon_usage = (self.yaml["Algorithm"]["operators"]["mooda.operator.CodonUsageOperator"][
            "step_size"]) * 100

    def __swap_number_(self, codon_list):
        codons_to_swap = int(((len(codon_list)) * self.step_size_codon_usage) / 100)
        return codons_to_swap

    def initialise(self):
        self.__set_codon_usage_table()
        self.__set_step_size_codon_usage()

    def apply(self, ind):
        # reference sequence
        # checking that a CDS exists
        # picking a CDS at random
        cds_pt = random.randint(0, (len(ind.cds_list) - 1))
        # picking the CDS sequence
        # get codons for the CDS
        # number of indexes
        # evaluate the number of codons to edit
        codons_to_swap = self.__swap_number_(ind.cds_list[cds_pt].codons)
        codons_indexes = random.sample(
            range(0, (len(ind.cds_list[cds_pt].codons) - 1)), codons_to_swap
        )
        for index in codons_indexes:
            chosen_codon = ind.cds_list[cds_pt].codons[index]
            if len(chosen_codon) == 3 and index not in ind.cds_list[cds_pt].overlapping_codons_indexes:
                aa = chosen_codon.translate(table=ind.cds_list[cds_pt].translation_table_target)
                new_codon = np.random.choice(
                    list(self.codon_usage_table.codons[aa].keys()),
                    p=list(self.codon_usage_table.codons[aa].values()),
                )
                new_codon = Seq(new_codon, IUPAC.unambiguous_dna)
                ind.cds_list[cds_pt].codons[index] = new_codon

        ind.cds_list[cds_pt].build_sequence_from_codons()
        edited_cds =ind.cds_list[cds_pt].seq

        if ind.cds_list[cds_pt].strand != 1:
            edited_cds = edited_cds.complement()
            edited_cds = edited_cds[::-1]

        ind.sequence = (
            ind.sequence[: ind.cds_list[cds_pt].pt.location.start]
            + edited_cds
            + ind.sequence[ind.cds_list[cds_pt].pt.location.end:]
        )

"""
    Class Repetition
"""


class MotifOperator(Operator):

    def __set_repetition_table(self):
        self.repetition_table = RepetitionTable()
        self.repetition_table.motives_dict = self.repetition_table.load_config(
            self.yaml["Algorithm"]["operators"]["mooda.operator.MotifOperator"][
                "motif_table"
            ]
        )
        self.repetition_table.get_motives_list()

    def __set_codon_usage_table(self):
        self.codon_usage_table = CodonTable()
        self.codon_usage_table.codons = self.codon_usage_table.load_config(
            self.yaml["Algorithm"]["operators"]["mooda.operator.MotifOperator"][
                "codon_usage_table"
            ]
        )

    def __number_of_removal(self, repetition_list, cds):
        repetitions = 0
        for motive in repetition_list:
            n_rep = cds.count(motive)
            repetitions += n_rep
        removal_number = repetitions * self.step_size_rep
        return int(removal_number)

    def __remove_cds_without_rep(self, ind, motives_list):
        pt_list = []
        for motive in motives_list:
            for pt in range(0, len(ind.cds_list)):
                if pt not in pt_list and ind.cds_list[pt].seq.count(motive) != 0:
                    pt_list.append(pt)
        return pt_list

    def __select_repetion(self, repeats_list, ind, pt):
        motives = [
            motive for motive in repeats_list if ind.cds_list[pt].seq.count(motive) != 0
        ]
        return motives

    def __remove_motive(self, motive, ind, cds_pt):
        if ind.cds_list[cds_pt].seq.count(motive) > 0:
            cds_string = str(ind.cds_list[cds_pt].seq)
            start_rep = cds_string.index(motive)
            # getting the end index
            end_rep = start_rep + len(motive)
            # mapping this indexis on the list of codons
            start_codon = math.floor(start_rep / 3)
            end_codon = math.ceil(end_rep / 3) - 1
            # if the repetition is longer then 2 codon avoid the first and the last codon
            selected_codon_index = random.randint(start_codon, end_codon)
            selected_codon = ind.cds_list[cds_pt].codons[selected_codon_index]
            if len(selected_codon) == 3 and selected_codon_index not in ind.cds_list[cds_pt].overlapping_codons_indexes:
                aa = selected_codon.translate(table=ind.cds_list[cds_pt].translation_table_target)
                synonimous_codons = list(self.codon_usage_table.codons[aa].keys())
                if len(synonimous_codons) > 1:
                    synonimous_codons.remove(selected_codon)
                    new_codon = random.choice(synonimous_codons)
                    new_codon = Seq(new_codon, IUPAC.unambiguous_dna)
                    ind.cds_list[cds_pt].codons[selected_codon_index] = new_codon
                    ind.cds_list[cds_pt].build_sequence_from_codons()

    def initialise(self):
        self.__set_repetition_table()
        self.__set_codon_usage_table()
        self.step_size_rep = (self.yaml["Algorithm"]["operators"]["mooda.operator.MotifOperator"]
        ["step_size"]) * 100

    def apply(self, ind):
        pts_motives_cds = self.__remove_cds_without_rep(
            ind, self.repetition_table.motives
        )
        if len(pts_motives_cds) > 0:
            cds_pt = random.choice(pts_motives_cds)
            to_remove = self.__select_repetion(
                self.repetition_table.motives, ind, cds_pt
            )
            op_iter = self.__number_of_removal(to_remove, ind.cds_list[cds_pt].seq)
            for ip_iter in range(op_iter):
                motive = random.choice(to_remove)
                self.__remove_motive(motive, ind, cds_pt)

            if ind.cds_list[cds_pt].strand != 1:
                ind.cds_list[cds_pt].seq = ind.cds_list[cds_pt].seq.complement()
                ind.cds_list[cds_pt].seq = ind.cds_list[cds_pt].seq[::-1]
            ind.sequence = (
                ind.sequence[: ind.cds_list[cds_pt].pt.location.start]
                + ind.cds_list[cds_pt].seq
                + ind.sequence[ind.cds_list[cds_pt].pt.location.end:]
            )
