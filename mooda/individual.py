'''
    Class Individual

    Represents solution to a multi-objective DNA assembly problem

'''
import copy
from mooda.dna_component import Coding
from mooda.dna_component import NonCoding


class Individual:
    def __init__(self, gebank_record):
        # sequence information
        self.gebank_record = gebank_record
        self.sequence = self.gebank_record.seq
        self.features = self.gebank_record.features

        # list of blocks for assembly
        self.blocks = [[0, len(gebank_record.seq)]]
        # variables for the optimisation
        # objetctives
        self.objectives = []
        # rank of the solution
        self.rank = float("inf")
        # list of dominated solutions
        self.dominated_set = []
        # number of solutions that are dominated
        self.domination_count = 0
        # number of solutions that are better
        self.dominating_count = 0
        # crowding distance as defined by NSGAII
        self.crowding_distance = 0
        # list of CDS indexes included in the sequence
        self.cds_indexes_list = []
        # list of UTR indeces icluded in the sequence
        self.utr_indexes_list = []
        # list of codons
        self.cds_list = []
        # list of utr
        self.utr_list = []
        # fasta file
        self.fasta = None
        # genbank file
        self.genbank = None

    # representation of the solution
    def __str__(self):
        return str(self.objectives)

    def __eq__(self, other):
        if self.sequence == other.sequence \
            and self.blocks == other.blocks:
            return True
        else:
            return False

    # clone the instance with a deepcopy
    # useful when creating offsprings
    def clone(self):
        offspring = Individual(self.gebank_record)
        offspring.sequence = copy.deepcopy(self.sequence)
        offspring.features = copy.deepcopy(self.features)
        offspring.blocks = copy.deepcopy(self.blocks)
        offspring.cds_indexes_list = copy.deepcopy(self.cds_indexes_list)
        offspring.utr_indexes_list = copy.deepcopy(self.utr_indexes_list)
        offspring.cds_list = copy.deepcopy(self.cds_list)
        offspring.utr_list = copy.deepcopy(self.utr_list)
        return offspring

    def __get_feature_CDS_list(self):
        cds_index = []
        for feat in self.features:
            if feat.type == "CDS":
                cds_index.append(feat)
        self.cds_indexes_list = cds_index

    def __get_first_feature_UTR(self):
        if self.cds_indexes_list[0].location.start != 0:
            utr_pt = [0]
            utr_pt.append(self.cds_indexes_list[0].location.start - 1)
            self.utr_indexes_list.append(utr_pt)

    def __get_last_feature_UTR(self, pt):
        utr_pt = [None, None]
        utr_pt[0] = self.cds_indexes_list[pt].location.end + 1
        utr_pt[1] = len(self.sequence)
        if utr_pt[0] != utr_pt[1]:
            self.utr_indexes_list.append(utr_pt)

    def __get_feature_UTR_list(self):
        self.__get_first_feature_UTR()
        for pt in range(len(self.cds_indexes_list)):
            utr_pt = [None, None]
            if 0 <= pt < len(self.cds_indexes_list) - 1:
                utr_pt[0] = self.cds_indexes_list[pt].location.end + 1
                utr_pt[1] = self.cds_indexes_list[pt + 1].location.start - 1
                self.utr_indexes_list.append(utr_pt)
            elif pt == len(self.cds_indexes_list) - 1:
                self.__get_last_feature_UTR(pt)

    def __get_assembly(self, junction_size ):
        for block_index in range(len(self.blocks[:-1])):
            next_block = block_index + 1
            # if there is no overlap make it
            if self.blocks[block_index][1] == self.blocks[next_block][0]:
                self.blocks[block_index][1] = self.blocks[block_index][1] + junction_size

    def initialise(self):
        self.__get_feature_CDS_list()
        self.__get_feature_UTR_list()
        for pt in self.cds_indexes_list:
            cds = Coding()
            cds.intialise_cds(self, pt)
            self.cds_list.append(cds)
        for pt in self.utr_indexes_list:
            utr = NonCoding()
            utr.intialise_utr(self, pt)
            self.utr_list.append(utr)
