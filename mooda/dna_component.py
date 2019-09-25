from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import Bio.Data.CodonTable
from Bio.Alphabet import IUPAC
import random

class DNAComponent:
    def __init__(self):
        self.pt = None
        self.seq = None


class Coding(DNAComponent):
    def __init__(self,):
        super().__init__()
        self.codons = None
        self.strand = None
        self.translation_table_origin = None
        self.translation_table_target = None
        self.start_codon = None
        self.translation = None
        self.origin_table_dict = None
        self.target_table_dict = None
        # list of overlapping codon to not edit by sequence operators
        self.overlapping_codons_indexes = []

    def __get_CDS(self, ind):
        self.seq = ind.sequence[self.pt.location.start : self.pt.location.end]
        if 'transl_table' in self.pt.qualifiers:
             self.translation_table_origin = int(self.pt.qualifiers['transl_table'][0])
        else:
            self.translation_table_origin = 1

        if 'codon_start' in self.pt.qualifiers:
            self.codon_start = self.pt.qualifiers['codon_start'][0]
        else:
            self.codon_start = 1

        if 'translation' in self.pt.qualifiers:
            self.translation = self.pt.qualifiers['translation'][0]

        if self.pt.location.strand and self.pt.location.strand != 1:
            self.seq = self.seq.complement()
            self.seq = self.seq[::-1]




    def __setting_origin_and_target_genetic_code(self,genetic_code_table):
        # reading target genetic code
        self.translation_table_target = genetic_code_table
        # if target genetic code != origin genetic code
        if self.translation_table_origin != self.translation_table_target:
            origin_table = Bio.Data.CodonTable.generic_by_id[self.translation_table_origin]
            #Turn the the genetic code table into a dictionary
            self.origin_table_dict = origin_table.__dict__
            for stop_codon in self.origin_table_dict["stop_codons"]:
                self.origin_table_dict["forward_table"][stop_codon]= '*'
            target_table = Bio.Data.CodonTable.generic_by_id[self.translation_table_target]
            self.target_table_dict = target_table.__dict__
            for stop_codon in self.target_table_dict["stop_codons"]:
                self.target_table_dict["forward_table"][stop_codon]= '*'

    def __recoding_genetic_code(self):
        if self.translation_table_origin != self.translation_table_target:
            for codon_index in range(0,len(self.codons)):
                if len(self.codons[codon_index]) == 3:
                    codon =self.codons[codon_index]
                    amino_acid_origin =(self.origin_table_dict["forward_table"][codon])
                    amino_acid_target =(self.target_table_dict["forward_table"][codon])
                    if amino_acid_origin != amino_acid_target:
                        synonyms_codons =[]
                        for key, value in self.target_table_dict["forward_table"].items():
                            if amino_acid_origin == value:
                                if "U" not in key:
                                    synonyms_codons.append(key)
                        self.codons[codon_index] = Seq(random.choice(synonyms_codons),IUPAC.unambiguous_dna)
        self.build_sequence_from_codons()





    def __get_codon_list(self):
        self.codons = [self.seq[i : i + 3] for i in range(0, len(self.seq), 3)]


    def build_sequence_from_codons(self):
        self.seq = sum(self.codons, Seq("", DNAAlphabet()))



    def intialise_cds(self, ind, pt,genetic_code_table):
        self.pt = pt
        self.strand = self.pt.location.strand
        self.__get_CDS(ind)
        self.__get_codon_list()
        self.__setting_origin_and_target_genetic_code(genetic_code_table)
        self.__recoding_genetic_code()





class NonCoding(DNAComponent):
    def __get_UTR(self, ind):
        self.seq = ind.sequence[self.pt[0] : self.pt[1]]

    def intialise_utr(self, ind, pt):
        self.pt = pt
        self.__get_UTR(ind)

