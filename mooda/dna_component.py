from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet


class DNAComponent:
    def __init__(self):
        self.pt = None
        self.seq = None


class Coding(DNAComponent):
    def __init__(self):
        super().__init__()
        self.codons = None
        self.strand = None

    def __get_CDS(self, ind):
        self.seq = ind.sequence[self.pt.location.start : self.pt.location.end]
        if self.pt.location.strand != 1:
            self.seq = self.seq.complement()

    def __get_codon_list(self):
        self.codons = [self.seq[i : i + 3] for i in range(0, len(self.seq), 3)]

    def intialise_cds(self, ind, pt):
        self.pt = pt
        self.strand = self.pt.location.strand
        self.__get_CDS(ind)
        self.__get_codon_list()

    def build_sequence_from_codons(self):
        self.seq = sum(self.codons, Seq("", DNAAlphabet()))
        if self.strand != 1:
            self.seq = self.seq.complement()


class NonCoding(DNAComponent):
    def __get_UTR(self, ind):
        self.seq = ind.sequence[self.pt[0] : self.pt[1]]

    def intialise_utr(self, ind, pt):
        self.pt = pt
        self.__get_UTR(ind)

