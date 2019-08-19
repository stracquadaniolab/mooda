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
        self.translation_table = None
        self.start_codon = None
        self.translation = None

    def __get_CDS(self, ind):
        self.seq = ind.sequence[self.pt.location.start : self.pt.location.end]
        if 'transl_table' in self.pt.qualifiers:
             self.translation_table = self.pt.qualifiers['transl_table'][0]
        else:
            self.translation_table = 1
        if 'codon_start' in self.pt.qualifiers:
            self.codon_start = self.pt.qualifiers['codon_start'][0]
        else:
            self.codon_start = 1
        if 'translation' in self.pt.qualifiers:
            self.translation = self.pt.qualifiers['translation'][0]

        if self.pt.location.strand != 1:
            self.seq = self.seq.complement()
            self.seq = self.seq[::-1]
        # print(self.seq)
        # print(self.pt)
        # print(self.pt.location)
        # print('table',self.translation_table)
        # a = self.seq.translate(table=self.translation_table, to_stop=True)
        # b = self.translation
        # print(a)
        # print(len(a), len(b))
        # if a==b:
        #     print("everythin ok")
        # else :
        #     print('trnaslation error')




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
            self.seq = self.seq[::-1]



class NonCoding(DNAComponent):
    def __get_UTR(self, ind):
        self.seq = ind.sequence[self.pt[0] : self.pt[1]]

    def intialise_utr(self, ind, pt):
        self.pt = pt
        self.__get_UTR(ind)

