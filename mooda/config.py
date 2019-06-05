"""
    Class CodonUsageTable

    desc
"""
import yaml
from Bio.SeqUtils import GC


class YamlConfig:
    def __init__(self):
        pass

    def load_config(self, filename):
        with open(filename, "r") as config_file:
            config = yaml.load(config_file, Loader=yaml.FullLoader)
            return config


class CodonTable(YamlConfig):
    def __init__(self):
        super().__init__()
        self.codons = dict()

    def get_codons(self, aa, order_asc_gc=True):
        if order_asc_gc:
            return sorted(self.codons[aa], key=lambda k: GC(k))
        else:
            return sorted(self.codons[aa], key=lambda k: -GC(k))

    def get_codon_probabilities(self, aa):
        return self.codons[aa]


class RepetitionTable(YamlConfig):
    def __init__(self):
        super().__init__()
        self.motives_dict = dict()
        self.motives = []

    def get_repetition_len(self, rep):
        return self.motives_dict[rep]

    def get_motives_list(self):
        for unit in self.motives_dict:
            whole_motive = unit * self.motives_dict[unit]
            self.motives.append(whole_motive)

