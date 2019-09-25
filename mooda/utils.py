from Bio.SeqFeature import SeqFeature
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import copy
import random


class GenbankFile():
    def __init__(self, genbank_record):
        self.genbank_record = genbank_record
        self.overlap_check = None
        self.genetic_code_check = None
        self.cds_list = []
        self.dict_features_index = {}

    """
    Initialising a dictionary feature:index and a list of all CDS index
    """

    def __get_cds_features_index_list(self):
        feature_index = 0
        for feature in self.genbank_record.features:
            self.dict_features_index[feature] = feature_index
            feature_index += 1
            if feature.type == "CDS":
                self.cds_list.append(feature_index - 1)

    """
    Check if source and destination genetic code are the same, return True if they are, return False if
    they are not
    """

    def __genetic_code_check(self, options):
        source_genetic_code = None
        for cds_index in range(0, len(self.cds_list[:-1])):
            cds = self.genbank_record.features[self.cds_list[cds_index]]

            if 'transl_table' in cds.qualifiers:
                source_genetic_code = int(cds.qualifiers['transl_table'][0])
            else:
                source_genetic_code = 1
        if source_genetic_code != options.target_genetic_code:
            self.genetic_code_check = False
        else:
            self.genetic_code_check = True

    """
    Given a genbank cheks if there are overlapping CDS
    """

    def __overlap_check(self):
        self.overlap_check = True
        for cds_index in range(0, len(self.cds_list[:-1])):
            cds = self.genbank_record.features[self.cds_list[cds_index]]
            next_cds = self.genbank_record.features[self.cds_list[cds_index+1]]
            if next_cds.location.start < cds.location.end:
                self.overlap_check = False
                break

    """
    Generates a random DNA sequence without
    starting codons to insert between overlapping CDS
    """
    def random_sequence_generator(self):
        nucleotides = ["A", "T", "C", "G"]
        patch_dna_string = []
        for base in range(100):
            random_nucleotide = random.choice(nucleotides)
            patch_dna_string.append(random_nucleotide)
        patch_dna_string = "".join(patch_dna_string)
        stop_codons = ['TAA', 'TGA', 'TAG']
        if "ATG" in patch_dna_string:
            patch_dna_string = patch_dna_string.replace('ATG', random.choice(stop_codons))
        patch_dna_string = Seq(patch_dna_string, IUPAC.unambiguous_dna)
        return patch_dna_string

    """
    If source and target genetic codes are different and there are overlapping CDS 
    Remove overlaps within two consecutives CDS by adding the overlap at the beginning
    of the next CDS
    """
    def __remove_overlaps(self):
        for cds_index in range(0, len(self.cds_list[:-1])):
            cds = self.genbank_record.features[self.cds_list[cds_index]]
            next_cds = self.genbank_record.features[self.cds_list[cds_index + 1]]
            cds_end = cds.location.end
            next_cds_start = next_cds.location.start
            # if there is an overlap add between two CDS A and B, add patch DNA with length equal to the overlap between A and B
            # and copy the overlap from A to B
            if next_cds_start < cds_end:
                overlap = cds_end - next_cds_start
                patch_dna = self.random_sequence_generator()
                self.genbank_record.seq = (
                    self.genbank_record.seq[:cds.location.end]
                    + patch_dna
                    + self.genbank_record.seq[cds.location.end - overlap:]
                )
                # Uptating CDS locations
                self.__update_locations(cds_index, overlap, patch_dna)

    """
    Method to update genbank locations
    """
    def __update_locations(self,cds_index, overlap, patch_dna):
        for forward_cds_index in self.cds_list[cds_index + 1:]:
            forward_cds = (self.genbank_record.features[forward_cds_index])
            feature_to_change = copy.deepcopy(forward_cds)
            new_feature_location = SeqFeature.FeatureLocation(start=forward_cds.location.start + overlap+ len(patch_dna),
                                                              end=forward_cds.location.end + overlap+len(patch_dna),
                                                              strand=forward_cds.location.strand)
            feature_to_change.location = new_feature_location  # change old feature location
            del self.genbank_record.features[forward_cds_index]  # remove changed feature
            self.genbank_record.features.insert(forward_cds_index,
                                       feature_to_change)  # add identical feature with new location in place of previous feature

    """
    Putting all togheter checks if there are overlaps and different genetic codes
    
    """
    def genbank_verification(self, options):
        self.__get_cds_features_index_list()
        self.__genetic_code_check(options)
        self.__overlap_check()
        self.__remove_overlaps()
