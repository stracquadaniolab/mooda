import yaml
import pandas as pd
import os
from Bio.SeqRecord import SeqRecord
from mooda.evolutionary_algorithm import *
from pathlib import Path


class MoodaFile:


    def __init__(self, yaml_file):
        self.yaml_file = yaml_file


    def set_yaml_file(self, config):
        self.yaml_file = config

    def __get_path(self, output_dir):
        self.path = Path(output_dir)

    def __check_dir(self, output_dir):
        test_dir = os.path.isdir(output_dir)
        if not test_dir:
            os.makedirs(output_dir)




    def __file_check(self, output_dir, suffix):
        n = 1

        wildcard = self.path.parts[-1]
        filename = wildcard + suffix
        self.__check_dir(output_dir)
        while filename in os.listdir(output_dir):
            n += 1
            filename = (
                str(wildcard)
                + "_"
                + str(n)
                + suffix
            )
        filename = os.path.join(output_dir, filename)
        return filename

    def save_log(self, options, ag):
        suffix = "_logfile.yaml"
        self.__check_dir(options.output_dir)
        filename = self.path.parts[-1] + suffix
        output = os.path.join(options.output_dir,filename)
        parameters = dict()
        parameters["max_iterations"] = ag.max_iterations

        log = dict(
            {
                "Input file": options.input_file,
                "Outputs": options.output_dir,
                "Objective functions": ag.yaml["Algorithm"]["objective_functions"],
                "Operators": ag.yaml["Algorithm"]["operators"],
                "Max iterations": ag.max_iterations,
                "Output genbank and fasta": options.sequence_output,
                "GenBank data": str(ag.gebank_record),
                "Sequence": str(
                    ag.gebank_record.seq[0:100] + "...." + ag.gebank_record.seq[-100:]
                ),
                "start": str(ag.time_start),
                "time process": str(ag.time_end - ag.time_start),
                "end": str(ag .time_end),
            }
        )
        if options.algorithm == "mc":
            log["Algorithm"] = "Montecarlo"
        elif options.algorithm == "mo":
            log["Algorithm"] = "Multi-Objective"
            log["Population"] = ag.population.size
            log["Archive size"] = ag.archive.length
        yaml.dump(
            log, open(output, "w"), default_flow_style=False, allow_unicode=True
        )

    def save_mooda_result(self, options, ag):
        suffix = "_mooda_result.csv"
        self.__check_dir(options.output_dir)
        filename = self.path.parts[-1] + suffix
        output = os.path.join(options.output_dir, filename)
        objective_list = []
        for objective_function in ag.objective_functions:
            label = objective_function.__repr__()
            objective_list.append(label)
        fieldnames = ["Crowding distance"] + objective_list + ["Fragments", "Record"]

        individual_list = None
        if options.algorithm == "mo":
            individual_list = ag.population.individuals
        elif options.algorithm == "mc":
            individual_list = ag.archive.individuals

        dict_list = []
        for individual in individual_list:
            individual.objectives = [round(objective, 2) for objective in individual.objectives]
            objective_dictionary = dict(zip(objective_list, individual.objectives))
            objective_dictionary["Crowding distance"] = round(individual.crowding_distance, 2)
            objective_dictionary["Fragments"] = individual.fasta
            objective_dictionary["Record"] = individual.genbank
            dict_list.append(objective_dictionary)
        table = pd.DataFrame(dict_list, columns=fieldnames)
        table.to_csv(output)

    def save_fasta(self, options,ag):
        suffix = ".fasta"
        path_dir = Path(options.output_dir)
        wildcard = path_dir.parts[-1]
        if ag.outputfasta_gebank:
            individual_list =[]
            if options.algorithm == "mo":
                individual_list = ag.population.individuals
            elif options.algorithm == "mc":
                individual_list = ag.archive.individuals
            for ind in individual_list:
                ind.fasta = self.__file_check(options.output_dir, suffix)
                with open(ind.fasta, "w") as fastaseq:
                    for block in ind.blocks:
                        seq_record = SeqRecord(
                            ind.sequence[block[0] : block[1]],
                            id=ind.gebank_record.id,
                            name=ind.gebank_record.name,
                            description=">" + str(wildcard) + " " + str(block) + "\n",
                        )
                        SeqIO.write(seq_record, fastaseq, "fasta")

    def save_genbank(self, options, ag):
        suffix = ".genbank"
        if ag.outputfasta_gebank:
            individual_list = []
            if options.algorithm == "mo":
                individual_list = ag.population.individuals
            elif options.algorithm == "mc":
                individual_list = ag.archive.individuals
            for ind in individual_list:
                ind.genbank = self.__file_check(options.output_dir, suffix)
                record = SeqRecord(
                    ind.sequence,
                    id=ind.gebank_record.id,
                    name=ind.gebank_record.name,
                    description=ind.gebank_record.id,
                )
                for feature in ind.features:
                    record.features.append(feature)
                output_file = open(ind.genbank, "w")
                SeqIO.write(record, output_file, "genbank")

    def save(self,options, ag):
        self.__get_path(options.output_dir)
        self.save_fasta(options, ag)
        self.save_genbank(options, ag)
        self.save_log(options, ag )
        self.save_mooda_result(options, ag)

