class Assembly:
    def __init__(self, yaml):
        self.yaml = yaml

    def apply(self, ind):
        pass


class Gibson(Assembly):
    def __set_junction_size(self):
        self.junction_size = self.yaml["Algorithm"]["assemblies"][
            "mooda.assembly.Gibson"
        ]["junction_size"]

    def initialise(self):
        self.__set_junction_size()

    """
    Takes as input all blocks if there is no overlap it
    makes it
    """

    def apply(self, population):
        # for each block
        for ind in population.individuals:
            for block_index in range(len(ind.blocks[:-1])):
                next_block = block_index + 1
                # if there is no overlap make it
                if ind.blocks[block_index][1] == ind.blocks[next_block][0]:
                    ind.blocks[block_index][1] = ind.blocks[block_index][1] + self.junction_size
