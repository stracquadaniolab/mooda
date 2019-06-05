'''
    Class Archive

    Represents a list of all non-dominated solutions for a multi-objective DNA assembly problem

'''


class Archive:
    def __init__(self):
        # instance variabe with list of solutions
        self.individuals = []

        # list of non-dominated fronts
        self.fronts = [[]]
        # Number of individuals allowed in the archive
        self.length = None

    def __len__(self):
        return len(self.individuals)

    def add_individual(self, ind):
        self.individuals.append(ind)

    def add_to_archive(self, population):
        # add to archive
        for individual in population.fronts[0]:
            if individual not in self.individuals:
                self.individuals.append(individual)

    def compute_ranking(self):
        self.fronts = [[]]
        for p in self.individuals:

            p.dominated_set = []
            p.dominating_count = 0
            for l in self.individuals:
                p_counter = 0
                for objective in range(0, (len(p.objectives))):
                    if p.objectives[objective] < l.objectives[objective]:
                        p_counter += 1
                    elif p.objectives[objective] > l.objectives[objective]:
                        p_counter = 0
                        break
                if p_counter >= 1:
                    p.dominated_set.append(l)
                l_counter = 0
                for objective in range(0, len(p.objectives)):
                    if l.objectives[objective] < p.objectives[objective]:
                        l_counter += 1
                    elif l.objectives[objective] > p.objectives[objective]:
                        l_counter = 0
                        break
                if l_counter >= 1:
                    p.dominating_count += 1
            if p.dominating_count == 0:
                p.rank = 1
                self.fronts[0].append(p)

    def initialise_fronts(self):
        counter = 1
        while len(self.fronts[counter - 1]) > 0:
            next_front = []
            for p in self.fronts[counter - 1]:
                for l in p.dominated_set:
                    l.dominating_count -= 1
                    if l.dominating_count == 0:
                        l.rank = counter + 1
                        next_front.append(l)
            counter += 1
            self.fronts.append(next_front)

    def compute_crowding_distance(self, objective_functions):
        if len(self.individuals) > 2:
            for individual in self.individuals:
                individual.crowding_distance = 0

            for obj in range(0, len(objective_functions)):
                self.individuals.sort(key=lambda individual: individual.objectives[obj])
                min = self.individuals[0].objectives[obj]
                max = self.individuals[-1].objectives[obj]
                for front in self.fronts:
                    if len(front) > 2:
                        front.sort(key=lambda individual: individual.objectives[obj])
                        front[0].crowding_distance = front[-1].crowding_distance = 1
                        for individual_index in range(1, len(front) - 1):
                            try:
                                front[individual_index].crowding_distance += abs(
                                    front[individual_index - 1].objectives[obj]
                                    - front[individual_index + 1].objectives[obj]
                                ) / (abs(min - max))
                            except ZeroDivisionError:
                                front[individual_index].crowding_distance += 0

    def remove_duplicates(self):
        new_individuals = []
        for ind in self.individuals:
            if ind not in new_individuals:
                new_individuals.append(ind)
        self.individuals = new_individuals

    def select_only_rank1(self):
        self.individuals = self.fronts[0]

    def random_deletion(self):
        if len(self.individuals) >= self.length:
            self.individuals.sort(
                key=lambda individual: (individual.rank, -individual.crowding_distance))
        while len(self.individuals) > self.length:
            del self.individuals[-1]
