import random
import string
import time

import cipher_common as cc


UPPERCASE_LETTERS = string.ascii_uppercase

class Genotype():
    def __init__(self, genes=None):
        # if genes is not specified, genes will be random
        if genes is None:
            possible_genes = UPPERCASE_LETTERS
            # shuffle possible genes
            self.genes = ''.join(random.sample(
                possible_genes, 
                len(possible_genes)))
        else:
            self.genes = genes

    def __crossover_gens(self, gen1, gen2, letters):
        if gen1 not in letters and \
                gen2 not in letters:
            return letters.pop()
        elif gen1 not in letters:
            return gen2
        elif gen2 not in letters:
            return gen1
        else:
            return random.choice((gen1,gen2))

    def crossover(self, indiv):
        letters = set(UPPERCASE_LETTERS)
        genes = list()
        for i in range(len(self.genes)):
            new_gen = self.__crossover_gens(
                indiv.genes[i], self.genes[i], letters.copy())
            letters.remove(new_gen)
            genes.append(new_gen)

        return Genotype(''.join(genes))

    def mutate(self):
        i, j = random.sample(range(len(self.genes)), 2)
        x, y = self.genes[i], self.genes[j]
        self.genes = self.genes[:j] + x + self.genes[j + 1:]
        self.genes = self.genes[:i] + y + self.genes[i + 1:]
        return self.genes

    def __repr__(self):
        return str(self.genes)


class NaturalSelection():
    def __init__(self, 
        ciphertext, popul_num, kill_percent=0.5, mutation_chance=0.1,
        exp_index=-4.202, enable_logging=False):
        # chance to mutate
        self.mutation_chance = mutation_chance
        # etalon quadgram index for english text
        self.exp_index = exp_index
        self.ciphertext = ciphertext
        self.popul_num = popul_num
        self.saved_lives_num = self.popul_num - int(self.popul_num * kill_percent)
        # init population with random genes
        self.popul = [Genotype() for _ in range(self.popul_num)]
        # create fitness list for population and init it
        self.fitness_list = list()

        # counter of generations
        self.generation = 0
        self.enable_logging = enable_logging


    def translate_cipher(self, key):
        # translate our cipher text by genes of individual
        key_map = str.maketrans(UPPERCASE_LETTERS, key)
        return self.ciphertext.translate(key_map)


    def fitness(self, indiv):
        decoded = self.translate_cipher(indiv.genes)
        index = cc.quadgram_index(decoded)
        fitness = abs(index - self.exp_index)
        return fitness


    def _update_fitness_list(self):
        # create a list of fitness values for population
        self.fitness_list = list(map(self.fitness, self.popul))


    def _choose_parents(self):
        # select two random and uniq parents from population
        parent1, parent2 = random.sample(self.popul, 2)
        return parent1, parent2


    def _make_selection(self):
        # sort fitness list and corresponding individuals
        self.fitness_list, self.popul = zip(
            *sorted(
                zip(self.fitness_list, self.popul), 
                key=lambda x: x[0]))
        
        # delete the worst part
        self.fitness_list = self.fitness_list[:self.saved_lives_num]
        self.popul = self.popul[:self.saved_lives_num]


    def _crossover_populaton(self):
        new_popul = list(self.popul)
        # crossover population until max number of individuals reached
        while len(new_popul) < self.popul_num:
            parent1, parent2 = self._choose_parents()
            new_popul.append(parent1.crossover(parent2))

        self.popul = new_popul


    def _mutate_popul(self):
        for indiv in self.popul:
            if random.choices(
                    [True, False], 
                    (self.mutation_chance, 1 - self.mutation_chance))[0]:
                indiv.mutate()


    def print_stats(self):
        print('Generation #', self.generation, 
                'Best fitness:', self.fitness_list[0])


    def solve(self, eps=0.05):
        start = time.time()

        self._update_fitness_list()
        while self.fitness_list[0] >= eps: 
            self._update_fitness_list()
            self._make_selection()
            self._crossover_populaton()
            self._mutate_popul()
            self.generation += 1
            if self.generation % 10 == 0:
                self.print_stats()

        key = self.popul[0].genes
        text = self.translate_cipher(key)
        print('Runtime:', time.time() - start)
        return text, key
