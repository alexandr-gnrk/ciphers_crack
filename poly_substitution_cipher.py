import random
import string
import time
import itertools
import threading
from pprint import pprint

import cipher_common as cc


UPPERCASE_LETTERS = string.ascii_uppercase

# >>> import cipher_common as cc
# >>> cc.monogram_index(text)
# - 1.2571025627576808
# >>> cc.bigram_index(text)
# -2.3542927234721676
# >>> cc.trigram_index(text)
# -3.341533485572007
# >>> cc.quadgram_index(text)
# -4.2221607812359405

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
        result_genes = ''
        letters = set(UPPERCASE_LETTERS)

        for i in range(len(self.genes)):
            new_gen = self.__crossover_gens(
                indiv.genes[i], self.genes[i], letters.copy())
            letters.remove(new_gen)
            result_genes += new_gen

        return Genotype(result_genes)

    def mutate(self):
        i, j = random.sample(range(len(self.genes)), 2)
        x, y = self.genes[i], self.genes[j]
        self.genes = self.genes[:j] + x + self.genes[j + 1:]
        self.genes = self.genes[:i] + y + self.genes[i + 1:]
        return self.genes

    def __repr__(self):
        return str(self.genes)


class ParallelGeneticAlgorithm():

    def  __init__(self,
            ciphertext, gen_pos, gens_num, 
            best_gens, best_gens_lock, comm_freq, 
            popul_num, kill_percent=0.5, mutation_chance=0.1):
        self.ciphertext = ciphertext
        # position of current alphabet in list of alphabets
        self.gen_pos = gen_pos
        # amount of alphabets
        self.gens_num = gens_num
        self.best_gens = best_gens
        self.best_gens_lock = best_gens_lock
        self.best_gens_copy = None
        # after how many generations GAs will communicate
        self.comm_freq = comm_freq
        self.popul_num = popul_num
        self.saved_lives_num = self.popul_num - int(self.popul_num * kill_percent)
        # list of two nei
        self.kill_percent = kill_percent
        self.mutation_chance = mutation_chance

        self.generation = 0
        self.is_warming_up = False
        # calc fitness function for monograms by default
        self.popul = [Genotype() for _ in range(self.popul_num)]
        self.fitness_list = list()

    def best_gens_push(self, indiv):
        self.best_gens_lock.acquire()
        self.best_gens[self.gen_pos] = indiv
        self.best_gens_lock.release()

    def best_gens_pull(self):
        self.best_gens_lock.acquire()
        cp = self.best_gens[:]
        self.best_gens_lock.release()
        return cp

    def warm_up(self):
        self.is_warming_up = True
        for i in range(self.comm_freq):
            self.__update_fitness_list()
            self.__make_selection()
            self.__crossover_populaton()
            self.__mutate_popul()
            self.generation += 1
            self.print_stats()
        self.is_warming_up = False

        self.best_gens_push(self.popul[0])

    def ngram_index(self, groups, group_size):
        return sum(cc.ngram_index(group, group_size) for group in groups) / len(groups)
        
    def translate_cipher(self, key):
        key_maps = list()
        for i in range(self.gens_num):
            if i == self.gen_pos:
                key_maps.append(str.maketrans(UPPERCASE_LETTERS, key))
            else:
                key_maps.append(
                    str.maketrans(
                        UPPERCASE_LETTERS,
                        self.best_gens_copy[i].genes))

        res_str = ''
        maps_iter = itertools.cycle(key_maps)
        
        key_map = next(maps_iter)
        for char in self.ciphertext:
            res_str += char.translate(key_map)
            key_map = next(maps_iter)

        return res_str

    def get_mono_groups(self):
        groups = list()
        for i in range(self.gen_pos, len(self.ciphertext), self.gens_num):
            groups.append(self.ciphertext[i])
        return groups

    def translate_mono_groups(self, groups, key):
        res = list()
        key_map = str.maketrans(UPPERCASE_LETTERS, key)
        for group in groups:
            res.append(group.translate(key_map))
        return res

    def fitness(self, indiv):
        if self.is_warming_up:
            groups = self.get_mono_groups()
            decoded_groups = self.translate_mono_groups(groups, indiv.genes)
            index = self.ngram_index(decoded_groups, 1)
            fitness = abs(index + 1.2571025627576808)
            return fitness

        decoded = self.translate_cipher(indiv.genes)
        index = cc.quadgram_index(decoded)
        fitness = abs(index + 4.2221607812359405)
        # 2221607812359405
        return fitness

    def __update_fitness_list(self):
        self.fitness_list = list()

        # create a list of fitness values for population
        for indiv in self.popul:
            fitness = self.fitness(indiv)
            self.fitness_list.append(fitness)

    def __choose_parents(self):
        # select two random and uniq parents from population
        parent1, parent2 = random.sample(self.popul, 2)
        return parent1, parent2

    def __make_selection(self):
        # sort fitness list and corresponding individuals
        self.fitness_list, self.popul = zip(
            *sorted(
                zip(self.fitness_list, self.popul), 
                key=lambda x: x[0]))
        
        # delete the worst part
        self.fitness_list = self.fitness_list[:self.saved_lives_num]
        self.popul = self.popul[:self.saved_lives_num]

    def __crossover_populaton(self):
        new_popul = list(self.popul)
        # crossover population until max number of individuals reached
        while len(new_popul) < self.popul_num:
            parent1, parent2 = self.__choose_parents()
            new_popul.append(parent1.crossover(parent2))

        self.popul = new_popul

    def __mutate_popul(self):
        # start = int(0.1 * len(self.popul))
        # if start == 0 if len(self.po)
        for indiv in self.popul:
            if random.choices(
                    [True, False], 
                    (self.mutation_chance, 1 - self.mutation_chance))[0]:
                indiv.mutate()

    def print_stats(self):
        print('Alphabet #{} Generation #{} Best fitness: {}'.format(
            self.gen_pos,
            self.generation, 
            self.fitness_list[0]))

    def solve(self, eps=0.17):
        self.best_gens_copy = self.best_gens_pull()
        self.__update_fitness_list()
        # input()
        while self.fitness_list[0] >= eps:
            self.__update_fitness_list()
            self.__make_selection()
            self.__crossover_populaton()
            self.__mutate_popul()

            if self.generation % self.comm_freq == 0:
                self.best_gens_push(self.popul[0])
                self.best_gens_copy = self.best_gens_pull()

            self.generation += 1
            self.print_stats()

        self.best_gens_lock.acquire()
        self.best_gens = self.best_gens_copy[:]
        self.best_gens[self.gen_pos] = self.popul[0]
        self.best_gens_lock.release()

        key = self.popul[0].genes
        text = self.translate_cipher(key)
        print("""========== GA#{} ended with ff={} ==========
            ------ TEXT ------
            {}
            ------ KEYS ------
            best_gens_copy
            {}
            self.key = {}""".format(
                self.gen_pos,
                self.fitness_list[0],
                text,
                self.best_gens_copy,
                self.popul[0]))
        return text, key


class NaturalSelection():
    def __init__(self, 
            ciphertext, popul_num, alphabets_num, comm_freq, kill_percent=0.5, mutation_chance=0.1):
        # etalon quadgram index for english text
        self.ciphertext = ciphertext
        
        self.best_gens = [None for _ in range(alphabets_num)]
        self.best_gens_lock = threading.Lock()

        self.genetic_algorithms = list()
        for i in range(alphabets_num):
            self.genetic_algorithms.append(
                ParallelGeneticAlgorithm(
                    self.ciphertext,
                    i,
                    alphabets_num,
                    self.best_gens,
                    self.best_gens_lock,
                    comm_freq,
                    popul_num,
                    kill_percent,
                    mutation_chance))

    def print_stats(self):
        print('Generation #', self.generation, 
                'Best fitness:', self.fitness_list[0])

    def warm_up(self):
        threads = list()
        for ga in self.genetic_algorithms:
            threads.append(
                threading.Thread(
                    target=ParallelGeneticAlgorithm.warm_up,
                    args=(ga,)))
            threads[-1].start()

        for thread in threads:
            thread.join()


    def solve(self):
        threads = list()
        for ga in self.genetic_algorithms:
            threads.append(
                threading.Thread(
                    target=ParallelGeneticAlgorithm.solve,
                    args=(ga,)))
            threads[-1].start()

        for thread in threads:
            thread.join()

        best_ga = self.genetic_algorithms[0]
        for ga in self.genetic_algorithms[1:]:
            if ga.fitness_list[0] < best_ga.fitness_list[0]:
                best_ga = ga

        text = best_ga.translate_cipher(best_ga.popul[0].genes)
        # print("========== KEY ==========")
        # pprint(best_ga.best_gens_copy)
        # print(best_ga.popul[0])
        # print("========== TEXT ==========")
        # print(best_ga.translate_cipher(best_ga.popul[0].genes))
        return text


if __name__ == '__main__':
    CT = 'KZBWPFHRAFHMFSNYSMNOZYBYLLLYJFBGZYYYZYEKCJVSACAEFLMAJZQAZYHIJFUNHLCGCINWFIHHHTLNVZLSHSVOZDPYSMNYJXHMNODNHPATXFWGHZPGHCVRWYSNFUSPPETRJSIIZSAAOYLNEENGHYAMAZBYSMNSJRNGZGSEZLNGHTSTJMNSJRESFRPGQPSYFGSWZMBGQFBCCEZTTPOYNIVUJRVSZSCYSEYJWYHUJRVSZSCRNECPFHHZJBUHDHSNNZQKADMGFBPGBZUNVFIGNWLGCWSATVSSWWPGZHNETEBEJFBCZDPYJWOSFDVWOTANCZIHCYIMJSIGFQLYNZZSETSYSEUMHRLAAGSEFUSKBZUEJQVTDZVCFHLAAJSFJSCNFSJKCFBCFSPITQHZJLBMHECNHFHGNZIEWBLGNFMHNMHMFSVPVHSGGMBGCWSEZSZGSEPFQEIMQEZZJIOGPIOMNSSOFWSKCRLAAGSKNEAHBBSKKEVTZSSOHEUTTQYMCPHZJFHGPZQOZHLCFSVYNFYYSEZGNTVRAJVTEMPADZDSVHVYJWHGQFWKTSNYHTSZFYHMAEJMNLNGFQNFZWSKCCJHPEHZZSZGDZDSVHVYJWHGQFWKTSNYHTSZFYHMAEDNJZQAZSCHPYSKXLHMQZNKOIOKHYMKKEIKCGSGYBPHPECKCJJKNISTJJZMHTVRHQSGQMBWHTSPTHSNFQZKPRLYSZDYPEMGZILSDIOGGMNYZVSNHTAYGFBZZYJKQELSJXHGCJLSDTLNEHLYZHVRCJHZTYWAFGSHBZDTNRSESZVNJIVWFIVYSEJHFSLSHTLNQEIKQEASQJVYSEVYSEUYSMBWNSVYXEIKWYSYSEYKPESKNCGRHGSEZLNGHTSIZHSZZHCUJWARNEHZZIWHZDZMADNGPNSYFZUWZSLXJFBCGEANWHSYSEGGNIVPFLUGCEUWTENKCJNVTDPNXEIKWYSYSFHESFPAJSWGTYVSJIOKHRSKPEZMADLSDIVKKWSFHZBGEEATJLBOTDPMCPHHVZNYVZBGZSCHCEZZTWOOJMBYJSCYFRLSZSCYSEVYSEUNHZVHRFBCCZZYSEUGZDCGZDGMHDYNAFNZHTUGJJOEZBLYZDHYSHSGJMWZHWAFTIAAY'
    
    ns = NaturalSelection(CT, 300, 4, 10)
    print('+++++++++ Warming up +++++++++')
    ns.warm_up()
    print('+++++++++ Solving +++++++++')
    text = ns.solve()
    print(text)