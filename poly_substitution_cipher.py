import random
import string
import time
import itertools
import threading
from pprint import pprint

import cipher_common as cc
from substitution_cipher import Genotype, NaturalSelection


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

class ParallelGeneticAlgorithm(NaturalSelection):

    def  __init__(self,
            ciphertext, gen_pos, gens_num, 
            best_gens, best_gens_lock, solved_lock,
            comm_freq, popul_num, 
            kill_percent=0.5, mutation_chance=0.1):
        super().__init__(
            ciphertext, popul_num, kill_percent, 
            mutation_chance, exp_index=-4.22216)
       
        self.exp_mono_index = -1.25710
        # position of current alphabet in list of alphabets
        self.gen_pos = gen_pos
        # amount of alphabets
        self.gens_num = gens_num
        # list of best genes and lock from natural selection class
        self.best_gens = best_gens
        self.best_gens_lock = best_gens_lock
        self.solved_lock = solved_lock
        # copy of best gens, which updates according to comm_freq
        self.best_gens_copy = None
        # after how many generations GAs will communicate
        self.comm_freq = comm_freq

        self.fitness_func = self.fitness_poly

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
        print('GenAlg #{} Warming up...'.format(self.gen_pos))
        self.fitness_func = self.fitness_mono
        for i in range(5*self.comm_freq):
            self._update_fitness_list()
            self._make_selection()
            self._crossover_populaton()
            self._mutate_popul()
        
        self.fitness_func = self.fitness_poly
        self.best_gens_push(self.popul[0])

    def ngram_index(self, groups, group_size):
        return sum(cc.ngram_index(group, group_size) for group in groups) / len(groups)
        
    def translate_cipher(self, indiv):
        self.best_gens_copy[self.gen_pos] = indiv

        key_maps = list()
        for i in range(self.gens_num):
            key_maps.append(str.maketrans(
                UPPERCASE_LETTERS,
                self.best_gens_copy[i].genes))

        maps_iter = itertools.cycle(key_maps)
        res_str = ''.join(
            (char.translate(next(maps_iter)) for char in self.ciphertext))

        return res_str

    def get_mono_text(self):
        start, end, step = self.gen_pos, len(self.ciphertext), self.gens_num
        text = ''.join(
            (self.ciphertext[i] for i in range(start, end, step)))
        return text

    def translate_mono_text(self, mono_text, key):
        key_map = str.maketrans(UPPERCASE_LETTERS, key)
        return mono_text.translate(key_map)

    def fitness_mono(self, indiv):
        mono_text = self.get_mono_text()
        trans_mono_text = self.translate_mono_text(mono_text, indiv.genes)
        index = cc.ngram_index(trans_mono_text, 1)
        fitness = abs(index - self.exp_mono_index)
        return fitness

    def fitness_poly(self, indiv):
        decoded = self.translate_cipher(indiv)
        index = cc.quadgram_index(decoded)
        fitness = abs(index - self.exp_index)
        return fitness

    def fitness(self, indiv):
        return self.fitness_func(indiv)

    def print_stats(self):
        print('GenAlg #{} Generation #{} Best fitness: {}'.format(
            self.gen_pos,
            self.generation, 
            self.fitness_list[0]))

    def solve(self, eps=0.17):
        self.best_gens_copy = self.best_gens_pull()
        self._update_fitness_list()

        while self.fitness_list[0] >= eps:
            if self.solved_lock.locked():
                return

            self._update_fitness_list()
            self._make_selection()
            self._crossover_populaton()
            self._mutate_popul()

            if self.generation % self.comm_freq == 0:
                self.print_stats()
                self.best_gens_push(self.popul[0])
                self.best_gens_copy = self.best_gens_pull()

            self.generation += 1

        self.solved_lock.acquire()
        self.best_gens_lock.acquire()
        self.best_gens = self.best_gens_copy[:]
        self.best_gens[self.gen_pos] = self.popul[0]
        self.best_gens_lock.release()
        print('GA #{} found solution with fitness {}'.format(
            self.gen_pos, self.fitness_list[0]))


class PolyalphabeticNaturalSelection():
    def __init__(self, 
            ciphertext, popul_num, alphabets_num, comm_freq, kill_percent=0.5, mutation_chance=0.1):
        self.ciphertext = ciphertext
        
        self.best_gens = [None] * alphabets_num
        self.best_gens_lock = threading.Lock()
        self.solved_lock = threading.Lock()

        self.genetic_algorithms = list()
        for i in range(alphabets_num):
            params = {
                'ciphertext': self.ciphertext,
                'gen_pos': i,
                'gens_num': alphabets_num, 
                'best_gens': self.best_gens,
                'best_gens_lock': self.best_gens_lock,
                'solved_lock': self.solved_lock,
                'comm_freq': comm_freq, 
                'popul_num': popul_num, 
                'kill_percent': kill_percent, 
                'mutation_chance': mutation_chance
            }
            self.genetic_algorithms.append(ParallelGeneticAlgorithm(**params))

    def translate(self, ciphertext, alphabets):
        key_maps = list()
        for alphabet in alphabets:
            key_maps.append(str.maketrans(UPPERCASE_LETTERS, alphabet))

        text = ''
        maps_iter = itertools.cycle(key_maps)
        
        key_map = next(maps_iter)
        for char in ciphertext:
            text += char.translate(key_map)
            key_map = next(maps_iter)

        return text

    def warm_up(self):
        threads = list()
        for ga in self.genetic_algorithms:
            threads.append(threading.Thread(
                target=ParallelGeneticAlgorithm.warm_up,
                args=(ga,)))
            threads[-1].start()

        for thread in threads:
            thread.join()

    def solve(self):
        threads = list()
        for ga in self.genetic_algorithms:
            threads.append(threading.Thread(
                target=ParallelGeneticAlgorithm.solve,
                args=(ga,)))
            threads[-1].start()

        for thread in threads:
            thread.join()

        best_ga = self.genetic_algorithms[0]
        for ga in self.genetic_algorithms[1:]:
            if ga.fitness_list[0] < best_ga.fitness_list[0]:
                best_ga = ga

        text = best_ga.translate_cipher(best_ga.popul[0])

        self.best_gens_lock.acquire()
        alphabets = [ga.genes for ga in self.best_gens]
        self.best_gens_lock.release()

        text = self.translate(self.ciphertext, alphabets)
        return text, alphabets


if __name__ == '__main__':
    CT = 'KZBWPFHRAFHMFSNYSMNOZYBYLLLYJFBGZYYYZYEKCJVSACAEFLMAJZQAZYHIJFUNHLCGCINWFIHHHTLNVZLSHSVOZDPYSMNYJXHMNODNHPATXFWGHZPGHCVRWYSNFUSPPETRJSIIZSAAOYLNEENGHYAMAZBYSMNSJRNGZGSEZLNGHTSTJMNSJRESFRPGQPSYFGSWZMBGQFBCCEZTTPOYNIVUJRVSZSCYSEYJWYHUJRVSZSCRNECPFHHZJBUHDHSNNZQKADMGFBPGBZUNVFIGNWLGCWSATVSSWWPGZHNETEBEJFBCZDPYJWOSFDVWOTANCZIHCYIMJSIGFQLYNZZSETSYSEUMHRLAAGSEFUSKBZUEJQVTDZVCFHLAAJSFJSCNFSJKCFBCFSPITQHZJLBMHECNHFHGNZIEWBLGNFMHNMHMFSVPVHSGGMBGCWSEZSZGSEPFQEIMQEZZJIOGPIOMNSSOFWSKCRLAAGSKNEAHBBSKKEVTZSSOHEUTTQYMCPHZJFHGPZQOZHLCFSVYNFYYSEZGNTVRAJVTEMPADZDSVHVYJWHGQFWKTSNYHTSZFYHMAEJMNLNGFQNFZWSKCCJHPEHZZSZGDZDSVHVYJWHGQFWKTSNYHTSZFYHMAEDNJZQAZSCHPYSKXLHMQZNKOIOKHYMKKEIKCGSGYBPHPECKCJJKNISTJJZMHTVRHQSGQMBWHTSPTHSNFQZKPRLYSZDYPEMGZILSDIOGGMNYZVSNHTAYGFBZZYJKQELSJXHGCJLSDTLNEHLYZHVRCJHZTYWAFGSHBZDTNRSESZVNJIVWFIVYSEJHFSLSHTLNQEIKQEASQJVYSEVYSEUYSMBWNSVYXEIKWYSYSEYKPESKNCGRHGSEZLNGHTSIZHSZZHCUJWARNEHZZIWHZDZMADNGPNSYFZUWZSLXJFBCGEANWHSYSEGGNIVPFLUGCEUWTENKCJNVTDPNXEIKWYSYSFHESFPAJSWGTYVSJIOKHRSKPEZMADLSDIVKKWSFHZBGEEATJLBOTDPMCPHHVZNYVZBGZSCHCEZZTWOOJMBYJSCYFRLSZSCYSEVYSEUNHZVHRFBCCZZYSEUGZDCGZDGMHDYNAFNZHTUGJJOEZBLYZDHYSHSGJMWZHWAFTIAAY'
    
    ns = PolyalphabeticNaturalSelection(CT, 100, 4, 10)
    print('+++++++++ Warming up +++++++++')
    ns.warm_up()
    print('+++++++++ Solving +++++++++')
    text = ns.solve()
    print(text)