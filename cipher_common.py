import itertools
import collections
import string
import math

LETTER_FREQ = {
    'a': 0.0855, 'b': 0.0160, 'c': 0.0316, 'd': 0.0387, 'e': 0.1210,
    'f': 0.0218, 'g': 0.0209, 'h': 0.0496, 'i': 0.0733, 'j': 0.0022,
    'k': 0.0081, 'l': 0.0421, 'm': 0.0253, 'n': 0.0717, 'o': 0.0747,
    'p': 0.0207, 'q': 0.0010, 'r': 0.0633, 's': 0.0673, 't': 0.0894,
    'u': 0.0268, 'v': 0.0106, 'w': 0.0183, 'x': 0.0019, 'y': 0.0171, 
    'z': 0.0011
}


def __read_ngrams(filename):
    quadgrams = dict()
    for line in open(filename):
        quadgram, number = line.split()
        quadgrams[quadgram.lower()] = int(number)

    summa = sum(quadgrams.values())
    for key, val in quadgrams.items():
        quadgrams[key] = math.log10(val / summa)

    # to floor other quadgrams
    quadgrams['floor'] = math.log10(0.01 / summa)
    return quadgrams


NGRAMS_LIST = (
    __read_ngrams('./data/english_monograms.txt'),
    __read_ngrams('./data/english_bigrams.txt'),
    __read_ngrams('./data/english_trigrams.txt'),
    __read_ngrams('./data/english_quadgrams.txt'))


def chi_squared(obs_freq, letters, exp_freq):
    chi2 = 0

    # return infinity if there aren't any alphabet character
    if letters == 0:
        return float('inf')

    for obs_key, obs_val in obs_freq.items():
        # if there is a non printable character, return infinity
        if (chr(obs_key) not in string.printable) and \
                (bytes([obs_key]) not in (b'\x80', b'\x99', b'\xe2')):
            return float('inf')
        elif chr(obs_key) not in exp_freq:
            continue
        else:
            exp_val = letters * exp_freq[chr(obs_key)]
            chi2 += (obs_val - exp_val)**2 / exp_val

    return chi2


def get_freq(bytes_):
    return dict(collections.Counter(bytes_))


def count_alphabet_bytes(bytes_):
    res = len(list(itertools.filterfalse(
        lambda byte: not bytes([byte]).isalpha(), 
        bytes_)))
    return res


def bytes_xor(bytes_, key):
    res = bytes(x ^ y for (x,y) in zip(bytes_, itertools.cycle(key)))
    return res


def coincidence_index(text):
    ci = 0.0
    text_len = len(text)
    text_freq = get_freq(text)

    for amount in text_freq.values():
        ci += amount**2

    ci /= text_len**2
    return ci


def ngram_index(text, group_size):
    ngrams = NGRAMS_LIST[group_size - 1] 
    text = text.lower()
    index = 0
    for i in range(len(text) - group_size + 1):
        ngram = text[i:i + group_size]
        index += ngrams.get(ngram, ngrams['floor'])
    return index / (len(text) - group_size + 1)


def monogram_index(text):
    return ngram_index(text, 1)

def bigram_index(text):
    return ngram_index(text, 2)

def trigram_index(text):
    return ngram_index(text, 3)

def quadgram_index(text):
    return ngram_index(text, 4)