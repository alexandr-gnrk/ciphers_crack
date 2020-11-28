import itertools

import cipher_common as cc


def single_xor_cipher_decode(bytes_):
    res = bytes()
    min_chi2 = float('inf')
    key = None

    # for byte in bytes(string.printable, encoding='ascii'):
    for byte in range(256):
        # xor the string
        str_xor = cc.bytes_xor(bytes_, [byte])

        # count number of alphabet characters
        letters = cc.count_alphabet_bytes(str_xor)

        # get chi squared value
        chi2 = cc.chi_squared(
            cc.get_freq(str_xor.lower()), 
            letters, 
            cc.LETTER_FREQ)
        # update min chi2 and appropriate key
        if chi2 <= min_chi2:
            min_chi2 = chi2
            key = byte

    return cc.bytes_xor(bytes_, [key]), key