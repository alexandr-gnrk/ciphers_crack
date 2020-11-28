import cipher_common as cc
import single_xor_cipher as sxc


def get_vigenere_cipher_key_len(bytes_, eps=0.01):
    # expected coincidence index for English language
    EXP_IND = 0.065

    for key_len in range(1, len(bytes_)):
        # make slice
        sli = bytes_[::key_len]
        ind = cc.coincidence_index(sli)
        # if accuracy fits, return key_len
        if abs(EXP_IND - ind) <= eps:
            return key_len

    # otherwise return len of bytes_
    return len(bytes_)


def vigenere_cipher_decode(bytes_):
    key_len = get_vigenere_cipher_key_len(bytes_)
    key = str()
    parts = list()
    length = 0

    for i in range(key_len):
        part, key_char = sxc.single_xor_cipher_decode(bytes_[i::key_len])
        key += chr(key_char)
        parts.append(part)
        length += len(part)

    # concat all slices in one list
    res = bytes()
    for char_num in range(length):
        part_num = char_num % key_len
        byte_num = char_num // key_len
        res += bytes([parts[part_num][byte_num]])

    return res, key


import itertools
def vigenere_cipher_decode_by_key(bytes_, key):
    key_byte = itertools.cycle(key)
    res = bytes()
    for byte in bytes_:
        res += cc.bytes_xor(bytes([byte]), bytes([next(key_byte)]))

    return res