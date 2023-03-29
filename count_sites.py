import numpy as np


def read_end_frequency(read_end_inds, seq_length):
    freq = np.zeros(seq_length)
    for ind in read_end_inds:
        freq[ind] += 1

    return freq


def digestion_frequency(read_end_freq, issite):
    return np.multiply(read_end_freq, issite)