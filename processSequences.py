import numpy as np


# WARNING: Unsorted
def find_sub_ends(str, sub):
    """Returns a boolean array containing true at the start and end indices of each instance of sub in str.
    NOTE: instances of sub which start inside a previous instance will not be counted"""
    is_end = np.full(len(str), False)
    # This feels like an awful way of doing this but of all the things I've tried it's the fastest by far
    start = 0
    while True:
        start = str.find(sub, start)
        if start == -1: return is_end
        is_end[start] = True
        start += len(sub)
        is_end[start - 1] = True


def site_instances(seq, site_seq):
    """Returns a boolean array containing true at the start and end positions of each restriction site"""
    forward = find_sub_ends(seq, site_seq)
    backward = find_sub_ends(seq, site_seq[::-1])
    is_site = forward | backward
    return is_site, forward, backward


def load_ref_by_id(fasta_iterator, id):
    for fasta in fasta_iterator:
        if fasta.id == id:
            return fasta