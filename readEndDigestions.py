import numpy as np
import re
# path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Accessibility/Data/Filetype Samples/sample.sam"
path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Accessibility/Data/Filetype Samples/sorted_4DNFIQVF4H3Z.sam"
motif = "GATC"

headerRows = 0
digestions = {}
with open(path) as f:
    for line in f:
        if not line.startswith("@"):
            break
        headerRows += 1
        if line.startswith("@SQ"):
            digestions[
                ''.join(line.split("SN:")[1].split("\t")[0])
            ] = np.empty(0, dtype=int)

cigar_delimiters = "M|I|D|N|S|H|P"
def cigar_shift(cigar, read_length):
    beforeM = re.split(cigar_delimiters, cigar) [:-2]
    return read_length - sum([int(i) for i in beforeM])

print("loading fields")

flags = np.loadtxt(path, usecols=1, dtype=int, skiprows=headerRows)
refs = np.loadtxt(path, usecols=2, dtype=str, skiprows=headerRows)
inds = np.loadtxt(path, usecols=3, dtype=int, skiprows=headerRows)
cigars = np.loadtxt(path, usecols=5, dtype=str, skiprows=headerRows)
seqs = np.loadtxt(path, usecols=9, dtype=str, skiprows=headerRows)

print("processing")

for i in range(len(flags)):
    reverse = not flags[i] < 16 and bin(flags[i])[-5] == "1"
    print(str(flags[i]) + "\t" + str(reverse))
    if (reverse and not seqs[i].endswith(motif)) or (not reverse and not seqs[i].startswith(motif)):
        continue

    ind = inds[i]
    if reverse:
        ind += cigar_shift(cigars[i], len(seqs[i]))

    digestions[refs[i]] = np.append(
        digestions[refs[i]],
        ind
    )

print(digestions)
# with f = open(path):
