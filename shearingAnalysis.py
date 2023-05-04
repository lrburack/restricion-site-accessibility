import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Accessibility/Data/4DNESQWI9K2F/readEnds_plus_4DNFIHBDN5SX_chr7.txt"
print("Reading file")
read_end_seqs = np.loadtxt(path, dtype=str)
print("Total read ends: " + str(len(read_end_seqs)))

print("Counting ends")
counter = Counter(read_end_seqs)

unique = np.array(list(counter.keys())) # equals to list(set(words))
freq = np.array(list(counter.values())) # counts the elements' frequency
keepend = np.array(["N" not in val for val in unique])
unique = unique[keepend]
freq = freq[keepend]

unique = [x for _, x in sorted(zip(freq, unique))]
freq = np.sort(freq) / np.sum(freq)

print(str(len(freq)) + " unique sequences")

for i in range(len(freq)):
    print(unique[i] + ": " + str(freq[i]))

counts, bins = np.histogram(freq, bins=100)

fig, [ax1, ax2] = plt.subplots(2)
ax1.stairs(counts, bins)

ax2.bar(unique, freq)
plt.xticks(rotation='vertical')

plt.show()

print("Total read ends after filtering: " + str(sum(list(counter.values()))))
print("Mean (w/o GATC): " + str(np.mean(freq[:-1])))
print("Std (w/o GATC) " + str(np.std(freq[:-1])))
print("GATC read end count: " + str(counter['GATC']))
print("CTAG read end count (use as estimate for GATC shearing): " + str(counter['CTAG']))
print("Estimated percentage of GATC at read end due to shearing: " + str(100 * counter['CTAG'] / counter['GATC']) + "%")
print("Percentage of read ends with a GATC: " + str(100 * counter['GATC'] / len(read_end_seqs)) + "%")