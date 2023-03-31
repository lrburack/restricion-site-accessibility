from processSequences import *
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from Bio import SeqIO
from scipy.optimize import curve_fit

# Get actual digestions from Hi-C Experiment
path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Accessibility/Data/4DNESPXW8XHY/4DNFIOJB13JQ_chr7_digestions.txt"
digestions = np.loadtxt(path, dtype=int)

digestions_counter = Counter(digestions)
digestions_binned = np.array(sorted(digestions_counter.items()))
digestion_inds = digestions_binned[:, 0]
digestion_freqs = digestions_binned[:, 1]

print(digestion_inds[1:100])
print(digestion_freqs[1:100])

digestion_count = len(digestions)
print("Digestions: " + str(digestion_count))
digested_site_count = len(digestions_counter.keys())
print("Digested site count: " + str(digested_site_count))

# Get the indices of restriction sites in the reference genome
path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Reference Genomes/ncbi-genomes-2023-03-14/GCA_000001405.15_GRCh38_genomic.fna"
print("Loading reference genome...")
fasta_iterator = SeqIO.parse(open(path), 'fasta')
reference = str(load_ref_by_id(fasta_iterator, 'CM000663.2').seq)

print("Locating restriction sites")
restriction_seq = "GATC"
site_inds = site_instances(reference, restriction_seq)

site_count = int(len(site_inds))
print("Sites in reference genome: " + str(site_count))
sites_digested_percentage = round(100 * digested_site_count / site_count, 2)
print("Percentage of total sites digested: " + str(sites_digested_percentage) + "%")
digestions_per_site = round(digestion_count / site_count, 2)
print("Mean digestions per site: " + str(digestions_per_site))
digestions_per_accessible = round(digestion_count / digested_site_count, 2)
print("Mean digestions for sites digested at least once: " + str(digestions_per_accessible))

digestions_not_in_ref = np.setdiff1d(digestion_inds - 1, site_inds)
undigested_sites = np.setdiff1d(site_inds, digestion_inds - 1)

pct_odd_digestions = 100 * round(len(digestions_not_in_ref) / digestion_count, 2)
print("Percentage of digestions that did not appear in the reference genome: " + str(pct_odd_digestions))
undigested_site_percentage = 100 * round(len(undigested_sites) / site_count, 2)
print("Percentage of sites in the reference genome that were not digested: " + str(undigested_site_percentage))

# mask = digestion_freqs > 50
# print(mask)
# filtered_digestion_freqs = digestion_freqs[not mask]
counts, bins = np.histogram(digestion_freqs, bins=max(digestion_freqs) - min(digestion_freqs))

# counts = np.insert(counts, 0, len(undigested_sites))
# bins = np.insert(bins, 0, 0)


print(len(undigested_sites))
print(bins[0:5])
print(counts[0:5])

fig, ax = plt.subplots(1)

# ax.bar(range(1,100), digestion_freqs[1:100])
ax.stairs(counts, bins)

# def exponential(x, alpha, beta):
#     return beta * np.exp(-x/alpha)
#
#
# x = bins[1:]
# popt, pcov = curve_fit(exponential, x, counts, maxfev=6000000)
# y = exponential(x, *popt)
# ax.plot(x, y, color="navy", linestyle="dashed")

# ax.set_xscale("log")
ax.set_yscale("log")

plt.show()
# print(site_inds[1:100])
# print(list(digestions_counter.keys())[1:100])
print(digestion_freqs[1:100])


