from processSequences import *
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from insights import *
from Bio import SeqIO
from scipy.optimize import curve_fit

# Get actual digestions from Hi-C Experiment
path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Accessibility/Data/4DNESPXW8XHY/4DNFIOJB13JQ_chr1_digestions.txt"
print("Loading file...")
digestions = np.loadtxt(path, dtype=int)

print("Counting digestions...")
digestions_counter = Counter(digestions)
digestions_binned = np.array(sorted(digestions_counter.items()))
digestion_inds = digestions_binned[:, 0]
digestion_freqs = digestions_binned[:, 1]

digestion_count = len(digestions)
print("Digestions: " + str(digestion_count))
digested_site_count = len(digestions_counter.keys())
print("Digested site count: " + str(digested_site_count))

# Get the indices of restriction sites in the reference genome
path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Reference Genomes/ncbi-genomes-2023-03-14/GCA_000001405.15_GRCh38_genomic.fna"
print("Loading reference genome...")
fasta_iterator = SeqIO.parse(open(path), 'fasta')
reference = str(load_ref_by_id(fasta_iterator, 'CM000663.2').seq).upper()

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
print(digestions_not_in_ref[1:10])
print([find_nearest(site_inds, x) for x in digestions_not_in_ref[1:10]])


pct_odd_digestions = 100 * round(len(digestions_not_in_ref) / digestion_count, 2)
print("Percentage of digestions that did not appear in the reference genome: " + str(pct_odd_digestions))
undigested_site_percentage = 100 * round(len(undigested_sites) / site_count, 2)
print("Percentage of sites in the reference genome that were not digested: " + str(undigested_site_percentage))

counts, bins = np.histogram(digestion_freqs, bins=max(digestion_freqs) - min(digestion_freqs))


def exp3(x, a, b, c, d, e, f):
    return a * np.exp(-b * x) + c * np.exp(-d * x) + e * np.exp(-f * x)


fig, [ax, ax2] = plt.subplots(2)

ax.stairs(counts, bins)
ax.set_xlabel("Number of digestions")
ax.set_ylabel("Instances")
ax.set_yscale("log")
ax.set_xscale("log")

ax2.stairs(counts, bins)
ax2.set_xlabel("Number of digestions")
ax2.set_ylabel("Instances")
ax2.set_yscale("log")

x = bins[:-1]
popt, pcov = curve_fit(exp3, x, counts, maxfev=6000000)
y = exp3(x, *popt)
print(popt)
ax2.plot(x, y, color="navy", linestyle="dashed")
ax2.set_ylim([1, max(counts)])

plt.show()