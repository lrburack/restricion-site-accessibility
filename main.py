from processSequences import *
from count_sites import *
import matplotlib.pyplot as plt
from insights import *
import numpy as np
from Bio import SeqIO

# Load the file containing read end positions (with duplicates)
path = "Raw Hi-C/Processing Bam Files/chr1_plus.txt"
read_end_positions = np.loadtxt(path, dtype=int)

# Load the file containing digestion positions
path = "Raw Hi-C/Processing Bam Files/chr1_plus_digestedOnly.txt"
digestion_positions = np.loadtxt(path, dtype=int)

bincount = 300

fig, [ax1, ax2, ax3] = plt.subplots(3)
fig.tight_layout(pad=2)

site_distribution(ax1, read_end_positions, bincount=bincount)
ax1.set_title("Read End Positions")

site_distribution(ax2, digestion_positions, bincount=bincount)
ax2.set_title("Digestion positions")

# Normalize both arrays, divide one by the other, and plot the result
ax3.set_title("Normalized")
read_ends_binned, bins = np.histogram(read_end_positions, bins=bincount)
digestions_binned, bins = np.histogram(digestion_positions, bins=bincount)
y = np.divide(np.divide(read_ends_binned, np.max(read_ends_binned)), np.divide(digestions_binned, np.max(digestions_binned)))
ax3.plot(bins[1:], y)
ax3.set_ylim([0, 5])
ax3.plot([bins[0], bins[len(bins)-1]], [1, 1], linestyle="dashed", color="black")

print("Read end count: " + str(len(read_end_positions)))
print("Digestion count: " + str(len(digestion_positions)))

frequency = read_end_frequency(digestion_positions, max(read_end_positions))
print("Max: " + str(np.max(frequency)))

# Load the reference genome
reference_path = "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Reference Genomes/ncbi-genomes-2023-03-14/GCA_000001405.15_GRCh38_genomic.fna"

print("Loading reference genome")
fasta_iterator = SeqIO.parse(open(reference_path), 'fasta')
reference = str(load_ref_by_id(fasta_iterator, 'CM000663.2').seq)

print("Locating restriction sites")
restriction_seq = "GATC"
site_mask = site_instances(reference, restriction_seq)

read_end_freq = read_end_frequency(read_end_positions, len(reference))
freq = digestion_frequency(read_end_freq, site_mask)

print("Number of Digestions: " + str(np.sum(freq)))
print(len(digestion_positions))

# # Create a mask of restriction site positions
#

# Get the frequency that each position in the reference genome was a read end
#

# # Filter only for read ends at restriction sites
# digestion_freq = digestion_frequency(read_end_freq, site_mask)