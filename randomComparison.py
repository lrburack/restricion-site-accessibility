import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.stats import binom

dummy_sites = np.zeros(100000)
p = 1 / len(dummy_sites)
print(len(dummy_sites))
iterations = 1000000

for i in range(0, iterations):
    dummy_sites[random.randrange(0, len(dummy_sites))] += 1


counts, bins = np.histogram(dummy_sites, bins=int(np.max(dummy_sites)) - int(np.min(dummy_sites)))
counts = counts / np.sum(counts)
fig, ax = plt.subplots(1)

# ax.bar(range(1,100), digestion_freqs[1:100])
ax.plot(bins[1:], counts)
# ax.set_xscale("log")
# ax.set_yscale("log")
x = range(0, 100)
ax.set_xlabel("Number of Digestions")
ax.set_ylabel("Frequency")

# ax.plot(x, binom.pmf(x, iterations, p), color='black')

plt.show()

