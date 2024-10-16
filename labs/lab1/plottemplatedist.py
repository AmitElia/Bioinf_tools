import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def PlotTlenDist(plot=False):
    median_tlen = 0
    mean_tlen = 0
    
    # Read template lengths into a list (don't modify this path!!)
    tlens = open("alignment/child_template_lengths.txt").readlines()
    # Convert items in tlens list to integers
    tlens = [abs(int(item.strip())) for item in tlens]
    # Remove "0"
    tlens = [item for item in tlens if item > 0]
    
    # Compute median and mean of tlens below
    mean_tlen = np.mean(tlens)
    median_tlen = np.median(tlens)
    
    # Plot. See https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist.html
    print("Mean template length: %d"%mean_tlen)
    print("Median template length: %d"%median_tlen)

    mu, std = norm.fit(tlens)

    if plot:
        fig = plt.figure()
        fig.set_size_inches((8, 3))
        axs = fig.add_subplot(131)
        axs.hist(tlens, bins=25, color="g", density=True, alpha=0.6)
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        axs.plot(x, p, 'k', linewidth=2)
        title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
        axs.set_title(title)
        axs.set_xlabel("Fragment length", size=15)
        axs.set_ylabel("Number of fragments", size=15)
        axs.axvline(x=mean_tlen, color="red")
        plt.savefig('out/template_distribution_histogram.png')
        
        
    # Return stats    
    return median_tlen, mean_tlen

PlotTlenDist(plot=True)