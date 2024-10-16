import os
import numpy as np
import matplotlib.pyplot as plt

def PlotCovDist(plot=False):
    # Compute these values and return below
    mean_child_cov = 0
    mean_father_cov = 0
    mean_mother_cov = 0
    
    # Load coverage info
    # In trio.mpileup:
    #   4th column: child coverage
    #   7th column: father coverage
    #   10th column: mother coverage
    cov_child = []
    cov_father = []
    cov_mother = []
    with open("alignment/trio.mpileup") as f:
        for line in f:
            items = line.strip().split() # split columns of each line on whitespace
            cov_child.append(int(items[3]))
            cov_father.append(int(items[6]))
            cov_mother.append(int(items[9]))
    
    # Plot coverage distributions
    if plot:
        fig = plt.figure()
        fig.set_size_inches((8, 3))
        ax = fig.add_subplot(131)
        ax.hist(cov_child)
        ax.set_xlabel("# reads - child")
        ax.set_ylabel("# positions")
        ax = fig.add_subplot(132)
        ax.hist(cov_father)
        ax.set_xlabel("# reads - father")
        ax.set_ylabel("# positions")
        ax = fig.add_subplot(133)
        ax.hist(cov_mother)
        ax.set_xlabel("# reads - mother")
        ax.set_ylabel("# positions")    
        fig.tight_layout()
        plt.savefig('out/coverage_distribution_histogram.png')
    
    # Compute mean coverate of child, father, and mother.
    # Hint, see the example in plotting template 
    # length distribution at the end of Part 1.
    # your code here
    mean_child_cov = np.mean(cov_child)
    mean_father_cov = np.mean(cov_father)
    mean_mother_cov = np.mean(cov_mother)

    # Print out results
    print("Mean child coverage: %d"%mean_child_cov)
    print("Mean father coverage: %d"%mean_father_cov)
    print("Mean mother coverage: %d"%mean_mother_cov)

    # Return results
    return mean_child_cov, mean_father_cov, mean_mother_cov

PlotCovDist(plot=True);