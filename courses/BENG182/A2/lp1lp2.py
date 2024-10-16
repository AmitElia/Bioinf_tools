# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import copy
import math
import sys
import argparse
import random
import matplotlib.pyplot as plt
import numpy
import randomDNA
import locAL






if __name__ == '__main__':
    sys.setrecursionlimit(1500)
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=int, required=True)
    parser.add_argument('-s', type=int, required=True)
    parser.add_argument('-a', action='store_true')
    args = parser.parse_args()

    wdir = "/Users/Amit Elia/Documents/BENG182/datasets/"
    print(args)

    lengths = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]
    lp1 = dict()
    lp2 = dict()
    for currlen in lengths:
        P1 = []
        for i in range(500):
            sequences = randomDNA.get_sequences(2, currlen)
            seq1 = sequences[0]
            seq2 = sequences[1]
            bt, S, tup, max_score = locAL.regular_local_alignment(seq1, seq2, args.m, args.s, 0)
            out = locAL.outputLocalAlignment(bt, S, seq1, seq2, tup[0], tup[1])
            length = locAL.format_alignment(out, args.m, args.s, 0)[1]
            P1 += [length]

        P2 = []
        for i in range(500):
            sequences = randomDNA.get_sequences(2, currlen)
            seq1 = sequences[0]
            seq2 = sequences[1]
            bt, S, tup, max_score = locAL.regular_local_alignment(seq1, seq2, args.m, args.s, -20)
            out = locAL.outputLocalAlignment(bt, S, seq1, seq2, tup[0], tup[1])
            length = locAL.format_alignment(out, args.m, args.s, -20)[1]
            P2 += [length]

        print(currlen)
        print(P1)
        print(P2)
        lp1[currlen] = numpy.mean(P1)
        lp2[currlen] = numpy.mean(P2)

    print(lp1)
    print(lp2)

    fig, axs = plt.subplots(1, 2, tight_layout=True)
    axs[0].plot(lengths, lp1.values(), '-o', alpha=0.5, label='P1', color="red")
    axs[0].legend(loc='upper right')
    axs[1].plot(lengths, lp2.values(), '-o', alpha=0.5, label='P2')
    axs[1].legend(loc='upper right')
    if(args.a):
        plt.show()

    '''
    lp1: {10: 10.046, 20: 23.668, 30: 37.118, 40: 50.618, 50: 64.016, 60: 77.78, 70: 91.45, 80: 104.854, 90: 118.134, 100: 131.758, 110: 145.374, 120: 159.066, 130: 172.668, 140: 186.388, 150: 199.802, 160: 213.392, 170: 226.558, 180: 240.378, 190: 253.74, 200: 267.388}
    lp2: {10: 2.578, 20: 3.7, 30: 4.324, 40: 4.814, 50: 5.154, 60: 5.454, 70: 5.718, 80: 5.912, 90: 6.146, 100: 6.2, 110: 6.386, 120: 6.508, 130: 6.628, 140: 6.778, 150: 6.9, 160: 6.982, 170: 7.05, 180: 7.146, 190: 7.22, 200: 7.282}
    '''







# See PyCharm help at https://www.jetbrains.com/help/pycharm/
