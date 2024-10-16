# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import copy
import math
import sys
import argparse
import random
import matplotlib.pyplot as plt
import locAL
import randomDNA






if __name__ == '__main__':
    sys.setrecursionlimit(1500)
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=int, required=True)
    parser.add_argument('-s', type=int, required=True)
    parser.add_argument('-a', action='store_true')
    args = parser.parse_args()

    wdir = "/Users/Amit Elia/Documents/BENG182/datasets/"
    print(args)

    P1 = []
    for i in range(500):
        sequences = randomDNA.get_sequences(2, 1000)
        seq1 = sequences[0]
        seq2 = sequences[1]
        bt, S, tup, max_score = locAL.regular_local_alignment(seq1, seq2, args.m, args.s, 0)
        out = locAL.outputLocalAlignment(bt, S, seq1, seq2, tup[0], tup[1])
        length = locAL.format_alignment(out, args.m, args.s, 0)[1]
        P1 += [length]
        print(i)

    P2 = []
    for i in range(500):
        sequences = randomDNA.get_sequences(2, 1000)
        seq1 = sequences[0]
        seq2 = sequences[1]
        bt, S, tup, max_score = locAL.regular_local_alignment(seq1, seq2, args.m, args.s, -20)
        out = locAL.outputLocalAlignment(bt, S, seq1, seq2, tup[0], tup[1])
        length = locAL.format_alignment(out, args.m, args.s, -20)[1]
        P2 += [length]
        print(i)

    print(P1)
    print(P2)
    fig, axs = plt.subplots(1, 2, tight_layout=True)
    axs[0].hist(P1, alpha=0.5, label='P1', color="red")
    axs[0].legend(loc='upper right')
    axs[1].hist(P2, alpha=0.5, label='P2')
    axs[1].legend(loc='upper right')
    if(args.a):
        plt.show()








# See PyCharm help at https://www.jetbrains.com/help/pycharm/
