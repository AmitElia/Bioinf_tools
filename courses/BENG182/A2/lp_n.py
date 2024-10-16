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
    parser.add_argument('-d', type=int, required=True)
    parser.add_argument('-a', action='store_true')
    args = parser.parse_args()

    wdir = "/Users/Amit Elia/Documents/BENG182/datasets/"
    print(args)

    # mismatch = indel
    scores = [-4,-3,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0]
    #scores = [-30, -25, -20, -15, -10, -5, -2, -1.5, -1, -0.5, -0.4, -0.3, -0.2, -0.1 , 0]
    lp = dict()
    for score in scores:
        P = []
        for i in range(200):
            sequences = randomDNA.get_sequences(2, 1000)
            seq1 = sequences[0]
            seq2 = sequences[1]
            bt, S, tup, max_score = locAL.regular_local_alignment(seq1, seq2, args.m, -math.inf, score)
            #out, length = locAL.locAL(seq1, seq2, args.m, score, score)
            out = locAL.outputLocalAlignment(bt, S, seq1, seq2, tup[0], tup[1])
            length = locAL.format_alignment(out, args.m, score, score)[1]
            P += [length]


        print(score, P)
        lp[score] = numpy.mean(P)

    print(lp)

    plt.plot(scores, lp.values(), '-o', alpha=0.5, label='P1', color="red")
    plt.legend(loc='upper right')
    if(args.a):
        plt.show()

    '''
    Run 1: [-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,-1.5,-1,-0.5,-0.25,-0.125,0]
    {-30: 9.71, -25: 9.64, -20: 9.65, -15: 9.59, -10: 9.65, -5: 9.84, -2: 16.29, -1.5: 44.375, -1: 1065.59, -0.5: 1165.44, -0.4: 1179.5, -0.3: 1188.14, -0.2: 1198.37, -0.1: 1205.34, 0: 1229.72}
    
    Run 2: [-4,-3,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0]
    {-4: 9.895, -3: 10.495, -2.5: 12.58, -2.25: 14.625, -2: 16.245, -1.75: 24.04, -1.5: 49.705, -1.25: 302.83, -1: 1064.43, -0.75: 1141.93, -0.5: 1166.175, -0.25: 1192.935, 0: 1229.3}
    
    Final3: 
    
    {-4: 9.905, -3: 10.61, -2.5: 11.39, -2.25: 11.765, -2: 13.05, -1.75: 16.495, -1.5: 21.585, -1.25: 36.51, -1: 209.51, -0.75: 1247.695, -0.5: 1316.05, -0.25: 1339.64, 0: 1347.455}
    '''







# See PyCharm help at https://www.jetbrains.com/help/pycharm/
