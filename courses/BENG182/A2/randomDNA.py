import copy
import math
import sys
import argparse
import random



def random_seq(length):
    seq = ""
    for i in range(length):
        seq += random.choice("ACGT")
    return seq

def get_sequences(num, length):
    sequences = []
    for i in range(num):
        sequences += [random_seq(length)]

    return sequences

def get_summary(seqs):
    freqs = {'A':0,'C':0,'G':0,'T':0}
    count = len(seqs[0])*len(seqs)
    for seq in seqs:
        for char in seq:
            freqs[char] += 1
    for nuc in freqs.keys():
        freqs[nuc] /= count

    return freqs

if __name__ == '__main__':
    sys.setrecursionlimit(1500)
    parser = argparse.ArgumentParser()
    parser.add_argument('num', type=int)
    parser.add_argument('len', type=int)
    args = parser.parse_args()

    wdir = "/Users/Amit Elia/Documents/BENG182/datasets/"
    print(args)
    seqs = get_sequences(args.num, args.len)
    summ = get_summary(seqs)
    for seq in seqs:
        print(seq)
    #print(seqs)
    print(summ)