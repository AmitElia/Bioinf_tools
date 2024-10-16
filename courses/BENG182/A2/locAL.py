# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import copy
import math
import sys
import argparse



# Press the green button in the gutter to run the script.

def read_fasta_file(file_name):
    file = open(file_name, "r")
    seqs = dict()
    name = ""
    for line in file:
        if line[0] == '>':
            name = line.strip()[1:]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs

def locAL2(seq1, seq2, m, s, d):

    S = linear_space_alignment(seq1, seq2, 0, len(seq1), 0, len(seq2), m, s, d)
    out, length = format_alignment(S, m, s, d)
    return out, length



def regular_local_alignment(seq1, seq2, m, s, d):



    row = [0] * (len(seq2) + 1)
    backtrack_row = ['E'] * (len(seq2) + 1)
    backtrack_row[0] = 'S'
    score_matrix = []
    backtrack_matrix = []
    for i in range(len(seq1) + 1):
        score_matrix.append(copy.deepcopy(row))
        backtrack_matrix.append(copy.deepcopy(backtrack_row))
        score_matrix[i][0] = 0
        #score_matrix[i][0] = d * i
    for j in range(len(seq2) + 1):
        score_matrix[0][j] = 0
        #score_matrix[0][j] = d * j
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) +1):
            match = s
            if seq1[i-1] == seq2[j-1]:
                match = m
            score_matrix[i][j] = max(0,score_matrix[i-1][j-1] + match, score_matrix[i-1][j] + d, score_matrix[i][j-1] + d)

            if score_matrix[i][j] == score_matrix[i-1][j-1] + match:
                backtrack_matrix[i][j] = 'D'
            elif score_matrix[i][j] == score_matrix[i-1][j] + d:
                backtrack_matrix[i][j] = 'S'
            elif score_matrix[i][j] == score_matrix[i][j-1] + d:
                backtrack_matrix[i][j] = 'E'

    max_score = 0
    tup = (len(seq1) +1, len(seq2) + 1)
    for i in range(len(seq1) + 1):
        for j in range(len(seq2) +1):
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                tup = (i, j)

    return backtrack_matrix, score_matrix, tup, max_score

def linear_local_alignment(seq1, seq2, m, s, d):

    row = [0] * (len(seq2) + 1)
    #backtrack_row = ['N'] * (len(seq2) + 1)
    #backtrack_row[0] = 'N'
    bt = dict()
    score_matrix = []
    for i in range(2):
        score_matrix.append(copy.deepcopy(row))
    #score_matrix = 2* [row]
    #S = []
    #backtrack_matrix = []
    for i in range(len(seq1) + 1):
        #S.append(copy.deepcopy(row))
        #backtrack_matrix.append(copy.deepcopy(backtrack_row))
        bt[i] = dict()
        if i == 0:
            for j in range(len(seq2) + 1):
                bt[i][j] = 'N'
        else:
            bt[i][0] = 'N'
        #S[i][0] = 0
    #print(score_matrix, seq1, seq2)

    max_score = 0
    tup = (len(seq1) + 1, len(seq2) + 1)
    for i in range(1, len(seq1) + 1):
        #print(i)
        for j in range(1, len(seq2) +1):
            i_2 = i%2
            i_1 = (i-1)%2
            #print(i, "***", i_2, i_1)
            match = s
            if seq1[i-1] == seq2[j-1]:
                match = m
            #S[i][j] = max(0,S[i-1][j-1] + match, S[i-1][j] + d, S[i][j-1] + d)
            score_matrix[i_2][j] = max(0, score_matrix[i_1][j - 1] + match, score_matrix[i_1][j] + d,
                                     score_matrix[i_2][j - 1] + d)
            #print(score_matrix[i_2][j])

            if score_matrix[i_2][j] > max_score:
                max_score = score_matrix[i_2][j]
                tup = (i, j)
                #print(tup)
            if score_matrix[i_2][j] == 0:
                bt[i][j] = 'N'
            if score_matrix[i_2][j] == score_matrix[i_1][j - 1] + match:
                #backtrack_matrix[i][j] = 'D'
                bt[i][j] = 'D'
            elif score_matrix[i_2][j] == score_matrix[i_1][j] + d:
                #backtrack_matrix[i][j] = 'S'
                bt[i][j] = 'S'
            elif score_matrix[i_2][j] == score_matrix[i_2][j - 1] + d:
                #backtrack_matrix[i][j] = 'E'
                bt[i][j] = 'E'

        #print(score_matrix)

    #for line in bt.keys():
    #    print(bt[line])
    #for line in backtrack_matrix:
    #    print(line)

    print(tup)
    #bt = bt[:tup[0] + 1][:tup[1] + 1]
    return bt, tup, max_score

def outputLinearLocalAlignment(backtrack_matrix, seq1, seq2, i, j):

    if backtrack_matrix[i][j] == 'S':
        #print(backtrack_matrix[i][j] ,i, j)
        out = outputLinearLocalAlignment(backtrack_matrix, seq1, seq2, i-1, j)
        out.append((seq1[i-1],'-'))

        return out
    elif backtrack_matrix[i][j] == 'E':
        #print(backtrack_matrix[i][j] ,i, j)
        out = outputLinearLocalAlignment(backtrack_matrix, seq1, seq2, i, j-1)
        out.append(('-',seq2[j-1]))

        return out
    elif backtrack_matrix[i][j] == 'D':
        #print(backtrack_matrix[i][j] ,i, j)
        out = outputLinearLocalAlignment(backtrack_matrix, seq1, seq2, i-1, j-1)
        out.append((seq1[i-1],seq2[j-1]))
        return out
    else:
        #print(backtrack_matrix[i][j] ,i, j)
        return []

def outputLocalAlignment(backtrack_matrix, score_matrix, seq1, seq2, i, j):
    if score_matrix[i][j] == 0:
        return []
    if backtrack_matrix[i][j] == 'S':
        out = outputLocalAlignment(backtrack_matrix, score_matrix, seq1, seq2, i-1, j)
        out.append((seq1[i-1],'-'))
        return out
    elif backtrack_matrix[i][j] == 'E':
        out = outputLocalAlignment(backtrack_matrix, score_matrix, seq1, seq2, i, j-1)
        out.append(('-',seq2[j-1]))
        return out
    else:
        out = outputLocalAlignment(backtrack_matrix, score_matrix, seq1, seq2, i-1, j-1)
        out.append((seq1[i-1],seq2[j-1]))
        return out

def format_alignment(tuple_list, m, s, d):
    a = ""
    b = ""
    mid = ""
    score = 0
    length = len(tuple_list)
    while len(tuple_list):
        tupl = tuple_list.pop(0)
        if(tupl[0] == '-') or (tupl[1] == '-'):
            score += d
            mid += " "
        else:
            match = s
            if tupl[0] == tupl[1]:
                match = m
                mid += "|"
            else:
                mid += " "
            score += match
        a += tupl[0]
        b += tupl[1]
    out ="Score:" + str(score) + '\n' + "Length:" +  str(length) +'\n' + a + '\n'+ mid + '\n' + b
    return out, length


def get_mid_col(seq1, seq2, m, s, d):
    prev = []
    bt = []
    S = []
    for j in range(len(seq2) + 1):
        S = []
        bt = []
        for i in range(len(seq1) + 1):
            if j == 0:
                #S.insert(i, d*i)
                S.insert(i, 0)
                bt.insert(i, "S")
            elif i == 0:
                #S.insert(i,d *j)
                S.insert(i, 0)
                bt.insert(i, "E")
            else:
                match = s
                if seq1[i - 1] == seq2[j - 1]:
                    match = m
                decision = max(prev[i-1] + match, S[i-1] + d, prev[i] + d)
                S.insert(i, decision)
                if decision == prev[i-1] + match:
                    bt.insert(i, "D")
                elif decision == prev[i] + d:
                    bt.insert(i, "E")
                elif decision == S[i-1] + d:
                    bt.insert(i, "S")
        prev = copy.deepcopy(S)
    return S, bt

def middleEdge(seq1, seq2, top, bottom, left, right, m, s, d):
    Lengths = []
    midJ = int((left + right)/2)
    fromSource = get_mid_col(seq1[top:bottom],seq2[left:midJ], m, s, d)
    toSink = get_mid_col(seq1[top:bottom][::-1],seq2[midJ:right][::-1], m, s, d)

    max_length = -math.inf
    imax = -1
    for i in range(bottom - top + 1):
        Lengths += [fromSource[0][i] + toSink[0][::-1][i]]
        if Lengths[i] > max_length:
            max_length = Lengths[i]
            imax = i

    #The node in the middle column with maximum score
    max_node = (imax + top, midJ)

    #the node right after max_node in the optimal path
    next_node = max_node

    #the edge from max_node to next_node
    max_edge = toSink[1][::-1][imax]
    if max_edge == "D":
        next_node = (top + imax + 1, midJ + 1)
    elif max_edge == "E":
        next_node = (top + imax, midJ +1)
    elif max_edge == "S":
        next_node = (top + imax + 1, midJ)
    return max_node, max_edge, next_node


def linear_space_alignment(seq1, seq2, top, bottom, left, right, m, s, d):
    S = []

    if left == right:
        for i in range(bottom - top):
            S += [(seq1[i + top],"-")]
        return S

    if bottom == top:
        for j in range(right - left):
            S += [("-",seq2[j + left])]
        return S

    middle = int((left + right)/2)
    mid_edge = middleEdge(seq1, seq2, top, bottom, left, right, m, s, d)
    mid_node_row = mid_edge[0][0]

    S += linear_space_alignment(seq1, seq2, top, mid_node_row, left, middle, m, s, d)

    if mid_edge[1] == "D":
        S += [(seq1[mid_node_row], seq2[middle])]
        middle += 1
        mid_node_row += 1
    elif mid_edge[1] == "E":
        S += [("-",seq2[middle])]
        middle += 1
    elif mid_edge[1] == "S":
        S += [(seq1[mid_node_row],"-")]
        mid_node_row += 1

    S += linear_space_alignment(seq1, seq2, mid_node_row, bottom, middle, right, m, s, d)

    return S

def format_alignment2(S, m, s, d):
    a = ""
    b = ""
    score = 0

    while(len(S)):
        tuple = S.pop(0)

        if (tuple[0] == "-") or (tuple[1] == "-"):
            score += d
        else:
            match = s
            if tuple[0] == tuple[1]:
                match = m
            score += match
        a += tuple[0]
        b += tuple[1]

    length = len(S)
    output = str(length) + '\n' + str(score) + '\n' + a + '\n' + b

    return output


if __name__ == '__main__':
    sys.setrecursionlimit(1500)

    parser = argparse.ArgumentParser()
    parser.add_argument('seq_files', type=str)
    parser.add_argument('-m', type=int, required=True)
    parser.add_argument('-s', type=int, required=True)
    parser.add_argument('-d', type=int, required=True)
    parser.add_argument('-a', action='store_true')
    args = parser.parse_args()

    wdir = "/Users/Amit Elia/Documents/BENG182/datasets/"
    print(args)
    fasta_file = wdir + args.seq_files
    sequences = read_fasta_file(fasta_file)
    seq1 = list(sequences.values())[0]
    seq2 = list(sequences.values())[1]
    #print(seq1)
    #print(seq2)
    #bt, S, tup, max_score = regular_local_alignment(seq1, seq2, args.m, args.s, args.d)
    bt, tup, max_score = linear_local_alignment(seq1, seq2, args.m, args.s, args.d)
    out = outputLinearLocalAlignment(bt, seq1, seq2, tup[0], tup[1])
    """
    for i in range(len(bt)):
        print(bt[i])
    for i in range(len(S)):
        print(S[i])
    """
    output,length = format_alignment(out, args.m, args.s, args.d)
    if(args.a == True):
        print(output)





# See PyCharm help at https://www.jetbrains.com/help/pycharm/
