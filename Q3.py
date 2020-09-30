"""
Question 3 (10 points): Generate 100 random 1KB sequences (as if they are upstream of
some genes, you can use real data from UCSC for some co-regulated genes as well). Plant a
motif of length 10 after introducing 0,1,2 random mutations (uniformly distributed). Use
Gibb’s sampler to identify the motif locations and the consensus motif.
"""


import random
import math
from numpy.random import choice


def gibbs(seq_with_motif, len_seq_with_motif, len_of_motif):
    """
        1) Choose 1 sequence at random and remove it from the list of sequences
        2) Create profile matrix picking random k-mers from remaining sequences (probability matrix, so divide all numbers
                                                                            by 4 so all numbers in a column add up to 1)
        3) Find consensus string
        4) Scan every k-mer from the sequence selected in step 1 and find probability (by multiplying the probability of
                                                                                            nucleotides being the same)
        5) Normalize all the probabilities by dividing them by the minimum probability > 0
        6) Divide all these probabilities by the sum of probabilities
        7) Select one of these k-mers probabilistically (random function)

        Keep following 1-7 till the score can't be increased anymore
    """
    last_score = -1
    this_score = 0
    separate_seq = None

    ans_lmer = ""
    loop_variable = 0
    motif_pos_to_return = []
    for i in range(len(seq_with_motif)):
        motif_pos_to_return.append(random.randint(0, len_seq_with_motif - len_of_motif))

    while this_score > last_score:
        loop_variable += 1
        last_score = max(this_score, last_score)
        prob_for_this_run = dict()

        sequence_number_to_not_include = random.randint(0, len(seq_with_motif) - 1)
        separate_seq = seq_with_motif[sequence_number_to_not_include]  # because we'll pick a diff seq every time

        supposed_motifs = []  # lmers
        for i in range(len(seq_with_motif)):
            sequence = seq_with_motif[i]
            if sequence != separate_seq:
                supposed_motifs.append(sequence[motif_pos_to_return[i]: motif_pos_to_return[i] + len_of_motif])

        consensus_str, profile_mat = find_consensus(set(supposed_motifs), len(supposed_motifs))
        # Make probability matrix
        for el in profile_mat["A"]:
            el = el/4
        for el in profile_mat["T"]:
            el = el/4
        for el in profile_mat["C"]:
            el = el/4
        for el in profile_mat["G"]:
            el = el/4

        for i in range(0, len_seq_with_motif - len_of_motif):
            lmer_in_sep_seq = separate_seq[i:i+len_of_motif]
            probability = 1
            """find probability step 4"""
            for pos in range(len(lmer_in_sep_seq)):
                probability *= profile_mat[lmer_in_sep_seq[pos]][pos] # probability of this lmer to be motif
            prob_for_this_run[lmer_in_sep_seq] = probability
            """Normalize probability"""
            divisor = sum(prob_for_this_run.values())
            for key in prob_for_this_run.keys():
                prob_for_this_run[key] /= divisor  # Normalization

        """ Select score """
        # print(prob_for_this_run)
        # print(prob_for_this_run.keys())  # iske index karte hain
        keys = []
        values = []
        for k in prob_for_this_run.keys():
            keys.append(k)
            values.append(prob_for_this_run[k])
        ans_lmer = choice(keys, p=values)
        this_score = prob_for_this_run[ans_lmer]
        pos_to_change = separate_seq.find(ans_lmer)
        motif_pos_to_return[sequence_number_to_not_include] = pos_to_change
    # print("loop ran", loop_variable, "times")
    return ans_lmer, seq_with_motif, motif_pos_to_return


def find_consensus(set_seq, len_seq):
    consensus_str = ""

    profile_matrix = dict()
    profile_matrix["A"] = [0] * len_seq
    profile_matrix["T"] = [0] * len_seq
    profile_matrix["G"] = [0] * len_seq
    profile_matrix["C"] = [0] * len_seq

    for seq in set_seq:
        for pos, nucleotide in (enumerate(seq)):
            profile_matrix[nucleotide][pos] += 1
    for i in range(len_seq):
        maxi = 0
        consensus_char = ""
        if profile_matrix["A"][i] > maxi:
            maxi = profile_matrix["A"][i]
            consensus_char = "A"
        if profile_matrix["C"][i] > maxi:
            maxi = profile_matrix["C"][i]
            consensus_char = "C"
        if profile_matrix["T"][i] > maxi:
            maxi = profile_matrix["T"][i]
            consensus_char = "T"
        if profile_matrix["G"][i] > maxi:
            maxi = profile_matrix["G"][i]
            consensus_char = "G"
        consensus_str += consensus_char

    return consensus_str, profile_matrix
# from Q1 creates profile and finds consensus,  returns consensus str and profile


def make_mutation_set(seq, len_seq):
    sequences = set()
    N = ["A", "T", "G", "C"]
    # 0 mutations
    sequences.add(seq)

    # 1 mutation
    for i in range(len_seq):
        for e1 in N:
            temp = seq[:i] + e1 + seq[i + 1:]
            sequences.add(temp)

    s2 = set()
    for s in sequences:
        for i in range(len_seq):
            for e1 in N:
                temp = s[:i] + e1 + s[i + 1:]
                s2.add(temp)

    for s in s2:
        sequences.add(s)

    return sequences


if __name__ == "__main__":
    sequences = []
    n = ["A", "T", "G", "C"]
    no_of_seq = 100
    len_of_seq = 1000
    len_of_motif = 10
    for _ in range(no_of_seq):
        seq = ""
        for __ in range(len_of_seq):
            x = math.floor(random.randint(0, 3))
            # print(x)
            seq += n[x]
        sequences.append(seq)

    motif = ""
    for _ in range(len_of_motif):
        motif += n[math.floor(random.randint(0, 3))]
    print("Planted motif: ", motif)

    motifs = make_mutation_set(motif, len_of_motif)

    """now plant the motif"""
    # for each sequence, find random pos to plant a random motif
    seq_with_motif = []
    for s in sequences:
        motif_pos = random.randint(0, len_of_seq)
        s = s[:motif_pos] + random.sample(motifs, 1)[0] + s[motif_pos:]
        seq_with_motif.append(s)
    # print(len(seq_with_motif), len(seq_with_motif[0]))
    len_seq_with_motif = len(seq_with_motif[0])  # 1010

    """Gibbs sampling"""
    print("Gibbs sampling done!")
    gib, seq_with_motif, motif_pos_to_return = gibbs(seq_with_motif, len_seq_with_motif, len_of_motif)
    print("Consensus motif:", gib)

    print("Sequences: ")
    print(seq_with_motif)
    print("Positions of motif: ")
    print(motif_pos_to_return)

    # print(len(seq_with_motif), len(motif_pos_to_return))

    # print("Motifs:")
    # for i in range(len(seq_with_motif)):
    #     print(seq_with_motif[i][motif_pos_to_return[i]:motif_pos_to_return[i]+10])
