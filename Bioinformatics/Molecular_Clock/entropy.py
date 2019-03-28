# Entropy is a measure of the uncertainty of a probability distribution
# entropy = - sum(p[i] * log2(p[i]))
import math

profile = {
    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}

def calc_matrix_entropy(profile):
    '''
    Calculates the entropy of the a motif matrix

    INPUT:
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
    
    OUTPUT:
        matrix_entropy(float): the entropy of the motif matrix
    '''
    matrix_entropy = 0
    for i in range(len(profile["A"])):
        entropy = 0
        for nuc in "ACGT":
            p = profile[nuc][i]
            if p != 0.0:
                entropy += (profile[nuc][i]*math.log2(profile[nuc][i]))
        matrix_entropy += -entropy
    return matrix_entropy

print(calc_matrix_entropy(profile))