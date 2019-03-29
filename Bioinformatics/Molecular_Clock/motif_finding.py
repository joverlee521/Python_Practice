# Transcription factors regulates gene expression by binding to regulatory motifs 
# in the gene's upstream region, a 600-1000 nucleotide-long region preceding start of gene

# Unlike DnaA box, which clumps within DNA string, a regulatory motif appears once, 
# in several different regions scattered throughout the genome

# We do not know the "ideal" motif, so we find k-mers and score them depending on how similar they are to each other
# Construct large numbers of motif matrices and find collection that minimizes motif score, i.e. the most "conserved" motif matrix
def count_motif_nuc(motifs, pseudo = False):
    '''
    Calculates the count of each nucleotide at each index of the motifs,
    include pseudo = True if want to use pseudocounts in calculation

    INPUT:
        motifs(lst): a list of strings that are potential motifs of the same length
        pseudo(bool): True = include pseudocounts in calculation
    
    OUTPUT:
        count(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of counts of nucleotide at each index of motif
    '''
    count = {}
    # Setting up the key value pairs in count, so each key has a list of 0's with the length of the motifs
    for char in "ATGC":
        init = 0
        # If pseudo = True, initiate with 1 to add pseudocount to the final count
        if pseudo:
            init = 1
        count[char] = [init] * len(motifs[0])
    # Loop through all nucleotides of all motifs and add 1 to count for corresponding nucleotide
    for i in range(len(motifs)):
        for j in range(len(motifs[0])):
            count[motifs[i][j]][j] += 1
    return count

# To get a profile motif, divide all elements in count by length of motifs matrix
# Note that the elements of any column in the profile matrix sum to 1
def profile_matrix(motifs, pseudo = False):
    '''
    Finds the profile motif of the motifs matrix,
    include pseudo = True if want to use pseudocounts in profile

    INPUT:
        motifs(lst): a list of strings that are potential motifs of the same length
        pseduo(bool): True = include pseudocounts in profile
    
    OUTPUT:
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
    '''
    profile = count_motif_nuc(motifs, pseudo)
    for key, value in profile.items():
        l = len(motifs)
        # If including pseudocount, add 4 to divisor to account for the pseudocounts in the total
        if pseudo:
            l += 4
        profile[key] = [num/l for num in value]
    return profile

# Form a consensus motif, which is a motif made of the most popular nucleotide in each column of the motif matrix
def consensus_motif(motifs):
    '''
    Finds the consensus motif of the motif matrix

    INPUT:
        motifs(lst): a list of strings that are potential motifs of the same length
    
    OUTPUT:
        consensus(str): the consensus motif generated from the motif matrix
    '''
    count = count_motif_nuc(motifs)
    consensus = ""
    for i in range(len(motifs[0])):
        max_count = 0
        nuc = ""
        for char in "ATGC":
            if count[char][i] > max_count:
                max_count = count[char][i]
                nuc = char
        consensus += nuc
    return consensus

# Give the motif matrix a score by calculating how different each motif is from the consensus motif
def motifs_matrix_score(motifs):
    '''
    Calculate the score for motif matrix by comparing motifs to consensus motif
    
    INPUT:
        motifs(lst): a list of strings that are potential motifs of the same length

    OUTPUT:
        score(int): score for motif matrix, lower is better
    '''
    score = 0
    consensus = consensus_motif(motifs)
    for i in range(len(motifs)):
        for j in range(len(motifs[0])):
            if motifs[i][j] != consensus[j]:
                score += 1
    return score

def profile_probability(text, profile):
    '''
    Calculate the probability that the text is generated from matrix profile

    INPUT:
        text(str): string to be tested
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
    
    OUTPUT:
        p(float): probability of text generated from profile
    '''
    p = 1
    for i in range(len(text)):
        p *= profile[text[i]][i]
    return p

def profile_most_probable(text, k, profile):
    '''
    Find the kmer within text that has the highest probabiliy to be generated from profile

    INPUT:
        text(str): text in which to find kmer
        k(int): the length of the substring to be found
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
    
    OUTPUT:
        kmer(str): the substring with the highest probability to be generated from profile
    '''
    max_p = -1
    kmer = ""
    for i in range(len(text) - k + 1):
        subtext = text[i:i+k]
        p = profile_probability(subtext, profile)
        if p > max_p:
            max_p = p
            kmer = subtext
    return kmer

def greedy_motif_search(dna, k, t, pseudo = False):
    '''
    Search list of DNA strings to find the best motif matrix of kmers,
    include pseudo = True if want to include pseudocount for search

    INPUT:
        dna(lst): list of DNA strings
        k(int): length of motifs to be found within DNA
        t(int): the length of the dna list
        pseudo(bool): True = include pseudocounts during search
    
    OUTPUT:
        best_motifs(lst): list of kmer strings from each DNA string that has the best motif matrix score
    '''
    best_motifs = []
    # Initializing the "best motif" to just be the first kmer in each DNA string
    for i in range(t):
        best_motifs.append(dna[i][:k])
    n = len(dna[0])
    # This loops through indexes of an individual DNA string
    for i in range(n-k+1):
        motifs = []
        # On each loop, add the next kmer from first DNA string as the first motif
        motifs.append(dna[0][i:i+k])
        # This loops through the remaining DNA strings of the list, starting at 1
        for j in range(1, t):
            # Create a profile based on current motifs in the matrix
            profile = profile_matrix(motifs, pseudo)
            # In the current DNA string, find the most probable kmer based on the matrix profile and add it to the motif matrix
            motifs.append(profile_most_probable(dna[j], k, profile))
        # After constructing the motif matrix from all strings of DNA, determine if it has a better score than the current best matrix
        # Replace best_motifs with new motif matrix if the score is better
        if motifs_matrix_score(motifs) < motifs_matrix_score(best_motifs):
            best_motifs = motifs
    return best_motifs

# ------------------------------------- RANDOMIZED MOTIF SEARCH ------------------------------------- 

def motif_search(dna, profile):
    '''
    Search for best motif matrix for DNA when given a motif profile

    INPUT:
        dna(lst): list of DNA strings
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
    
    OUTPUT:
        motifs(lst): list of motifs from each DNA string that has the best motif matrix score
    '''
    motifs = []
    for i in range(len(dna)):
        motifs.append(profile_most_probable(dna[i], len(profile["A"]), profile))
    return motifs

from random import randint
def random_motifs(dna, k, t):
    '''
    Generate a random motif matrix of kmers from the given DNA list
    
    INPUT:
        dna(lst): list of DNA strings
        k(int): length of motifs to find
        t(int): length of dna list
    
    OUTPUT:
        motifs(lst): a random list of motifs from each DNA string
    '''
    motifs = []
    for i in range(t):
        start = randint(0, len(dna[0])-k)
        motifs.append(dna[i][start:start+k])
    return motifs

def randomized_motif_search(dna, k, t):
    '''
    Find the best motif matrix by starting with a completely random matrix.
    Continue to generate new matrices from the profile until the matrix score stops improving.
    Returns the motif matrix with the best score. 

    INPUT:
        dna(lst): list of DNA strings
        k(int): length of motifs to find
        t(int): length of dna list
    
    OUTPUT:
        best_motifs(lst): the best motif matrix generated from random search
    '''
    motifs = random_motifs(dna, k, t)
    best_motifs = motifs
    while True:
        profile = profile_matrix(motifs, True)
        motifs = motif_search(dna, profile)
        if motifs_matrix_score(motifs) < motifs_matrix_score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs

#  ------------------------------------- GIBBS SAMPLING MOTIF SEARCH  ------------------------------------- 

def normalize(probabilities):
    '''
    Rescale a collection of probabilities such that the probabilities sum to 1

    INPUT & OUTPUT:
        probabilities(dict):
            key = k-mers
            value = floats representing their probabilities
    '''
    sum = 0
    for val in probabilities.values():
        sum += val
    for key, val in probabilities.items():
        probabilities[key] = val / sum
    return probabilities

from random import uniform
def weighted_die(probabilities):
    '''
    Choose a single kmer from the probabilities dict based on each kmer's normalized probability, 
    using a random float generator.

    INPUT:
        probabilities(dict):
            key = k-mers
            value = floats representing their probabilities
    
    OUTPUT:
        kmer(str): the kmer that was choosen based on random float and kmer probabilities
    '''
    kmer = ""
    random_float = uniform(0, 1)
    for key in probabilities.keys():
        random_float -= probabilities[key]
        if random_float <= 0:
            kmer = key
            return kmer

def profile_generated_kmer(text, profile, k):
    '''
    Randomly chooses a kmer from text based on the given profile

    INPUT:
        text(str): DNA string to be evaluated
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
        k(int): the length of the kmer to be returned
    
    OUTPUT:
        kmer(str): random substring with text found based on profile
    '''
    probabilities = {}
    kmer = ""
    for i in range(len(text) - k + 1):
        subtext = text[i:i+k]
        subtext_probability = profile_probability(subtext, profile)
        probabilities[subtext] = subtext_probability
    probabilities = normalize(probabilities)
    kmer = weighted_die(probabilities)
    return kmer

def gibbs_sampler(dna, k, t, n):
    best_motifs = [] 
    motifs = random_motifs(Dna, k, t)
    best_motifs = motifs.copy()
    for _ in range(1, n):
        i = randint(0, t-1)
        profile = profile_matrix(motifs[:i] + motifs[i+1:], True)
        motifs[i] = profile_generated_kmer(dna[i], profile, k)
        if motifs_matrix_score(motifs) < motifs_matrix_score(best_motifs):
            best_motifs = motifs
    return best_motifs


Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
t = 10
k = 15
matrix = greedy_motif_search(Dna, k, t, True)
print("Greedy search results: ")
print(matrix)
print(motifs_matrix_score(matrix))

best_matrix = randomized_motif_search(Dna, k, t)
N = 100
for i in range(N):
    matrix = randomized_motif_search(Dna, k, t)
    if motifs_matrix_score(matrix) < motifs_matrix_score(best_matrix):
        best_matrix = matrix

print("Random search results: ")
print(best_matrix)
print(motifs_matrix_score(best_matrix))

best_matrix = gibbs_sampler(Dna, k, t, N)
for i in range(1, 20):
    matrix = gibbs_sampler(Dna, k, t, N)
    if motifs_matrix_score(matrix) < motifs_matrix_score(best_matrix):
        print("replaced")
        best_matrix = matrix

print("Gibbs Sampler search results: ")
print(best_matrix)
print(motifs_matrix_score(best_matrix))