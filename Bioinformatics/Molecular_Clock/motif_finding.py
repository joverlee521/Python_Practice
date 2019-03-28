# Transcription factors regulates gene expression by binding to regulatory motifs 
# in the gene's upstream region, a 600-1000 nucleotide-long region preceding start of gene

# Unlike DnaA box, which clumps within DNA string, a regulatory motif appears once, 
# in several different regions scattered throughout the genome

# We do not know the "ideal" motif, so we find k-mers and score them depending on how similar they are to each other
# Construct large numbers of motif matrices and find collection that minimizes motif score, i.e. the most "conserved" motif matrix
def count_motif_nuc(motifs):
    '''
    Calculates the count of each nucleotide at each index of the motifs

    INPUT:
        motifs(lst): a list of strings that are potential motifs of the same length
    
    OUTPUT:
        count(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of counts of nucleotide at each index of motif
    '''
    count = {}
    # Setting up the key value pairs in count, so each key has a list of 0's with the length of the motifs
    for char in "ATGC":
        count[char] = [0] * len(motifs[0])
    # Loop through all nucleotides of all motifs and add 1 to count for corresponding nucleotide
    for i in range(len(motifs)):
        for j in range(len(motifs[0])):
            count[motifs[i][j]][j] += 1
    return count

# To get a profile motif, divide all elements in count by length of motifs matrix
# Note that the elements of any column in the profile matrix sum to 1
def profile_matrix(motifs):
    '''
    Finds the profile motif of the motifs matrix

    INPUT:
        motifs(lst): a list of strings that are potential motifs of the same length
    
    OUTPUT:
        profile(dict):
            key = nucleotides("A", "T", "C", "G")
            value = list of ratios of nucleotide at each index of motif
    '''
    l = len(motifs)
    profile = count_motif_nuc(motifs)
    for key, value in profile.items():
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

def greedy_motif_search(dna, k, t):
    '''
    Search list of DNA strings to find the best motif matrix of kmers

    INPUT:
        dna(lst): list of DNA strings
        k(int): length of motifs to be found within DNA
        t(int): the length of the dna list
    
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
            profile = profile_matrix(motifs)
            # In the current DNA string, find the most probable kmer based on the matrix profile and add it to the motif matrix
            motifs.append(profile_most_probable(dna[j], k, profile))
        # After constructing the motif matrix from all strings of DNA, determine if it has a better score than the current best matrix
        # Replace best_motifs with new motif matrix if the score is better
        if motifs_matrix_score(motifs) < motifs_matrix_score(best_motifs):
            best_motifs = motifs
    return best_motifs