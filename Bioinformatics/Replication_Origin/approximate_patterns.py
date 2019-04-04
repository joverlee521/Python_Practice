# Hamming Distance = the number of positions at which the corresponding chars are different between two strings of equal length
from collections import defaultdict
from reverse_complement import reverse_complement

def hamming_distance(string1, string2):
    '''
    Calculates the Hamming distance between two strings

    INPUT:
        Two strings of equal length

    OUTPUT:
        ham_distance(int)
    '''
    ham_distance = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            ham_distance += 1
    return ham_distance

def approx_pattern_matching(text, pattern, d):
    '''
    Finds all approximate occurrences of pattern in text with at most d mismatches

    INPUT:
        text(str): string to analyze
        pattern(str): pattern to look for in text
        d(int): max Hamming Distance between approx patterns and given pattern
    
    OUTPUT:
        positions(lst): list of starting index of approximate patterns
    '''
    positions = []
    for i in range(len(text) - len(pattern) + 1):
        sub_text = text[i:i+len(pattern)]
        if sub_text == pattern or hamming_distance(sub_text, pattern) <= d:
            positions.append(i)
    return positions 

def approx_pattern_count(text, pattern, d):
    '''
    Calculates the number of approximate occurrences of pattern in text with at most d mismatches

    INPUT:
        text(str): string to analyze
        pattern(str): pattern to look for in text
        d(int): max Hamming Distance between approx patterns and given pattern
    
    OUTPUT:
        count(int): the number of times the approx patttern appears in text with at most d mismatches
    '''
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        sub_text = text[i:i+len(pattern)]
        if sub_text == pattern or hamming_distance(sub_text, pattern) <= d:
            count += 1
    return count

def pattern_neightbors(pattern, d):
    '''
    Find all patterns that are at most d Hamming distanace from pattern

    INPUT:
        pattern(str): the original pattern string
        d(int): the max difference from the given pattern
    
    OUTPUT:
        neighbors(set): a set of all patterns d away from given pattern
    '''
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return set("ACGT")
    neighbors = set()
    suffix_neighbors = pattern_neightbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hamming_distance(pattern[1:], text) < d:
            for nuc in "ACGT":
                neighbors.add(nuc + text)
        else:
            neighbors.add(pattern[0] + text)
    return neighbors

def approx_frequency_map(text, k , d):
    '''
    Creates a frequency map of all possible k length patterns within a text,
    including counts of approximate patterns that are at most d Hamming distance away from pattern,
    and accounts for the reverse complement of the approximate patterns 

    INPUT:
        text(str): the string to map through
        k(int): the length of the substrings to find within text
        d(int): the max difference betweetn pattern and other strings to count as pattern
    
    OUTPUT:
        freq(dict): 
            key = k-length substrings of text
            value = count of substring within text
    '''
    freq = defaultdict(lambda: 0)
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        neighbors = pattern_neightbors(pattern, d)
        reverse_neighbors = pattern_neightbors(reverse_complement(pattern), d)
        for neighbor in neighbors:
            freq[neighbor] += 1 
        for reverse in reverse_neighbors:
            freq[reverse] += 1
    return freq
    

def most_frequent_approx_pattern(text, k, d):
    '''
    Finds the most frequent kmers with up to d mismatches in the text
    
    INPUT:
        text(str): string to analyze
        k(int): the length of patterns to find within text
        d(int): the number of mismatches two strings can have to still count as the same pattern
    
    OUTPUT:
        kmers(lst): list of most frequent kmers with up to d mismatches
    '''
    kmers = []
    freq = approx_frequency_map(text, k, d)
    max_freq = max(freq.values())
    for pattern in freq:
        if freq[pattern] == max_freq:
            kmers.append(pattern)
    return kmers

if __name__ == "__main__":
    text = "TCCTTTCCTCCCACTTTCCTCCCACTTTCCTCCACTCCTCCGCGTTTCCTCCTTTCTCCTCTCCTCCCACCACGCGTCTCCTTTCGCGGCGTCCGCGTCCTTGCGGCGCACTCCTTGCGTTTCCCACCACTTGCGTCCTCCTCCCACTCCTCCTCCTCCGCGCACTTCACGCGGCGTTTCCTTGCGTCCTCCACGCGGCGTTTCCGCGTCCACTCCCACTCCGCGCACTCTCCAC"
    k = 5
    d = 2
    print(" ".join(x for x in most_frequent_approx_pattern(text, k, d)))