import urllib.request
from collections import defaultdict

def pattern_count(text, pattern):
    '''
    Counts the number of times a pattern appears in a text.
    Includes overlapping occurrences of the pattern

    INPUT:
        text(str): the string to check for pattern
        pattern(str): the pattern/substring to look for in text

    OUTPUT:
        count(int): the number of times the pattern appears in text
    '''
    count = 0
    # Loop through all indicies of text, stopping at the last index that still fits the length of the pattern
    for i in range(len(text) - len(pattern) + 1):
        # Look at substring in text at index to see it is equal to the pattern
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

def frequency_map(text, k):
    '''
    Creates a frequency map of all possible k length patterns within a text

    INPUT:
        text(str): the string to map through
        k(int): the length of the substrings to find within text
    
    OUTPUT:
        freq(dict): 
            key = k-length substrings of text
            value = count of substring within text
    '''
    freq = defaultdict(lambda: 0)
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        freq[pattern] += 1
    return freq

def most_frequent_pattern(text, k):
    '''
    Find the most frequent pattern of length k within text

    INPUT:
        text(str): the string to look for most frequent pattern
        k(int): the length of the pattern to find within text
    OUTPUT:
        patterns(lst): a list of k-length patterns with the highest frequency
    '''
    patterns = []
    freq = frequency_map(text, k)
    max_freq = max(freq.values())
    for pattern in freq:
        if freq[pattern] == max_freq:
            patterns.append(pattern)
    return patterns

def pattern_matching(text, pattern):
    '''
    Find all occurrences of a pattern in text

    INPUT:
        text(str): the string to check for pattern
        pattern(str): the pattern/substring to look for in text
    OUTPUT:
        positions(lst): list of starting indices of pattern within text
    '''
    positions = []
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions

def symbol_to_number(symbol):
    '''
    Convert nucleotides("A", "C", "G","T") to numbers(0, 1, 2, 3)
    '''
    for i, char in enumerate("ACGT"):
        if symbol == char:
            return i


def pattern_to_number(pattern):
    '''
    Convert pattern to a number that represents the patterns index in a list of lexicographically ordered patterns

    INPUT:
        pattern(str): pattern to convert to number
    
    OUTPUT:
        (int): integer representing the pattern
    '''
    if len(pattern) == 0:
        return 0
    symbol = pattern[len(pattern)-1]
    prefix = pattern[:len(pattern)-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)

def number_to_symbol(number):
    '''
    Convert number(0,1,2,3) to nucleotide symbols("A", "C", "G", "T")
    '''
    for i, char in enumerate("ACGT"):
        if number == i:
            return char

def number_to_pattern(number, k):
    '''
    Convert a number to a nucleotide pattern
    
    INPUT:
        number(int)
        k(int): the length of the pattern
    
    OUTPUT:
        (str): pattern represented by number
    '''
    if k == 1:
        return number_to_symbol(number)
    prefix_number = number // 4
    remainder = number % 4
    prefix_pattern = number_to_pattern(prefix_number, k-1)
    return prefix_pattern + number_to_symbol(remainder)

def computing_frequencies(text, k):
    '''
    Creates a frequency list of all possible k length patterns within a text

    INPUT:
        text(str): the string to map through
        k(int): the length of the substrings to find within text
    
    OUTPUT:
        freq_array(lst): a list of all the counts of k-length patterns withing text,
        with each index pointing to the number representing the pattern
    '''
    freq_array = [0] * (4**k)
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        j = pattern_to_number(pattern)
        freq_array[j] += 1
    return freq_array

def pattern_clump_finder(genome, k, L, t):
    '''
    Find kmers that form (L,t)-clumps in the genome. 
    Meaning, we are finding kmer that appear at least t times within substring of length L in the genome

    INPUT:
        genome(str): the genome to be analyzed
        k(int): the length of patterns to find within genome
        L(int): the length of the substring that encompasses a clump
        t(int): the number of times a kmer needs to appear within substring
    
    OUTPUT:
        kmers(set): a set of kmers that form clumps within genome
    '''
    kmers = set()
    substring = genome[:L]
    freq_map = frequency_map(substring, k)
    for key, value in freq_map.items():
        if value >= t:
            kmers.add(key)
    for i in range(1, len(genome) - L + 1):
        first_pattern = genome[i-1:i-1+k]
        freq_map[first_pattern] -= 1
        last_pattern = genome[i+L-k:i+L]
        freq_map[last_pattern] += 1
        if freq_map[last_pattern] >= t:
            kmers.add(last_pattern)
    return kmers

if __name__ == "__main__":
    text = urllib.request.urlopen("http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt").read()
    k = 9
    L = 500
    t = 3
    # pattern_clump = pattern_clump_finder(text, k, L, t)
    # print(len(pattern_clump))