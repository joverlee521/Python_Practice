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
    freq = {}
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        if pattern in freq:
            freq[pattern] += 1
        else:
            freq[pattern] = 1
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

