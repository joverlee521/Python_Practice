# Hamming Distance = the number of positions at which the corresponding chars are different between two strings of equal length

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

