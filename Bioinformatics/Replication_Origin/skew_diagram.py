# DNA replication is asymmetrical due to DNA polymerase only able to move in (3'-> 5') direction
# Reverse half strands only need one primer, and proceeds quickly, double-stranded most of the time
# Forward half strands need to wait for replication forks to open, so another primer can be put down, and form Okazaki fragments
# During the wait for replication forks to open, the forward half strand is single-stranded!! 
    # This has a higher tendency to mutate by deamination of cytosine(C) to thymine(T)
    # The mismatched T-G pairs then get repaired to T-A pairs
    # Overall, this leads to a decrease of C in the forward strand and a decrease in G in the reverse strand

# Analyzing Genome Halfstrands: since bacterial DNA is circular, we have to account for windows that wrap around the end of genome
from pattern_count_frequency import pattern_count

def symbol_array(genome, symbol):
    '''
    INEFFICIENT VERSION!
    Counts how many times a symbol appears within a half the genome starting at each index

    INPUT:
        genome(str): string of bases of genome
        symbol(str): the specific base to count
    
    OUTPUT:
        array(dict): 
            key = index within genome
            value = count of symbol within half of genome starting at key index
    '''
    array = {}
    n = len(genome)
    extended_genome = genome + genome[:n//2]
    for i in range(n):
        array[i] = pattern_count(extended_genome[i:i+(n//2)], symbol)
    return array

def faster_symbol_array(genome, symbol):
    '''
    EFFICIENT VERSION!
    Counts how many times a symbol appears within a half the genome starting at each index

    INPUT:
        genome(str): string of bases of genome
        symbol(str): the specific base to count
    
    OUTPUT:
        array(dict): 
            key = index within genome
            value = count of symbol within half of genome starting at key index
    '''
    array = {}
    n = len(genome)
    extended_genome = genome + genome[:n//2]
    array[0] = pattern_count(genome[:n//2], symbol)
    for i in range(1, n):
        # set new value equal to previous value
        array[i] = array[i-1]
        # new value can differ from previous value by losing a symbol and/or gaining symbol
        if extended_genome[i-1] == symbol:
            # if previous char == symbol, then count decreases by 1
            array[i] -= 1
        if extended_genome[i+(n//2)-1] == symbol:
            # if char added at the end of the half genome is equal to symbol, then count increase by 1
            array[i] += 1
    return array

def skew_array(genome):
    '''
    Denotes the occurrences of C and G within genome up to each index.
    skew[i] is equal to the occurrence of G subtracted by occurrence of C
    
    INPUT:
        genome(str): a string of bases to analyze
    
    OUTPUT:
        skew(lst): list of array Skew for genome, starts with 0 at skew[0]
    '''
    skew = [0]
    for i in range(len(genome)):
        change = 0
        if genome[i] == "C":
            change = -1
        elif genome[i] == "G":
            change = 1
        skew.append(skew[i] + change)
    return skew

# The skew decreases along the reverse half strand (ter -> ori) but increases along forward half strand (ori -> ter):
    # SO: we can assume that the ori is located near the minimum of the skew!!
def minimum_skew(genome):
    '''
    Find all positions within genome where skew reaches minimum
    
    INPUT:
        genome(str): genome to analyze
    
    OUTPUT:
        positions(lst): list of all indices where skew is minimum
    '''
    positions = []
    array = skew_array(genome)
    min_skew = min(array)
    for i in range(len(array)):
        if array[i] == min_skew:
            positions.append(i)
    return positions
