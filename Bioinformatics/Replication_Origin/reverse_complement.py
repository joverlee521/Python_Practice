def reverse_complement(dna):
    '''
    Find the reverse complement of a DNA string
    
    INPUT: 
        dna(str): a string of bases of a single strand of DNA
    OUTPUT:
        complement(str): a string of bases of reseverse complement
    '''
    dna = dna[::-1]
    base_pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complement = ""
    for nuc in dna:
        complement += base_pairs[nuc]
    return complement

print(reverse_complement("AAAACCCGGT"))