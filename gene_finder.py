# -*- coding: utf-8 -*-
# MINI PROJECT 1: GENE FINDER
#@author: Claire Kincaid

from load import load_seq
dna = load_seq('./data/X73525.fa')

import random
from amino_acids import aa, codons, aa_table   #you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G') #extra doctest, just to cover all bases
    'C'
    >>> get_complement('T') #extra doctest, just to cover all bases
    'A'
    """
    if nucleotide == 'A':
        complement == 'T'
    elif nucleotide == 'C':
        complement == 'G'
    elif nucleotide == 'T':
        complement == 'A'
    elif nucleotide == 'G':
        complement == 'C'
    return complement

    

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("AAAAAAAAA") #control doctest of all one letter
    'TTTTTTTTT'
    >>> get_reverse_complement("ACT") #additional simple doctest
    "AGT"
    """
    i = 0
    Reversedna = ''
    while i<= range(len(dna)):
        Reversedna += get_complement(dna[i])
        i = i + 1
    Reversedna = dna[::-1]
    return Reversedna

    


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    for i in range(len(dna)/3):
        if dna[i*3:i*3+3] == 'TGA':
            ORF = dna[:i*3]
        elif dna[i*3:i*3+3] == 'TAG':
            ORF = dna[:i*3]
        elif dna[i*3:i*3+3] =='TAA':
            ORF = dna[:i*3]
        else:
            ORF = dna
    return ORF
    


def find_all_ORFs_oneframe(dna):

    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC'] 
    >>> find_all_ORFs_oneframe("CATGAATGTAGTTAG") #tests find_all_ORFs w/ no start codon at beginning
    ['ATG']
    """
    allORf = []
    i = 0
    while i < (len(dna)/3):
        if dna[i*3:i*3+3] == 'ATG':
            newORF = rest_of_ORF(dna[i*3:])
            allORf.append(newORF)
        else:
            i = i+1
            i += len(newORF)/3

    return allORf
    


def find_all_ORFs(dna):
    """ 
    Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG'] #I thought this was fine
    """
    allORF = []
    i = 0
    while 0 <= i <= 2:
        ORF = find_all_ORFs_oneframe(dna[i:])
        allORF.append(ORF)
        i = i + 1
    return allORF
    


def find_all_ORFs_both_strands(dna):
    """ 
     Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allORF = []
    allORF_reverse = []
    i = 0
    while 0 <= i <= 2:
        ORF = find_all_ORFs_oneframe(dna[i:])
        allORF.append(ORF)
        reversedna = get_reverse_complement(dna)
        reverseORF = find_all_ORFs_oneframe(reversedna[i:])
        allORF_reverse.append(reverseORF)
        i = i + 1
    both_strands_ORF = allORF + allORF_reverse
    return both_strands_ORF
    


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest_ORF = ''
    for ORF in find_all_ORFs_both_strands(dna):
        if len(ORF) > len(longest_ORF):
            longest_ORF = ORF
    return longest_ORF

    # TODO: implement this
    


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    >>> longest_ORF_noncoding("ATGCGAATGTAGCATCAAA")
    15
    >>> longest_ORF_noncoding("ATGCAATGT")
    6
    """
    longest_ORFs = []
    for i in range(num_trials):
        shuffle_string(dna)
        longest_ORF(dna)
        longest_ORFs.append(longest_ORF)
    longest_ORF_noncoding = max(longest_ORFs, key=len)
    return longest_ORF_noncoding

    


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    #list comprehension, 3 = optional argument, step size for indices
    #list comprehensions are the shit
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)] 
    aa_coding = ''.join([aa_table[codon] for codon in codons if len(codon)==3])
    return aa_coding # I could have made this a list comprehension but I wanted the variable name to make nesting easier
    


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    >>> gene_finder('TAG')
    None
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    find_all_ORFs_both_strands(longest_ORF(dna))
    coding_strand = ''.join([coding_strand_to_AA[orf] for orf in both_strands_ORF if len(orf) >= threshold])
    return coding_strand




    # TODO: implement this
    
#list comprehension, 3 = optional argument, step size for indices
if __name__ == "__main__":
    import doctest
    doctest.testmod()