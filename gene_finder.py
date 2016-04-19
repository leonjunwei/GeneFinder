"""
@author: Leon Lam (@leonjunwei)
Last updated: 01/30/16

This code analyzes a DNA sequence and outputs snippets of DNA that are likely to be protein-coding genes.
"""


import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
from load import load_contigs
contigs = load_contigs()

mycontig = contigs[2][1]
newcontig = "".join(['A' if x not in 'ATGC' else x for x in mycontig])
# print newcontig

def shuffle_string(s):
    """Shuffles the characters in the input string
        # NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

# compl = {"A":"T","T":"A","C":"G","G":"C"}
# def get_complement(nucleotide):
#     """ Returns the complementary nucleotide

#         nucleotide: a nucleotide (A, C, G, or T) represented as a string
#         returns: the complementary nucleotide
#     >>> get_complement('A')
#     'T'
#     >>> get_complement('C')
#     'G'
#     """

#     try:
#         return compl[nucleotide]
#     except KeyError:
#         return "X"

    # if nucleotide == 'A':
    #     return 'T'
    # elif nucleotide == 'T':
    #     return 'A'
    # elif nucleotide == 'C':
    #     return 'G'
    # elif nucleotide == 'G':
    #     return 'C'

#print get_complement('C')

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    result = ""
    compl = {"A":"T","T":"A","C":"G","G":"C"}
    for i in dna:
        result+=compl[i]
    return result[::-1]


# print get_reverse_complement("CCGCGTTCA")

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
    stopcodon = ['TAG','TAA','TGA']
    for i in range(0,len(dna)-2,3):
        if dna[i:i+3] in stopcodon: #dna is formatted as "ATGAGATAGG", so list slicing on dna gives a concatenated string
            return dna[0:i]
    return dna


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
    """
    result = []
    i = 0
    stopcodon = ['TAG','TAA','TGA']
    while i<len(dna)-2:
        if dna[i:i+3] == "ATG": #Be very careful about the variable type! ["ATG"] is not the same as "ATG".
            result.append(rest_of_ORF(dna[i:]))
            i = next((m for m in range(i+3,len(dna)-2, 3) if dna[m:m+3] in stopcodon), len(dna)-2)
            # i = y
        else:
            i+=3
    return result

# print find_all_ORFs_oneframe("CATATGCATATGTGTAGATAGATGTGCCC")


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    result = []
    for i in range(0,3):
        result+=(find_all_ORFs_oneframe(dna[i:])) #find_all_ORFs outputs a list, and we want to add the elements and not append a list to another list.
    return result

# print find_all_ORFs("ATGCATGAATGTAG")

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    result = (find_all_ORFs(dna)) + (find_all_ORFs(get_reverse_complement(dna)))
    result = [x for x in result if x != []]
    return result

# print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    store = {}
    largestlength = 0
    for i in find_all_ORFs_both_strands(dna):
        if len(i)>largestlength:
            store[len(i)] = i
            largestlength = len(i)
    return store[largestlength]

# print longest_ORF("ATGCGAATGTAGCATCAAA")

# dna = load_seq("/home/leon/GeneFinder/data/X73525.fa")
# print dna

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    count = 0
    for i in range(num_trials):
        x = len(longest_ORF(shuffle_string(dna)))
        if x > count:
            count = x
    return count

# print longest_ORF_noncoding(newcontig, 1)


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
    aaList = [aa_table[dna[i:i+3]] for i in range(0,len(dna)-2,3)]
    return "".join(aaList)
    
# print coding_strand_to_AA("ATGCCCGCTTT")


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = 1500#longest_ORF_noncoding(dna, 1)
    orfList = find_all_ORFs_both_strands(dna)
    aaList = [coding_strand_to_AA(c) for c in orfList if len(c)>threshold]
    return aaList

print gene_finder(newcontig)


# if __name__ == "__main__":
#     import doctest
#     doctest.testmod()
