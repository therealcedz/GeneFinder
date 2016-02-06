# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Cedric Kim

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')		#testing if T/ G works
    'A'
    >>> get_complement('G')
    'C'
    >>> get_complement('H')		#testing weird cases
    'none'

    """
    if nucleotide == 'A':
    	return 'T'
    elif nucleotide == 'C':
    	return 'G'
    elif nucleotide == 'T':
    	return 'A'
    elif nucleotide == 'G':
    	return 'C'
    else:
    	return 'none'


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
    reverse_string = ''
    for i in range(len(dna)):										# for each index,
    	reverse_string += get_complement(dna[len(dna) - 1 - i])		# take the complement of the len-index
    return reverse_string											# and put into this new reverse_string

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
    stop_codon = ['TGA', 'TAG', 'TAA']				#initiate stop codons
    found = False 									# run if substring is not found 
    i = 0											# initiate index
    while not found and len(dna) > (i-1)*3:			#while substring is not found, or index goes over length
        sub_string = dna[i*3: (3*(i+1))]			#place every 3 letters into substring
    	if (sub_string == stop_codon[0] or sub_string == stop_codon[1] or sub_string == stop_codon[2]):	#check if substring is found
   			found = True
    	i += 1
    return dna[0: 3*(i-1)]

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

    >>> find_all_ORFs_oneframe("ATGTAGAAAATGAAA")			#checks if it works for small sections
    ['ATG', 'ATGAAA']
    >>> find_all_ORFs_oneframe("ATGCATTAGAATGAAATAG")			#checks if it works for offset ATG
    ['ATGCAT']
    """
    dna_list = [];
    temp_dna = dna
    finished = False
    while not finished:											#while not finished
    	if(len(temp_dna) == 0):									#if the length is equal to 0, exit loop
    		finished = True
    	elif(temp_dna[:3] == 'ATG'):
            rest_of_ORF_list = rest_of_ORF(temp_dna)			#if we find starting ATG at multiples of 3,
            dna_list.append(rest_of_ORF_list)			          #add the ORF string to the list
            temp_dna = temp_dna[len(rest_of_ORF_list)+3:]         #cut off the ORF string
    	else:													
    		temp_dna = temp_dna[3:]								#else, cut off first 3 letters
    return dna_list
def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequenc
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    dna_list = []
    for i in range(0,3):										#run 3 times
    	dna_list.extend(find_all_ORFs_oneframe(dna[i:]))			#find the first frame
    return dna_list

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    dna_list = []
    dna_list.extend(find_all_ORFs(dna))							#add the list
    dna_list.extend(find_all_ORFs(get_reverse_complement(dna)))	#add the complement
    return dna_list

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest_string = ''
    ORFs_both_strands = find_all_ORFs_both_strands(dna)
    for string in ORFs_both_strands:
        if len(string) > len(longest_string):
            longest_string = string
    return longest_string
#print(longest_ORF("ATGCGAATGTAGCATCAAA"))
#print(longest_ORF("AATATTGAGCAAGACATGC"))
def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: they maximum length longest ORF """
    num_char = 0
    length = 0
    for i in range(num_trials):
        dna = shuffle_string(dna)
        length = len(longest_ORF(dna))
        if num_char < length:
            num_char = length
    return num_char

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
    amino_acid_string = ''
    temp_dna = dna
    finished = False
    while not finished:                                        
        if(len(temp_dna) < 3):                                 
            finished = True
        else:
            amino_acid_string += aa_table[temp_dna[:3]]
            temp_dna = temp_dna[3:]
    return amino_acid_string

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    amino_acid_list = []
    threshold = longest_ORF_noncoding(dna, 1500)
    for string in find_all_ORFs_both_strands(dna):
        if len(string) > threshold:
            amino_acid_list.append(coding_strand_to_AA(string))
    return amino_acid_list

print(gene_finder(dna))
aa = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
      '|', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R',
      'G']

codons = [['TTT', 'TTC'],
          ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
          ['ATT', 'ATC', 'ATA'],
          ['ATG'],
          ['GTT', 'GTC', 'GTA', 'GTG'],
          ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
          ['CCT', 'CCC', 'CCA', 'CCG'],
          ['ACT', 'ACC', 'ACA', 'ACG'],
          ['GCT', 'GCC', 'GCA', 'GCG'],
          ['TAT', 'TAC'],
          ['TAA', 'TAG', 'TGA'],
          ['CAT', 'CAC'],
          ['CAA', 'CAG'],
          ['AAT', 'AAC'],
          ['AAA', 'AAG'],
          ['GAT', 'GAC'],
          ['GAA', 'GAG'],
          ['TGT', 'TGC'],
          ['TGG'],
          ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
          ['GGT', 'GGC', 'GGA', 'GGG']]

# create a dictionary lookup table for mapping codons into amino acids
aa_table = {}
for i in range(len(aa)):
    for codon in codons[i]:
        aa_table[codon] = aa[i]


if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(find_all_ORFs, globals())