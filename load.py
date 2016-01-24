# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 22:02:04 2014

@author: pruvolo
"""

from os import path


def load_seq(fasta_file):
    """ Reads a FASTA file and returns the DNA sequence as a string.

    fasta_file: the path to the FASTA file containing the DNA sequence
    returns: the DNA sequence as a string
    """
    retval = ""
    f = open(fasta_file)
    lines = f.readlines()
    for l in lines[1:]:
        retval += l[0:-1]
    f.close()
    return retval