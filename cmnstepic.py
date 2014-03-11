import numpy as np




def revcomp(t):
    '''Generates the reverse complement of a DNA sequence'''
    
    rcdict = {
    'A':'T', 
    'T':'A', 
    'C':'G', 
    'G':'C'}
    
    base = 4
    comp = ''
    for b in t:
        c = rcdict[b]
        comp = comp+c
    return(comp[::-1])
