import re
import numpy as np
import itertools
import math
import cProfile
import pstats

#t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
#d = 1
#k = 4

#t = 'CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC'
#k = 10
#d = 2

#t = 'GCTGCGACAAGCGAGCGAGCTCAAGCGAGTGGCTGTGGCTGTGGCGATACCAAGCGATACCAATACTACGCTCAAGTGTACGCGAGCGAGCGATACGCTGCGACAAGTGGCTCAAGCGACAAGCGACAAGTGGCGAGTGGCGAGCTGCTGTGGTGGCGATACGTGGCTGCTGTGGCGACAACAATACGCTCAAGTGGTGGCTGCTGCTTACTACCAACAATACCAACAACAACAAGCGACAAGTGCAATACGTGGCTGCGAGCGAGTGTACGCTGCGACAAGCTGCTCAAGTGGCTTACTACTACGTGCAATACGCTTACGCTGTGGCTTACTACCAAGTG'
#k = 8
#d = 2

#t = 'CCGAAGCCGACACACACCGACAAAGAAGCACACCGGCCCACACACACCGCACAAAGAAGACAGCCGCCAAGAAGCCGGCCCCGGCCACACACACACAAAGCCGGCCGCCGCCACACACACCGAAGGCCCACAGCCCCGAAGGCCGCCGCCGCCAAGACACACACCGAAGAAGGCCCCGACAAAGCCGGCCACACACACCGAAGCCGGCCCCGCACACACACCGCACAGCCGCCAAGCCGCACACACAAAGAAGCACAACACCGACACCGCCGACACCGACAACACCGACAGCCCCGACACCGCCGACACACAACACACACACAGCCCACAACAAAGAAGACACCGCACAAAGCACAAAGGCCAAGAAGCCGCCG'
#k = 10
#d = 2
#ans = 'CACACACACA'

#t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
#k = 4
#d = 1

t = 'CTTGCCGGCGCCGATTATACGATCGCGGCCGCTTGCCTTCTTTATAATGCATCGGCGCCGCGATCTTGCTATATACGTACGCTTCGCTTGCATCTTGCGCGCATTACGTACTTATCGATTACTTATCTTCGATGCCGGCCGGCATATGCCGCTTTAGCATCGATCGATCGTACTTTACGCGTATAGCCGCTTCGCTTGCCGTACGCGATGCTAGCATATGCTAGCGCTAATTACTTAT'
k = 9
d = 3

#ans = ['AGCGCCGCT', 'AGCGGCGCT']

#t = 'CTATTCGACATCTATTCTTCCGATCTAGACCGCTAATGACCTAATGACGACATATCTACGCTAATATTTCATGACGACCGCGTTCGACGACTTCTTCATTTCTTCGACTTCTTCATGACATGACATATCGATCGTTCCGCTAGACATATCTATTCATGACCGCGTTCCGATGACATCGCGTTCGACCGATCTAATCTAGACTTCATATGACGAC'
#k = 9
#d = 2

b = {
'0': 'A', 
'1': 'T', 
'2': 'C', 
'3': 'G'
}

rcdict = {
'A':'T', 
'T':'A', 
'C':'G', 
'G':'C'}
base = 4

def approxpm(p, t, d):
    '''For Text t and Pattern p, finds sequences within t that approximately 
    match p with at most d mismatches. Returns a list of positions which are given 
    as the starting index of the matching sequence within t.

    Algorithm works by sliding a window of length(p) along t and comparing each 
    window to p. If there are greater than d mismatches, it stops counting and 
    proceeds to next window.
    '''
    indices = []

    for n in range(len(t)):
        window = t[n:n+len(p)]
        if len(window) < len(p):
            break
        counter = 0
        for o in range(len(p)):
            if p[o] != window[o]:
                counter = counter+1
            if counter > d:
                break
        if counter <= d:
            indices.append(n)
        
        #count = hamdist(window, p)
        #if hamdist <= d:
            #indices.append(n)
        #print('counter', counter)
        
    return indices

#def genkmers (k):
    #'''Returns a list of all possible DNA sequences of length k.

    #Algorithm works by listing all the possible sequences as numbers in base 4. 
    #The number of kmers is 4**k, and each number up to this one is converted to 
    #base 4. Then the numbers are exchanged with bases as listed in the dictionary 
    #'b'.
    #'''

    #kmers = []
    #numkmers = 4**k

    #for i in range(numkmers):
        #numform = num2base(i, 4)
        #lendiff = len(num2base(numkmers, 4)) - len(numform)-1
        #if lendiff > 0:
            #numform = '0'*lendiff + numform
            
        #for k, v in b.iteritems():
            #numform = numform.replace(k, v)
        ##numform = numform.replace('0', 'A')
        #kmers.append(numform)

    #return kmers


def writetofile(anslist, ansfile):
    with open(ansfile, 'w') as f:
        for ans in anslist:
            f.write('{0} '.format(ans))
    
#freqkmers = freqapproxpm(t, d, k, base)
#writetofile(freqkmers, 'approxpmfreqans.txt')

bases = ['C', 'G', 'A', 'T']

#ex = 'CAC'

def changebases(listk, iset, iseti, newkmers):
    '''Given a list of bases (listk), finds all possible sequences where the 
    bases at position iset is substituted by every base. iseti is a counter that 
    should begin at 0. newkmers is a dictionary which the new sequences will be 
    added to.
    '''
    
    #listk = list(kmer)
    #print 'firstkmer', listk
    #print listk[iset[iseti]]
    #print(iset[iseti])
    #inds.append(iset)
    #print iset
    for b in bases:
        #print iseti
        #print 'newbase', b
        #print 'base to be changed', listk[iset[iseti]]
        listk[iset[iseti]] = b
        newkmer = ''.join(listk)
        newkmers.append(newkmer)
        #print newkmers
        #print listk
        if iseti != len(iset)-1:
            changebases(listk, iset, iseti+1, newkmers)
        if iseti == len(iset)-1:
            continue
            
        #count = count-1
        #print count
        #if count > 0:
            #iseti = iseti+1
            #changebases(kmer, iset, iseti, count)
        #newkmer = ''.join(listk)
        #newkmers.append(newkmer)
    #print 'end newkmers', newkmers
    #print 'len newkmers', len(newkmers)
    return newkmers

def genkmers(kmer, poskmers, d):
    '''For a sequence kmer, finds all sequences that differ from kmer by at 
    most d mismatches. Adds all results to and then returns the dictionary 
    poskmers.
    
    Algortihm works by using itertools.combinations() to find every 
    combination of indices for the positions that must be substituted in the 
    kmer. Then it substitues every combination of bases at those positions 
    and adds them as keys to the dictionary poskmers.
    '''
    
    #inds = []
    isets = itertools.combinations(range(len(kmer)), d)
    for iset in isets:
        #print iset
        listk = list(kmer)
        tn = []
        newkmers = changebases(listk, iset, 0, tn)
        #print 'final newkmers', newkmers
        for newkmer in newkmers:
            if newkmer not in poskmers:
                poskmers[newkmer] = []
            
    
    
    #for i in range(len(kmer)):
        #for b in bases:
            ##print b
            #listk = list(kmer)
            ##if listk[i] != b:
            #listk[i] = b
            #newkmer = ''.join(listk)
            #if newkmer not in poskmers:
                #poskmers[newkmer] = 0
            
            #if counter != 0:
                #genkmers(newkmer, poskmers, d, counter-1)
    #print(inds)
    return poskmers

def revcomp(rcdict, t):
    '''Generates the reverse complement of a DNA sequence'''
    
    comp = ''
    for b in t:
        c = rcdict[b]
        comp = comp+c
    return(comp[::-1])
    
def findfreq(t, k, d):

    allposkmers = {}
    c = {}
    for n in range(len(t)):
        window = t[n:n+k]
        #print(window)
        if len(window) < k:
            break
        poskmerstemp = {}
        poskmers = genkmers(window, poskmerstemp, d)
        #if n == 0:
            ##print(poskmers)
        
        #for pk in poskmers.keys():
            #if pk not in c:
                #c[pk] = 1
        
        for pk in poskmers.keys():
            if pk not in allposkmers:
                allposkmers[pk] = 1

    for apk in allposkmers.keys():
        apmatches = len(approxpm(apk, t, d))
        c[apk] = apmatches
            #apmatches = len(approxpm(pk, t, d)) 
            ##print(inds)
            #c[pk] = c[pk]+apmatches
    print len(c)
    #print c['']
    
    e = {}
    for key, val in c.iteritems():
        if val not in e:
            e[val] = []
        e[val].append(key)
    print max(e.keys())
    print e[max(e.keys())]


def findfreq2(t, k, d):
    '''For Text t, finds the most frequent kmers of length k that appear in 
    the text with at most d mismatches.
    
    Algorithm works by sliding a window along the text of length k, and then 
    finding all the possible sequences that differ from the window sequence 
    by at most d mismatches. These sequences are added as keys to a dictionary and as 
    the window slides along, each occurence of the sequence is added to the 
    value of the respective key. Then a new dictionary is generated with the 
    keys as values of the first dictionary and the values as keys, so that 
    the keys are the # of times each sequence is represented in the first 
    dictionary. Then it finds the key with the highest value in this new 
    dictionary and returns the sequences associated with that key.
    '''

    allposkmers = {}
    c = {}
    for n in range(len(t)):
        window = t[n:n+k]
        #print(window)
        if len(window) < k:
            break
        poskmerstemp = {}
        poskmers = genkmers(window, poskmerstemp, d)
        #if n == 0:
            ##print(poskmers)
        
        #for pk in poskmers.keys():
            #if pk not in c:
                #c[pk] = 1
        
        for pk in poskmers.keys():
            if pk not in allposkmers:
                allposkmers[pk] = 1
            allposkmers[pk] = allposkmers[pk]+1

    #for apk in allposkmers.keys():
        #apmatches = len(approxpm(apk, t, d))
        #c[apk] = apmatches
            #apmatches = len(approxpm(pk, t, d)) 
            ##print(inds)
            #c[pk] = c[pk]+apmatches
    print len(c)
    #print c['']
    
    e = {}
    for key, val in allposkmers.iteritems():
        if val not in e:
            e[val] = []
        e[val].append(key)
    print max(e.keys())
    print e[max(e.keys())]


def findfreqrc(t, k, d, rcdict):
    ''''
    Identail to findfreq2 except that the analysis is also run on the reverse 
    complement.
    '''
    
    allposkmers = {}
    trc = revcomp(rcdict, t)
    #print 't', t
    #print 'trc', trc
    #c = {}
    for n in range(len(t)):
        window = t[n:n+k]
        #print(window)
        if len(window) < k:
            break
        poskmerstemp = {}
        poskmers = genkmers(window, poskmerstemp, d)
        
        for pk in poskmers.keys():
            if pk not in allposkmers:
                allposkmers[pk] = 1
            allposkmers[pk] = allposkmers[pk]+1
    
    for n in range(len(trc)):
        window = trc[n:n+k]
        #print(window)
        if len(window) < k:
            break
        poskmerstemp = {}
        poskmers = genkmers(window, poskmerstemp, d)
        
        for pk in poskmers.keys():
            if pk not in allposkmers:
                allposkmers[pk] = 1
            allposkmers[pk] = allposkmers[pk]+1
    
            
    #print allposkmers
    #for apk in allposkmers.keys():
        #apmatches = len(approxpm(apk, t, d)) + len(approxpm(apk, trc, d)) 

    ##print(inds)
        #c[apk] = apmatches
    #print c
    #for n in range(len(trc)):
        #window = t[n:n+k]
        ##print(window)
        #if len(window) < k:
            #break
        #poskmerstemp = {}
        #poskmers = genkmers(window, poskmerstemp, d)
        ##if n == 0:
            ###print(poskmers)
        
        #for pk in poskmers.keys():
            #if pk not in c:
                #c[pk] = 1

            #apmatches = len(approxpm(pk, t, d))
            ##print(inds)
            #c[pk] = c[pk]+apmatches

    #print c['']
    #print c['ACAT']
    e = {}
    for key, val in allposkmers.iteritems():
        if val not in e:
            e[val] = []
        e[val].append(key)
    print max(e.keys())
    print e[max(e.keys())]

def nchoosek(n, k):
    '''Implements the function n choose k.
    '''
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

if __name__ == "__main__":
    #myans = findfreqrc(t, k, d, rcdict)
    #myans == ans
    ##findfreq(t, k, d)
    #cProfile.run('findfreqrc(t, k, d, rcdict)', 'prof_approxpmfreqwords3')
    #p = pstats.Stats('prof_approxpmfreqwords3')
    #p.strip_dirs().sort_stats(-1).print_stats()
    #p.sort_stats('cumulative').print_stats(10)

    cProfile.run('findfreqrc(t, k, d, rcdict)', 'prof_approxpmfreqwords3')
    p = pstats.Stats('prof_approxpmfreqwords3')
    p.strip_dirs().sort_stats(-1).print_stats()
    p.sort_stats('cumulative').print_stats(10)
    #test = 'TCTG'
    #print revcomp(rcdict, test)
    #testd = {}
    #d = 1
    #testks = genkmers(test, testd, d)
    #print 'testks', testks
    #print 'len testks', len(testks.keys())
    #print 9*nchoosek(3, 2) + nchoosek(3,1)*3

