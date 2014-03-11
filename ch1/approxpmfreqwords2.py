import re
import numpy as np
import itertools
import math
import cProfile
import pstats

t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
d = 1
k = 4

#t = 'CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC'
#k = 10
#d = 2

#t = 'TAAGGAAAGTGAGAAAGAGAGAAAGCAGAGAGAGAGCAGAGAGTGAGCAGTGAGAGAGTGAGAGAGCAGAAAGCAGAAAGAAAGAGAGAAAGTGAGAGAGAAAGAGAGTGAGAAAGTGAGAAAGAGAGAAAGAGAGAAAGCAGAGATAAGTAAGGTGATAAGGTGATAAGGTGAGAGAGTGAGTGAGCAGAAAGCATAAGGCAGCATAAGTAAGGTGAGTGAGCAGAGATAAGGCAGAAATAAGGAAAGAGAGAGAGTGATAAGGCATAAGGTGAGTGAGAAAGAGAGTGAGCATAAGGAGAGCAGAGAGAGAGAAA '
#k = 10
#d = 3

b = {
'0': 'A', 
'1': 'T', 
'2': 'C', 
'3': 'G'
}

base = 4

def num2base(number,base):
    if number < 0:
        return '-' + num2base(-number, base)
    d, m = divmod(number, base)
    if d > 0:
        return num2base(d, base) + str(m)
    return(str(m))



def genkmers (k, base):
    kmers = []
    numkmers = 4**k

    for i in range(numkmers):
        numform = num2base(i, base)
        numform.zfill(num2base(numkmers, base)-1)
            
        for k, v in b.iteritems():
            numform = numform.replace(k, v)

        kmers.append(numform)

    return kmers


def approxpm(p, t, d):
    indices = []

    for n in range(len(t)):
        #print ('n', n)
        window = t[n:n+len(p)]
        #print(window)
        #print(p)
        #print('wind', len(window))
        if len(window) < len(p):
            #print(len(window), len(p))
            break
        counter = 0
        for o in range(len(p)):
            ##print('o', o)
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

#def hamdist(str1, str2):
   #"""Count the # of differences between equal length strings str1 and str2"""

   #diffs = 0
   #for ch1, ch2 in zip(str1, str2):
       #if ch1 != ch2:
           #diffs += 1
   #return diffs


def freqapproxpm(t, d, k, base):

    freq = 0
    freqkmers = []
    countfreq = []

    kmers = genkmers(k, base)
    for kmer in kmers:
        apmatches = len(approxpm(kmer, t, d))
        if apmatches >= freq:
            freq = apmatches
            freqkmers.append(kmer)
            countfreq.append(apmatches)
    
    maxi = maxindex(countfreq)
    return np.array(freqkmers)[maxi]
    #return freq, freqkmers, countfreq
    
def maxindex(anylist):
    mvlist = np.tile(max(anylist), len(anylist))
    indlist = np.arange(len(anylist))
    return(indlist[mvlist == anylist])
    
def writetofile(anslist, ansfile):
    with open(ansfile, 'w') as f:
        for ans in anslist:
            f.write('{0} '.format(ans))
    

def fast(k, base, t, d):
    freq = 0
    freqkmers = []
    countfreq = []
    
    
    numkmers = 4**k
    for i in range(numkmers):
        numform = num2base(i, base)
        #print(numform)
        #print(len(num2base(numkmers, base))-1)
        lendiff = len(num2base(numkmers, base)) - len(numform)-1
        if lendiff > 0:
            numform = '0'*lendiff + numform
        #print(numform)
        for k, v in b.iteritems():
            numform = numform.replace(k, v)
        
        p = numform
        #print(p)
        apmatches = len(approxpm(p, t, d))
        #print(apmatches)
        if apmatches >= freq:
            freq = apmatches
            freqkmers.append(p)
            countfreq.append(apmatches)
    #print(freqkmers)
    maxi = maxindex(countfreq)
    return np.array(freqkmers)[maxi]

#def fast(k, base, t, d):
    #freq = 0
    #freqkmers = []
    #countfreq = []
    
    
    #numkmers = 4**k
    #for i in range(numkmers):
        #numform = num2base(i, base)
        ##print(numform)
        ##print(len(num2base(numkmers, base))-1)
        #lendiff = len(num2base(numkmers, base)) - len(numform)-1
        #if lendiff > 0:
            #numform = '0'*lendiff + numform
        ##print(numform)
        #for k, v in b.iteritems():
            #numform = numform.replace(k, v)
        
        #p = numform
        ##print(p)
        #apmatches = len(approxpm(p, t, d))
        ##print(apmatches)
        #if apmatches >= freq:
            #freq = apmatches
            #freqkmers.append(p)
            #countfreq.append(apmatches)
    ##print(freqkmers)
    #maxi = maxindex(countfreq)
    #return np.array(freqkmers)[maxi]



#print fast(k, base)
cProfile.run('fast(k, base, t, d)', 'prof_approxpmfreqwords2')
p = pstats.Stats('prof_approxpmfreqwords2')
p.strip_dirs().sort_stats(-1).print_stats()
p.sort_stats('cumulative').print_stats(10)


