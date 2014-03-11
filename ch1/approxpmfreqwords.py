import re
import numpy as np
import itertools
import math

#t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
#d = 1
#k = 4

t = 'CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC'
k = 10
d = 2


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
        lendiff = len(num2base(numkmers, base)) - len(numform)-1
        if lendiff > 0:
            numform = '0'*lendiff + numform
            
        for k, v in b.iteritems():
            numform = numform.replace(k, v)
        numform = numform.replace('0', 'A')
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
        #print('counter', counter)
        
    return indices

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
    
#freqkmers = freqapproxpm(t, d, k, base)
#writetofile(freqkmers, 'approxpmfreqans.txt')

bases = ['C', 'G', 'A', 'T']


#ex = 'CAC'

def genkmers(kmer, poskmers, d, counter):
    
    for i in range(len(kmer)):
        for b in bases:
            #print b
            listk = list(kmer)
            #if listk[i] != b:
            listk[i] = b
            newkmer = ''.join(listk)
            if newkmer not in poskmers:
                poskmers[newkmer] = 0
            
            if counter != 0:
                genkmers(newkmer, poskmers, d, counter-1)

    return poskmers

def findfreq(t, k, d):

    c = {}
    for n in range(len(t)):
        window = t[n:n+k]
        #print(window)
        if len(window) < k:
            break
        poskmerstemp = {}
        poskmers = genkmers(window, poskmerstemp, d, d)
        #if n == 0:
            ##print(poskmers)
        
        for pk in poskmers.keys():
            if pk not in c:
                c[pk] = 1

            apmatches = len(approxpm(pk, t, d))
            #print(inds)
            c[pk] = c[pk]+apmatches

    #print c['']
    
    e = {}
    for key, val in c.iteritems():
        if val not in e:
            e[val] = []
        e[val].append(key)
    print max(e.keys())
    print e[max(e.keys())]

findfreq(t, k, d)
#d = 2
#p = 'GATG'
#isets = itertools.combinations(range(len(p)), d)
#for i in isets:
    #print i

#myans = approxpm(p, t, d)
##print myans[0:5]
##print myans == ans
##print myans
##with open('approxpmans.txt', 'w') as g:
    ##for x in myans:
        ##print x
        ##g.write('{0} '.format(x))

#print len(myans)

#ans = [int(x) for x in ansstr.strip('\n').split(' ')]
#t = t[:20]

#print 'p', p
#print 't', t
#print 'd', d
#print ans

#for n in range(len(p)):
    #s = re.compile(
#p = re.compile('ab*')

#test = t[6:14]
#print test
#s = re.compile('AATCCA[\w][\w]')
#res = s.search(t)
#print res.start()

#isets = itertools.combinations(range(len(p)), d)
#for iset in isets:
    

#check that number of sets matches n choose k
#x = 0
#for i in sets:
    #print i
    #x = x+1
#print x

#print math.factorial(len(p))/(math.factorial(d)*math.factorial((len(p)-d)))



