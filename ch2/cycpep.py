import numpy as np

MASSTABLEFILE = 'integer_mass_table.txt'



def nsubpep(n):
    return n*(n-1)

def load_masstable(masstablefile):
    d = {}
    with open(masstablefile, 'r') as f:
        for l in f:
            elist = l.strip('\n').split(' ')
            d[elist[0]] = int(elist[1])
    return d

def find_mass(subpep, dictmasstable):
    subpepmass = 0
    for sp in subpep:
        subpepmass = subpepmass + dictmasstable[sp]
    return subpepmass

def gencyclospectrum(pep, masstablefile):
    
    dmt = load_masstable(masstablefile)
    fullpep = pep+pep
    
    cspect = [0]
    
    subpeplengths = range(1, len(pep))
   
    for spl in subpeplengths:
        for n, p in enumerate(pep):
            subpep = fullpep[n:n+spl]            
            #print subpep
            subpepmass = find_mass(subpep, dmt)
            cspect.append(subpepmass)
            
    cspect.append(find_mass(pep, dmt))
    return sorted(cspect)



def pep_with_mass(m, masstablefile):
    
    aamasses = sorted(list(set(load_masstable(masstablefile).values())), 
    reverse=True)
    n_aamasses = len(aamasses)
    
    print aamasses
    
    for aa in aamasses:
        
        newm = m-aa
        
    
pep_with_mass(1024, MASSTABLEFILE)
    
#PEP = 'VQGVLLMGMRAGD'
##ANS = '0 71 71 99 101 103 113 113 114 128 128 131 147 163 170 172 184 199 215 227 227 231 244 259 260 266 271 286 298 298 310 312 328 330 330 372 385 391 394 399 399 399 401 413 423 426 443 443 470 493 498 502 513 519 526 527 541 554 556 557 564 569 590 598 616 626 640 654 657 658 665 670 682 697 697 703 711 729 729 753 753 771 779 785 785 800 812 817 824 825 828 842 856 866 884 892 913 918 925 926 928 941 955 956 963 969 980 984 989 1012 1039 1039 1056 1059 1069 1081 1083 1083 1083 1088 1091 1097 1110 1152 1152 1154 1170 1172 1184 1184 1196 1211 1216 1222 1223 1238 1251 1255 1255 1267 1283 1298 1310 1312 1319 1335 1351 1354 1354 1368 1369 1369 1379 1381 1383 1411 1411 1482'

#myans = ' '.join([str(x) for x in gencyclospectrum(PEP, MASSTABLEFILE)])
#print myans
#print myans == ANS

