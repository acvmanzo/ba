from approxpmfreqwords3 import *
# Ran skew.py: [3764856, 3764858]

with open('Salmonella_enterica.fasta', 'r') as f:
    f.next()
    g = ''
    for l in f:
        g = g + l.rstrip('\n')
        
oric = g[3764856-500: 3764856+500]
k = 9
d = 1
#ans = ['CTTTTTTAT', 'CGGATCATC', 'GATGATCCG', 'ATAAAAAAG'] for d =2
#ans = ['TTTTAAAAA', 'TTTTTAAAA'] for d = 3
#ans = ['TTATCCACA', 'TGTGGATAA'] for d = 1
#from the internet: TTATC[CA]A[CA]A, pretty good!!

cProfile.run('findfreqrc(oric, k, d, rcdict)', 'prof_approxpmfreqwords3')
p = pstats.Stats('prof_approxpmfreqwords3')
p.strip_dirs().sort_stats(-1).print_stats()
p.sort_stats('cumulative').print_stats(10)
