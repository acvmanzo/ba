with open('ans.txt', 'r') as f:
    with open('revans.txt', 'w') as g:
        g.write('[')
        for l in f:
            g.write('\''+l.replace('\n', '\'') + ',\n')
        g.write(']')
