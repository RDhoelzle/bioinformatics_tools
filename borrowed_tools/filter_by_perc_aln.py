import sys
in_file = sys.argv[1]
cutoff = float(sys.argv[2])
for line in open(in_file):
    if line.startswith('#'): continue
    sline = line.strip().split()
    hmmsize = float(sline[5])
    f = float(sline[15])
    t = float(sline[16])
    l = t - f
    if l/hmmsize>cutoff:
        print sline[0]