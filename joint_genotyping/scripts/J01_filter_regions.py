#!/usr/bin/env python3
"""Joint calls: keep records inside the 1,370 gene models (mRNA +/- 100 bp) and set FILTER
with the original SNP/INDEL hard filters (hardfilter.py). Records that are not pure SNP/INDEL get FILTER=TYPE.
Usage: python3 J01_filter_regions.py joint_raw.vcf.gz uni/universe_regions.txt vcf/joint.regions.vcf"""
import sys,gzip,collections,bisect
sys.path.insert(0,__import__('os').path.dirname(__file__) or '.')
from hardfilter import passes
src,regf,out=sys.argv[1:4]
R=collections.defaultdict(list)
for l in open(regf):
    if l.strip(): c,s,e=l.split()[:3]; R[c].append((int(s),int(e)))
for c in R: R[c].sort()
def inside(c,p):
    L=R.get(c,[]); i=bisect.bisect_right(L,(p,10**12))-1
    while i>=0 and L[i][0]<=p:
        if p<=L[i][1]: return True
        i-=1
    return False
cnt=collections.Counter()
with gzip.open(src,'rt') as fi, open(out,'w') as fo:
    for l in fi:
        if l.startswith('##'): fo.write(l); continue
        if l.startswith('#'):
            fo.write('##FILTER=<ID=HardFilter,Description="Failed original SNP/INDEL hard filter">\n##FILTER=<ID=TYPE,Description="Mixed/MNP/spanning-deletion record (not retained in original pipeline)">\n'); fo.write(l); continue
        f=l.split('\t')
        if not inside(f[0],int(f[1])): cnt['outside']+=1; continue
        p=passes(f[3],f[4].split(','),f[7])
        f[6]='PASS' if p else ('TYPE' if p is None else 'HardFilter'); cnt[f[6]]+=1
        fo.write('\t'.join(f))
print(dict(cnt))
