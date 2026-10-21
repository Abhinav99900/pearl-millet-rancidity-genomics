#!/usr/bin/env python3
"""Write the site lists for read-level checking: reliable alleles not called identically in both varieties."""
import csv,sys
tab=sys.argv[1] if len(sys.argv)>1 else 'uni/Coding_alleles.tsv'
S={'GHB-732':[],'GHB-538':[]}
for r in csv.DictReader(open(tab),delimiter='\t'):
    if r['reliable']=='True' and r['status_in_other']!='shared':
        k=(r['chrom'],r['pos'],r['ref'],r['alt'])
        if k not in S[r['variety']]: S[r['variety']].append(k)
for v,fn in (('GHB-732','bam/sites_732.tsv'),('GHB-538','bam/sites_538.tsv')):
    with open(fn,'w') as f:
        for k in S[v]: f.write('\t'.join(k)+'\n')
    print(fn,len(S[v]))
