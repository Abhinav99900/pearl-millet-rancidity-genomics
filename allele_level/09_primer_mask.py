import csv,json,gzip
CDS={'Pgl_GLEAN_10017888':('chr2','+',[(1467956,1468195),(1468341,1469207),(1469671,1469814),(1469903,1469963),(1470204,1470868)]),
     'Pgl_GLEAN_10036216':('chr6','-',[(234247242,234247360),(234245665,234245703),(234245484,234245554),(234245048,234245112),(234244883,234244972),(234244558,234244665),(234243736,234244152)]),
     'Pgl_GLEAN_10033116':('chr3','-',[(202355236,202355774),(202353506,202353695)])}
seq={}
for fn in ['uni/cds/RANCIDITY_538_CDS.tsv','uni/cds/RANCIDITY_732_CDS.tsv']:
    for r in csv.DictReader(open(fn),delimiter='\t'):
        if r['gene_id'] in CDS:
            s=r['cds_sequence'].upper()
            if r['gene_id'] in seq: assert seq[r['gene_id']]==s, 'CDS differs between files'
            seq[r['gene_id']]=s
def g2c(g,pos):
    ch,st,segs=CDS[g]; off=0
    order=segs if st=='+' else sorted(segs,reverse=True)
    for a,b in order:
        if a<=pos<=b: return off+(pos-a if st=='+' else b-pos)+1
        off+=b-a+1
    return None
# sanity checks with SnpEff HGVS
print(g2c('Pgl_GLEAN_10017888',1470216),'expect 1325'); print(g2c('Pgl_GLEAN_10036216',234243871),'expect 774'); print(g2c('Pgl_GLEAN_10033116',202353636),'expect 599')
for g,(ch,st,segs) in CDS.items(): print(g,len(seq[g]),sum(b-a+1 for a,b in segs))
var={g:{} for g in CDS}
for v,fn in (('732','vcf/u732.filt.vcf.gz'),('538','vcf/u538.filt.vcf.gz')):
    for l in gzip.open(fn,'rt'):
        if l[0]=='#': continue
        f=l.split('\t'); ch=f[0]; p=int(f[1]); ref=f[3]; alts=f[4].split(',')
        for g,(c,st,segs) in CDS.items():
            if c!=ch: continue
            span=range(p,p+len(ref)+(1 if any(len(a)!=len(ref) for a in alts) else 0))
            for q in span:
                i=g2c(g,q)
                if i: var[g].setdefault(i,set()).add(v+(':indel' if any(len(a)!=len(ref) for a in alts) else ''))
out={}
for g in CDS:
    s=list(seq[g])
    for i in var[g]: s[i-1]='N'
    out[g]=''.join(s); print(g,'variant CDS positions',len(var[g]))
json.dump({'seq':seq,'masked':out,'var':{g:{str(k):sorted(v) for k,v in d.items()} for g,d in var.items()}},open('primers/masked.json','w'))
for g in CDS: open(f'primers/{g}_masked.txt','w').write(out[g])
