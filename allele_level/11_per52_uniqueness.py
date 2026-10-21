import csv,json,collections,numpy as np
UP='uni/cds/'
cds={}
for fn in [UP+'RANCIDITY_538_CDS.tsv',UP+'RANCIDITY_732_CDS.tsv']:
    for r in csv.DictReader(open(fn),delimiter='\t'): cds.setdefault(r['gene_id'],r['cds_sequence'].upper())
g='Pgl_GLEAN_10033116'; s=cds[g]
rc=lambda x:x[::-1].translate(str.maketrans('ACGTN','TGCAN'))
k=10; km=set(s[i:i+k] for i in range(len(s)-k+1)); km|=set(rc(x) for x in km)
sc=collections.Counter()
for h,t in cds.items():
    if h==g: continue
    sc[h]=sum(t[i:i+k] in km for i in range(len(t)-k+1))
top=[h for h,n in sc.most_common(15)]; print(sc.most_common(15))
enc=lambda x:np.frombuffer(x.encode(),dtype=np.uint8)
T=[enc(cds[h]) for h in top]+[enc(rc(cds[h])) for h in top]
W=20; res=[]
for i in range(len(s)-W+1):
    w=enc(s[i:i+W]); best=99
    for t in T:
        if len(t)<W: continue
        v=np.lib.stride_tricks.sliding_window_view(t,W)
        best=min(best,int((v!=w).sum(1).min()))
    res.append(best)
line=''.join(str(min(x,9)) for x in res)
for j in range(0,len(line),60): print(j+1,line[j:j+60])
json.dump(res,open('primers/per52_uniq.json','w'))
