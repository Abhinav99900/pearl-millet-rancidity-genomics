# PER52-like peroxidase: local primer-candidate design (no Primer-BLAST pair met specificity);
# the selected pair was then checked in Primer-BLAST primer-check mode.
import csv,json,collections,numpy as np,math
UP='uni/cds/'
cds={}
for fn in [UP+'RANCIDITY_538_CDS.tsv',UP+'RANCIDITY_732_CDS.tsv']:
    for r in csv.DictReader(open(fn),delimiter='\t'): cds.setdefault(r['gene_id'],r['cds_sequence'].upper())
g='Pgl_GLEAN_10033116'; s=cds[g]; M=json.load(open('primers/masked.json')); ms=M['masked'][g]; var=M['var'][g]
rc=lambda x:x[::-1].translate(str.maketrans('ACGTN','TGCAN'))
k=8; km=set(s[i:i+k] for i in range(len(s)-k+1)); km|=set(rc(x) for x in km)
sc=collections.Counter({h:sum(t[i:i+k] in km for i in range(len(t)-k+1)) for h,t in cds.items() if h!=g})
top=[h for h,n in sc.most_common(60)]
enc=lambda x:np.frombuffer(x.encode(),dtype=np.uint8)
T=[enc(cds[h]) for h in top]+[enc(rc(cds[h])) for h in top]
def mm(p):
    w=enc(p); b=99
    for t in T:
        if len(t)<len(w): continue
        v=np.lib.stride_tricks.sliding_window_view(t,len(w)); b=min(b,int((v!=w).sum(1).min()))
    return b
# simple NN Tm (SantaLucia 1998 unified, 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM oligo; Owczarzy 2008 Mg correction approx via von Ahsen equivalent)
NN={'AA':(-7.9,-22.2),'TT':(-7.9,-22.2),'AT':(-7.2,-20.4),'TA':(-7.2,-21.3),'CA':(-8.5,-22.7),'TG':(-8.5,-22.7),'GT':(-8.4,-22.4),'AC':(-8.4,-22.4),
    'CT':(-7.8,-21.0),'AG':(-7.8,-21.0),'GA':(-8.2,-22.2),'TC':(-8.2,-22.2),'CG':(-10.6,-27.2),'GC':(-9.8,-24.4),'GG':(-8.0,-19.9),'CC':(-8.0,-19.9)}
def tm(p):
    H=0.2;S=-5.7
    for e in (p[0],p[-1]):
        if e in 'AT': H+=2.2;S+=6.9
    for i in range(len(p)-1): h,s_=NN[p[i:i+2]]; H+=h;S+=s_
    na=0.05+120*math.sqrt(0.0015-0.0006)/1000  # von Ahsen Na-equivalent
    S+=0.368*(len(p)-1)*math.log(na)
    return H*1000/(S+1.987*math.log(50e-9/4))-273.15
cands=[]
for L in range(18,23):
    for i in range(len(s)-L+1):
        f=s[i:i+L]
        if 'N' in ms[i:i+L]: continue
        gc=(f.count('G')+f.count('C'))/L
        if not 0.4<=gc<=0.6: continue
        for strand,p in (('F',f),('R',rc(f))):
            if strand=='F' and p[-1] not in 'GC' and p[-2:] .count('G')+p[-2:].count('C')==0: pass
            t=tm(p)
            if 58<=t<=62: cands.append((strand,i+1,L,p,round(t,2),gc))
print(len(cands))
good=[]
for c in cands:
    m=mm(c[3])
    if m>=4: good.append(c+(m,))
print('unique',len(good))
for c in good: print(c)
json.dump(good,open('primers/per52_cands.json','w'))
