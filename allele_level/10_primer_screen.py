import csv,json
UP='uni/cds/'
cds={}
for fn in [UP+'RANCIDITY_538_CDS.tsv',UP+'RANCIDITY_732_CDS.tsv']:
    for r in csv.DictReader(open(fn),delimiter='\t'): cds.setdefault(r['gene_id'],r['cds_sequence'].upper())
print('unique CDS',len(cds))
M=json.load(open('primers/masked.json'))
rc=lambda s:s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
def minmm(p,skip):
    best=99
    for g,s in cds.items():
        if g==skip: continue
        for q in (p,rc(p)):
            L=len(q)
            for i in range(len(s)-L+1):
                mm=0
                for a,b in zip(q,s[i:i+L]):
                    if a!=b:
                        mm+=1
                        if mm>=best: break
                if mm<best: best=mm
                if best==0: return 0
    return best
out={}
for name,g in (('SDP1','Pgl_GLEAN_10017888'),('MAGL','Pgl_GLEAN_10036216'),('PER52','Pgl_GLEAN_10033116')):
    s=cds[g]; var=M['var'][g]; res=[]
    for rec in open(f'primers/pb_{name}.txt').read().strip().split(';'):
        n,F,fs,ftm,fgc,R,rs,rtm,rgc,prod=rec.split(','); fs=int(fs); rs=int(rs)
        okF=s[fs-1:fs-1+len(F)]==F; okR=s[rs-len(R):rs]==rc(R); assert okF and okR,(name,n)
        assert rs-fs+1==int(prod)
        bind=set(range(fs,fs+len(F)))|set(range(rs-len(R)+1,rs+1))
        vb=[p for p in var if int(p) in bind]
        amp_indel=[p for p,v in var.items() if fs<=int(p)<=rs and any('indel' in x for x in v)]
        amp_var=[p for p in var if fs<=int(p)<=rs]
        res.append(dict(pair=int(n),F=F,R=R,start=fs,end=rs,FTm=float(ftm),RTm=float(rtm),FGC=float(fgc),RGC=float(rgc),prod=int(prod),
            mmF=minmm(F,g),mmR=minmm(R,g),var_in_primers=len(vb),indel_in_amplicon=len(amp_indel),snv_in_amplicon=len(amp_var)-len(amp_indel)))
    res.sort(key=lambda r:(-(min(r['mmF'],r['mmR'])>=4),r['indel_in_amplicon'],-min(r['mmF'],r['mmR']),r['pair']))
    out[name]=res
    for r in res[:6]: print(name,r)
json.dump(out,open('primers/screen.json','w'))
