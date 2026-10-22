#!/usr/bin/env python3
"""Allele table from the JOINT GHB-732 + GHB-538 calls (GATK GenotypeGVCFs on both gVCFs).
Same rules as 03_allele_table.py (reliable genotype, paralog-like cluster), but the status of an
allele in the other variety is read from the other sample's genotype at the same joint record:
  shared           other sample carries the allele, record PASS, other genotype reliable
  present_lowconf  other sample carries the allele but record failed filters or genotype not reliable
  absent_confirmed other sample does not carry the allele and has a called genotype with DP>=10
  absent_unknown_coverage  otherwise (no call or DP<10 in the other sample)
Input : vcf/joint.regions.ann.vcf  Output: uni/Coding_alleles.tsv
Optional args: K AF W (paralog cluster size, allele-fraction cut-off, window) for sensitivity runs."""
import sys,bisect,collections,csv
K,AFC,W=(int(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3])) if len(sys.argv)>3 else (4,0.35,100)
REG={}
for l in open('uni/universe_regions.txt'):
    if l.strip(): c_,s_,e_,g_=l.split()[:4]; REG[g_]=(int(s_)+100,int(e_)-100)
genes=set(REG); RANK={'HIGH':0,'MODERATE':1}
SAMPLE={'GHB-732':'SO_14026_A1','GHB-538':'538'}
recs=[]
for l in open('vcf/joint.regions.ann.vcf'):
    if l.startswith('##'): continue
    f=l.rstrip('\n').split('\t')
    if l[0]=='#': col={s:f.index(s) for s in SAMPLE.values()}; continue
    recs.append(f)
def geno(f,sname):
    fmt=f[8].split(':'); g=dict(zip(fmt,f[col[sname]].split(':')))
    gtraw=g.get('GT','./.'); gt=gtraw.replace('|','/')
    ad=[int(x) if x.isdigit() else 0 for x in g.get('AD','').split(',')] if g.get('AD') not in (None,'','.') else []
    dp=int(g['DP']) if g.get('DP','.').isdigit() else sum(ad)
    gq=int(g['GQ']) if g.get('GQ','.').isdigit() else 0
    return gtraw,gt.split('/'),ad,dp,gq
def rows_for(v):
    out=[]
    for f in recs:
        c,p,ref,alts,info=f[0],int(f[1]),f[3],f[4].split(','),f[7]
        gtraw,calls,ad,dp,gq=geno(f,SAMPLE[v])
        ann=[a.split('|') for a in info.split('ANN=')[1].split(';')[0].split(',')] if 'ANN=' in info else []
        for i,alt in enumerate(alts,1):
            if alt=='*': continue
            tot=sum(ad); af=(ad[i]/tot) if tot and i<len(ad) else 0.0
            best={}; coding=set()
            for a in ann:
                if a[0]!=alt or a[3] not in genes: continue
                if a[2] in ('HIGH','MODERATE','LOW'): coding.add(a[3])
                if a[2] not in RANK: continue
                if a[3] not in best or RANK[a[2]]<RANK[best[a[3]][2]]: best[a[3]]=a
            out.append(dict(rec=id(f),idx=i,chrom=c,pos=p,ref=ref,alt=alt,GT=gtraw,DP=dp,GQ=gq,AF=round(af,3),hom=calls.count(str(i))==2,
                            called=str(i) in calls,nocall='.' in calls,PASS=f[6]=='PASS',ann=best,coding=coding))
    return out
def paralog(rows):
    low=collections.defaultdict(set)
    for r in rows:
        if r['PASS'] and r['called'] and not r['hom'] and r['AF']<AFC:
            for g in r['coding']: low[(g,r['chrom'])].add(r['pos'])
    flag=set()
    for (g,c),ps in low.items():
        ps=sorted(ps)
        for i,p in enumerate(ps):
            j=bisect.bisect_right(ps,p+W)
            if j-i>=K: flag.update((g,c,x) for x in ps[i:j])
    return flag
def reliable(r,flagged,g):
    if (g,r['chrom'],r['pos']) in flagged: return False
    if r['DP']<10 or r['GQ']<20: return False
    return r['AF']>=0.85 if r['hom'] else 0.25<=r['AF']<=0.75
V={}
for v in SAMPLE:
    rows=rows_for(v); V[v]=dict(rows=rows,flag=paralog(rows),by={(r['rec'],r['idx']):r for r in rows})
out=[]
for v,o in (('GHB-732','GHB-538'),('GHB-538','GHB-732')):
    me=V[v]; ot=V[o]
    for r in me['rows']:
        if not (r['PASS'] and r['called']): continue
        x=ot['by'][(r['rec'],r['idx'])]
        for g,a in sorted(r['ann'].items()):
            rel=reliable(r,me['flag'],g)
            if not rel: st='not_reliable'
            elif x['called']: st='shared' if reliable(x,ot['flag'],g) else 'present_lowconf'
            elif not x['nocall'] and x['DP']>=10: st='absent_confirmed'
            else: st='absent_unknown_coverage'
            out.append([v,g,r['chrom'],r['pos'],r['ref'],r['alt'],a[2],a[1],a[9],a[10],r['GT'],r['DP'],r['GQ'],r['AF'],r['hom'],(g,r['chrom'],r['pos']) in me['flag'],rel,st,
                        x['GT'],x['DP'],x['GQ'],x['AF']])
with open('uni/Coding_alleles.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(['variety','gene','chrom','pos','ref','alt','impact','effect','hgvs_c','hgvs_p','GT','DP','GQ','AF','hom','paralog_flag','reliable','status_in_other',
                'other_GT','other_DP','other_GQ','other_AF'])
    w.writerows(out)
print(collections.Counter((r[0],r[17]) for r in out))
