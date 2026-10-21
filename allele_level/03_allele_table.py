#!/usr/bin/env python3
"""Stage 1 of the allele-level comparison: build the coding-allele table.

Usage:  python3 03_allele_table.py            (run from the allele_level/ directory)
Inputs: vcf/u{732,538}.{filt,raw}.ann.vcf  SnpEff 5.3-annotated VCF subsets of the
        1,370 candidate gene models (mRNA +/- 100 bp; made by 01_ and 02_)
        vcf/u538.gvcf_blocks_at_732_calls.tsv.gz  GHB-538 gVCF reference blocks
        uni/universe_regions.txt                   gene models (mRNA +/- 100 bp)
Output: uni/Coding_alleles.tsv  one row per HIGH/MODERATE coding allele per variety
"""
import sys,gzip,bisect,collections,csv
D='vcf'
REG={}
for l in open('uni/universe_regions.txt'):
    if l.strip():
        c_,s_,e_,g_=l.split()[:4]; REG[g_]=(int(s_)+100,int(e_)-100)   # file holds mRNA +/- 100 bp
genes=set(REG)
RANK={'HIGH':0,'MODERATE':1}
def parse(fn):
    out=[]
    for l in open(fn):
        if l[0]=='#': continue
        f=l.rstrip('\n').split('\t')
        c,p,ref,alts,info,fmt,smp=f[0],int(f[1]),f[3],f[4].split(','),f[7],f[8].split(':'),f[9].split(':')
        g=dict(zip(fmt,smp)); gtraw=g.get('GT','./.'); gt=gtraw.replace('|','/')
        ad=[int(x) if x.isdigit() else 0 for x in g.get('AD','').split(',')] if g.get('AD') else []
        dp=int(g['DP']) if g.get('DP','.').isdigit() else sum(ad)
        gq=int(g['GQ']) if g.get('GQ','.').isdigit() else 0
        ann=[a.split('|') for a in info.split('ANN=')[1].split(';')[0].split(',')] if 'ANN=' in info else []
        for i,alt in enumerate(alts,1):
            if alt=='*': continue
            tot=sum(ad); af=(ad[i]/tot) if tot and i<len(ad) else 0.0
            calls=gt.split('/'); hom=calls.count(str(i))==2
            best={}; coding=set()
            for a in ann:
                if a[0]!=alt or a[3] not in genes: continue
                if a[2] in ('HIGH','MODERATE','LOW'): coding.add(a[3])
                if a[2] not in RANK: continue
                if a[3] not in best or RANK[a[2]]<RANK[best[a[3]][2]]: best[a[3]]=a
            out.append(dict(chrom=c,pos=p,ref=ref,alt=alt,GT=gtraw,DP=dp,GQ=gq,AF=round(af,3),hom=hom,called=str(i) in calls,ann=best,coding=coding))
    return out
def paralog(rows):
    """Flag coding alleles lying in a cluster of >=4 heterozygous coding calls
    (SnpEff HIGH/MODERATE/LOW) with allele fraction <0.35 spanning <=100 bp in one gene."""
    low=collections.defaultdict(set)
    for r in rows:
        if r['called'] and not r['hom'] and r['AF']<0.35:
            for g in r['coding']: low[(g,r['chrom'])].add(r['pos'])
    flag=set()
    for (g,c),ps in low.items():
        ps=sorted(ps)
        for i,p in enumerate(ps):
            j=bisect.bisect_right(ps,p+100)
            if j-i>=4: flag.update((g,c,x) for x in ps[i:j])
    return flag
def reliable(r,flagged,g):
    if (g,r['chrom'],r['pos']) in flagged: return False
    if r['DP']<10 or r['GQ']<20: return False
    return r['AF']>=0.85 if r['hom'] else 0.25<=r['AF']<=0.75
V={}
for v,tag in (('GHB-732','u732'),('GHB-538','u538')):
    filt=parse(f'{D}/{tag}.filt.ann.vcf'); raw=parse(f'{D}/{tag}.raw.ann.vcf')
    V[v]=dict(filt=filt,flag=paralog(filt),raw={(r['chrom'],r['pos'],r['ref'],r['alt']):r for r in raw+filt if r['called']},
              callpos=collections.defaultdict(list))
    for r in filt:
        if r['called'] and r['DP']>=10: V[v]['callpos'][r['chrom']].append(r['pos'])
    for c in V[v]['callpos']: V[v]['callpos'][c].sort()
# GHB-538 gVCF reference blocks (coverage evidence for GHB-732 alleles)
blocks=collections.defaultdict(list)
with gzip.open(f'{D}/u538.gvcf_blocks_at_732_calls.tsv.gz','rt') as fh:
    for l in fh:
        f=l.rstrip('\n').split('\t')
        if len(f)<6 or not f[4].startswith('END='): continue
        dp=f[5].split(':')                      # GT:DP:GQ:MIN_DP:PL
        md=int(dp[1]) if len(dp)>1 and dp[1].isdigit() else 0
        blocks[f[0]].append((int(f[1]),int(f[4][4:]),md))
for c in blocks: blocks[c].sort()
def ref_cov_538(c,p):
    b=blocks.get(c,[]); i=bisect.bisect_right(b,(p,10**12,0))-1
    return i>=0 and b[i][0]<=p<=b[i][1] and b[i][2]>=10
def near_732(c,p):
    """A filtered GHB-732 call with DP>=10 lies within 300 bp (coverage proxy; no GHB-732 gVCF)."""
    ps=V['GHB-732']['callpos'].get(c,[]); i=bisect.bisect_left(ps,p-300)
    return i<len(ps) and ps[i]<=p+300
out=[]
for v,o in (('GHB-732','GHB-538'),('GHB-538','GHB-732')):
    me=V[v]; ot=V[o]
    for r in me['filt']:
        if not r['called']: continue
        for g,a in sorted(r['ann'].items()):
            rel=reliable(r,me['flag'],g); key=(r['chrom'],r['pos'],r['ref'],r['alt'])
            if not rel: st='not_reliable'
            elif key in ot['raw']:
                x=ot['raw'][key]
                st='shared' if any(x is y for y in ot['filt']) and reliable(x,ot['flag'],g) else 'present_lowconf'
            elif v=='GHB-732': st='absent_confirmed' if ref_cov_538(r['chrom'],r['pos']) else 'absent_unknown_coverage'
            else: st='absent_called_nearby' if near_732(r['chrom'],r['pos']) else 'absent_unknown_coverage'
            out.append([v,g,r['chrom'],r['pos'],r['ref'],r['alt'],a[2],a[1],a[9],a[10],r['GT'],r['DP'],r['GQ'],r['AF'],r['hom'],(g,r['chrom'],r['pos']) in me['flag'],rel,st])
with open('uni/Coding_alleles.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(['variety','gene','chrom','pos','ref','alt','impact','effect','hgvs_c','hgvs_p','GT','DP','GQ','AF','hom','paralog_flag','reliable','status_in_other'])
    w.writerows(out)
print(collections.Counter(r[0] for r in out))
