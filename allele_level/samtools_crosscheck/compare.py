import re,statistics as st,math,sys,csv
D='samtools_crosscheck/'
al={(r[0],int(r[1])):(r[2],r[3]) for r in (l.split() for l in open('samtools_crosscheck/sites30_alleles.tsv'))}
def parse(fn):
    out={}
    for l in open(fn):
        f=l.rstrip('\n').split('\t'); c,p,dp,b=f[0],int(f[1]),int(f[3]),f[4]
        b=re.sub(r'\^.','',b).replace('$','');b2='';i=0
        while i<len(b):
            if b[i] in '+-':
                m=re.match(r'[+-](\d+)',b[i:]); i+=len(m.group(0))+int(m.group(1)); continue
            b2+=b[i]; i+=1
        out[(c,p)]=(dp,b2.upper())
    return out
mp={'732':parse(D+'mpileupA_732.txt'),'538':parse(D+'mpileupA_538.txt')}
def load(fn): return {(f[0],int(f[1])):(int(f[5]),int(f[6]),int(f[7]),int(f[8])) for f in (l.split('\t') for l in open(fn))}
cust={'732':{**load('bam/self732.out'),**load('bam/x732.out')},'538':{**load('bam/x538.out'),**load('bam/self538.out')}}
own732=set((l.split()[0],int(l.split()[1])) for l in open('bam/sites_732.tsv'))
def corr(x,y):
    mx,my=st.mean(x),st.mean(y); return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
rows=[]
for v in ('732','538'):
    for k,(ref,alt) in sorted(al.items()):
        dp,b=mp[v].get(k,(0,'')); sR=b.count(ref); sA=b.count(alt)
        n,A,R,O=cust[v][k]
        sa=sA/(sR+sA) if sR+sA else 0; ca=A/(R+A) if R+A else 0
        rows.append(dict(bam='GHB-'+v,chrom=k[0],pos=k[1],ref=ref,alt=alt,allele_of='GHB-732' if k in own732 else 'GHB-538',
          samtools_ref=sR,samtools_alt=sA,samtools_AF=round(sa,3),custom_ref=R,custom_alt=A,custom_AF=round(ca,3)))
w=csv.DictWriter(open('samtools_crosscheck/crosscheck_table.tsv','w'),fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
cls=lambda a:'pres' if a>=0.15 else ('abs' if a<0.05 else 'mid')
for sub,name in ((rows,'all'),([r for r in rows if r['bam']=='GHB-732'],'732'),([r for r in rows if r['bam']=='GHB-538'],'538')):
    sa=[r['samtools_AF'] for r in sub]; ca=[r['custom_AF'] for r in sub]
    dep_s=[r['samtools_ref']+r['samtools_alt'] for r in sub]; dep_c=[r['custom_ref']+r['custom_alt'] for r in sub]
    exact=sum(r['samtools_ref']==r['custom_ref'] and r['samtools_alt']==r['custom_alt'] for r in sub)
    within2=sum(abs(a-b)<=2 for a,b in zip(dep_s,dep_c))
    d=[abs(a-b) for a,b in zip(sa,ca)]
    print(name,len(sub),'exact',exact,'depth within2',within2,'AF r',round(corr(sa,ca),4),'mean|dAF|',round(st.mean(d),4),'max',round(max(d),3),'depth r',round(corr(dep_s,dep_c),4),'median depth ratio',round(st.median([b/a for a,b in zip(dep_s,dep_c) if a]),3),'class agree',sum(cls(a)==cls(b) for a,b in zip(sa,ca)))
for r in rows:
    if cls(r['samtools_AF'])!=cls(r['custom_AF']) or abs(r['samtools_AF']-r['custom_AF'])>0.05: print(r)
