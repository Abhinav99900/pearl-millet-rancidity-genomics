import csv,collections,statistics,json
R=list(csv.DictReader(open('uni/Coding_alleles.tsv'),delimiter='\t'))
def load(fn):
    D={}
    for l in open(fn):
        c,p,ref,alt,k,n,A,Rr,O,C=l.rstrip('\n').split('\t'); D[(c,int(p),ref,alt)]=(k,int(n),int(A),int(Rr),int(O),int(C))
    return D
X={'GHB-732':(load('bam/self732.out'),load('bam/x538.out')),'GHB-538':(load('bam/self538.out'),load('bam/x732.out'))}
# per-BAM median depth over all checked sites
d732=[v[1] for v in list(X['GHB-732'][0].values())+list(X['GHB-538'][1].values())]
d538=[v[1] for v in list(X['GHB-538'][0].values())+list(X['GHB-732'][1].values())]
M={'732':statistics.median(d732),'538':statistics.median(d538)}
print('median MQ>=20 read depth at checked sites',M)
def altlike(v,ref,alt):
    k,n,A,Rr,O,C=v; L=abs(len(alt)-len(ref))
    a=A+(C if (k!='snv' and L>=15) else 0)
    return (a/n if n else 0),((A+O+C)/n if n else 0)
out=[];cnt=collections.Counter()
for r in R:
    r=dict(r)
    if r['reliable']!='True' or r['status_in_other']=='shared':
        r['read_check']='shared' if r['status_in_other']=='shared' else 'not_reliable'; out.append(r); continue
    me,oth=('732','538') if r['variety']=='GHB-732' else ('538','732')
    key=(r['chrom'],int(r['pos']),r['ref'],r['alt'])
    if key not in X[r['variety']][0] or key not in X[r['variety']][1]:
        r['read_check']='not_checked'; cnt[(r['variety'],r['impact'],'not_checked')]+=1; out.append(r); continue
    s=X[r['variety']][0][key]; o=X[r['variety']][1][key]
    sa,_=altlike(s,r['ref'],r['alt']); oa,onr=altlike(o,r['ref'],r['alt'])
    r.update(self_reads=s[1],self_alt_frac=round(sa,3),other_reads=o[1],other_alt_frac=round(oa,3),other_nonref_frac=round(onr,3))
    if s[1]>2*M[me] or o[1]>2*M[oth]: st='excess_depth'
    elif sa<0.15: st='not_supported_by_own_reads'
    elif o[1]<10: st='other_uncovered'
    elif oa>=0.15: st='present_in_other_reads'
    elif oa<0.05 and onr<0.25: st='specific'
    elif onr>=0.25: st='other_divergent_site'
    else: st='ambiguous'
    r['read_check']=st; cnt[(r['variety'],r['impact'],st)]+=1; out.append(r)
for k in sorted(cnt): print(k,cnt[k])
json.dump(M,open('uni/median_depth.json','w'))
with open('uni/coding_alleles_validated.tsv','w') as f:
    keys=list(out[0].keys())+['self_reads','self_alt_frac','other_reads','other_alt_frac','other_nonref_frac']
    keys=list(dict.fromkeys(keys))
    w=csv.DictWriter(f,keys,delimiter='\t',extrasaction='ignore'); w.writeheader(); [w.writerow(o) for o in out]
