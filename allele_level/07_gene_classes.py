import csv,collections,json,statistics as st
R=list(csv.DictReader(open('uni/coding_alleles_validated.tsv'),delimiter='\t'))
U=list(csv.DictReader(open('uni/Gene_classes.tsv'),delimiter='\t'))
med=json.load(open('uni/median_vcfdp.json'))
coords={l.split()[3]:l.split()[:3] for l in open('uni/universe_regions.txt')}
cat={}
for r in csv.DictReader(open('uni/af3_assign.tsv'),delimiter='\t'):
    c=r['functional_category']
    if r['gene_id'] not in cat or cat[r['gene_id']].startswith('Not'): cat[r['gene_id']]=c
gd=collections.defaultdict(lambda: collections.defaultdict(list)); par=collections.Counter()
spec=collections.defaultdict(lambda: collections.defaultdict(list)); sh=collections.Counter(); rel=collections.Counter()
for r in R:
    g=r['gene']; gd[g][r['variety']].append(int(r['DP']))
    if r['paralog_flag']=='True': par[g]+=1
    if r['reliable']=='True': rel[g]+=1
    if r['read_check']=='specific': spec[g]['LR' if r['variety']=='GHB-732' else 'HR'].append(r)
    if r['read_check']=='shared' and r['variety']=='GHB-732': sh[g]+=1
out=[]
for u in U:
    g=u['gene_id']
    dr={v:st.median(L)/med[v] for v,L in gd[g].items()}
    complex_=par[g]>0 or any(x>2 for x in dr.values())
    lr=len(spec[g]['LR']); hr=len(spec[g]['HR'])
    if complex_: c='Complex locus'
    elif lr and not hr: c='LR-specific'
    elif hr and not lr: c='HR-specific'
    elif lr and hr: c='Both varieties'
    elif sh[g]: c='Shared alleles only'
    else: c='No reliable coding difference'
    lrH=sum(x['impact']=='HIGH' for x in spec[g]['LR']); hrH=sum(x['impact']=='HIGH' for x in spec[g]['HR'])
    out.append(dict(gene_id=g,chrom=coords[g][0],start=coords[g][1],end=coords[g][2],in_732_list=u['in_732_list'],in_538_list=u['in_538_list'],
        functional_category=cat.get(g,'Not assigned (auxiliary/other)'),annotation=u['function'],allele_class=c,
        n_LR_specific=lr,n_LR_specific_HIGH=lrH,n_HR_specific=hr,n_HR_specific_HIGH=hrH,n_shared_reliable=sh[g],
        paralog_flagged_variants=par[g],depth_ratio_732=round(dr.get('GHB-732',float('nan')),2),depth_ratio_538=round(dr.get('GHB-538',float('nan')),2)))
with open('uni/final_gene_classes.tsv','w') as f:
    w=csv.DictWriter(f,list(out[0]),delimiter='\t'); w.writeheader(); [w.writerow(o) for o in out]
C=collections.Counter(o['allele_class'] for o in out); print(C)
print('HIGH: LR-spec genes with HIGH',sum(o['n_LR_specific_HIGH']>0 for o in out if o['allele_class'] in('LR-specific','Both varieties')),
      'HR',sum(o['n_HR_specific_HIGH']>0 for o in out if o['allele_class'] in('HR-specific','Both varieties')))
print('complex: paralog',sum(o['paralog_flagged_variants']>0 for o in out),'depth>2',sum((o['depth_ratio_732']>2) or (o['depth_ratio_538']>2) for o in out))
# list membership vs class
X=collections.Counter(); 
for o in out:
    grp='common' if o['in_732_list']=='yes' and o['in_538_list']=='yes' else ('732-list only' if o['in_732_list']=='yes' else '538-list only')
    X[(grp,o['allele_class'])]+=1
for k in sorted(X): print(k,X[k])
T=collections.Counter((o['functional_category'],o['allele_class']) for o in out)
cats=sorted(set(o['functional_category'] for o in out)); cls=['LR-specific','HR-specific','Both varieties','Shared alleles only','No reliable coding difference','Complex locus']
print('\t'.join(['cat']+cls))
for c in cats: print('\t'.join([c[:30]]+[str(T[(c,k)]) for k in cls]))
CH=collections.Counter((o['chrom'],o['allele_class']) for o in out)
for ch in sorted(set(o['chrom'] for o in out)): print(ch,[CH[(ch,k)] for k in cls])
# counts of validated specific alleles
A=collections.Counter((r['variety'],r['impact']) for r in R if r['read_check']=='specific' )
print(A)
