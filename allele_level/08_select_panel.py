import json,csv,collections
LEN=dict(LOX2_ORYSJ=870,LOX5_ORYSJ=899,LOX3_ORYSJ=None,LOX4_ORYSJ=877,MGLL_HUMAN=303,OBL1_TOBAC=572,SDP1_ARATH=825,LIP_RHIMI=363,SODC2_MAIZE=151,
 PER15_IPOBA=327,PER16_ARATH=323,PER18_ARATH=329,PER1_ARAHY=316,PER1_HORVU=315,PER1_MAIZE=367,PER1_ORYSJ=326,PER22_ORYSJ=317,PER24_ARATH=350,
 PER2_CUCSA=292,PER2_MAIZE=335,PER43_ARATH=326,PER47_ARATH=325,PER48_ARATH=404,PER52_ARATH=324,PER54_ARATH=358,PER5_VITVI=329,PER70_MAIZE=321,PER7_ARATH=349,PER66_MAIZE=320)
# documented seed/grain-storage function of the ortholog family member (tie-breaker ii)
SEEDFN={'SDP1_ARATH':'SUGAR-DEPENDENT1, major seed-oil TAG lipase (Eastmond 2006)'}
P=json.load(open('uni/panel184.json'))
G={r['gene_id']:r for r in csv.DictReader(open('uni/final_gene_classes.tsv'),delimiter='\t')}
R=list(csv.DictReader(open('uni/coding_alleles_validated.tsv'),delimiter='\t'))
IMP={'HIGH':4,'MODERATE':3}
spec=collections.defaultdict(list)
for r in R:
    if r['read_check']=='specific': spec[r['gene']].append(r)
rows=[]
for p in P:
    g=p['gene_id']; L=spec.get(g,[]); gc=G[g]
    row=dict(family=p['family'],gene_id=g,homolog=p['homolog'],relevance=int(p['relevance']),protein_aa=int(p['protein_aa']),
             hit_aa=LEN.get(p['homolog']),allele_class=gc['allele_class'],n_validated_specific=len(L),
             n_LR=sum(r['variety']=='GHB-732' for r in L),n_HR=sum(r['variety']=='GHB-538' for r in L))
    if L:
        b=max(L,key=lambda r:(IMP[r['impact']],float(r['self_alt_frac'])))
        row.update(best_allele=b['hgvs_p'] or b['hgvs_c'],best_allele_variety=b['variety'],impact=b['impact'],GT=b['GT'],own_read_alt_frac=float(b['self_alt_frac']),
                   other_reads=int(b['other_reads']),other_alt_frac=float(b['other_alt_frac']),score=2*IMP[b['impact']]+2*int(p['relevance']))
    cov=(row['protein_aa']/row['hit_aa']) if row['hit_aa'] else None
    row['model_coverage']=round(cov,2) if cov else None
    if not L: elig='No validated variety-specific coding allele'
    elif gc['allele_class']=='Complex locus': elig='Excluded: complex locus (paralog-like cluster or >2x depth)'
    elif cov is not None and cov<0.7: elig='Excluded: truncated gene model (<70% of best-hit length)'
    else: elig='Eligible'
    row['eligibility']=elig; row['seed_function']=SEEDFN.get(p['homolog'],'')
    rows.append(row)
sel={}
fam=collections.defaultdict(list)
for r in rows:
    if r['eligibility']=='Eligible': fam[r['family']].append(r)
for f,L in fam.items():
    L.sort(key=lambda r:(-r['score'],-r['n_validated_specific'],r['seed_function']=='',-r['own_read_alt_frac']))
    for i,r in enumerate(L): r['rank_eligible']=i+1
    sel[f]=L[0]['gene_id']; L[0]['selected']='YES'
    print('==',f); [print('   ',r['rank_eligible'],r['gene_id'],r['homolog'],r['score'],r['n_validated_specific'],r['best_allele'],r['best_allele_variety'],r['own_read_alt_frac'],r['model_coverage']) for r in L[:5]]
print(sel)
print(collections.Counter((r['family'],r['eligibility']) for r in rows))
json.dump(rows,open('uni/ranking2.json','w'))
