import csv,collections,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11})
R=list(csv.DictReader(open('uni/coding_alleles_validated.tsv'),delimiter='\t'))
G=list(csv.DictReader(open('uni/final_gene_classes.tsv'),delimiter='\t'))
def cat(r):
    rc=r['read_check']
    if rc=='not_reliable': return 'Paralog-like cluster' if r['paralog_flag']=='True' else 'Failed genotype QC'
    return {'shared':'Same allele reliable in both','present_in_other_reads':'Allele present in other variety\'s reads','excess_depth':'Excess read depth (>2× median)',
            'specific':'Validated variety-specific'}.get(rc,'Unresolved at read level')
order=['Validated variety-specific','Same allele reliable in both','Allele present in other variety\'s reads','Unresolved at read level','Excess read depth (>2× median)','Paralog-like cluster','Failed genotype QC']
cols=['#1b7837','#4575b4','#91bfdb','#bdbdbd','#fdae61','#d73027','#7f7f7f']
C={v:collections.Counter(cat(r) for r in R if r['variety']==v) for v in ('GHB-732','GHB-538')}
fig,(a,b)=plt.subplots(1,2,figsize=(16,6.4),gridspec_kw={'width_ratios':[1.05,1]})
labs=['GHB-732 (LR)','GHB-538 (HR)']
for i,v in enumerate(['GHB-732','GHB-538']):
    left=0; tot=sum(C[v].values())
    for k,cl in zip(order,cols):
        n=C[v][k]; a.barh(i,n,left=left,color=cl,edgecolor='black',lw=0.8,height=0.6,label=k if i==0 else None)
        if n/tot>0.045: a.text(left+n/2,i,f'{n:,}',ha='center',va='center',fontsize=10,fontweight='bold',color='white' if k not in("Allele present in other variety's reads",'Unresolved at read level','Excess read depth (>2× median)') else 'black')
        left+=n
    a.text(left+60,i,f'{tot:,}',va='center',fontsize=11,fontweight='bold')
a.set_yticks([0,1]); a.set_yticklabels(labs,fontweight='bold'); a.invert_yaxis(); a.set_xlim(0,6300)
a.set_xlabel('HIGH- and MODERATE-impact coding alleles in the combined set of 1,370 candidate genes',fontweight='bold')
a.legend(loc='upper center',bbox_to_anchor=(0.5,-0.16),ncol=2,frameon=False,fontsize=10)
a.set_title('a  Allele-level quality control and read-level validation',loc='left',fontweight='bold')
for s in ('top','right'): a.spines[s].set_visible(False)
gcl=['LR-specific','HR-specific','Both varieties','Shared alleles only','No reliable coding difference','Complex locus']
gcol=['#1b7837','#c51b7d','#8c6bb1','#4575b4','#d9d9d9','#fdae61']
CH=collections.Counter((g['chrom'],g['allele_class']) for g in G); chs=[f'chr{i}' for i in range(1,8)]
bottom=[0]*7
for k,cl in zip(gcl,gcol):
    vals=[CH[(c,k)] for c in chs]
    b.bar(chs,vals,bottom=bottom,color=cl,edgecolor='black',lw=0.8,width=0.7,label=f"{k} ({sum(vals)})")
    for j,(v,bt) in enumerate(zip(vals,bottom)):
        if v>=12: b.text(j,bt+v/2,str(v),ha='center',va='center',fontsize=9,fontweight='bold',color='white' if k in('LR-specific','HR-specific','Both varieties','Shared alleles only') else 'black')
    bottom=[x+y for x,y in zip(bottom,vals)]
for j,t in enumerate(bottom): b.text(j,t+4,str(t),ha='center',va='bottom',fontweight='bold')
b.set_ylabel('Candidate genes',fontweight='bold'); b.set_xlabel('Chromosome',fontweight='bold'); b.set_ylim(0,max(bottom)*1.12)
b.legend(loc='upper center',bbox_to_anchor=(0.5,-0.16),ncol=2,frameon=False,fontsize=10)
b.set_title('b  Gene-level allele classes (n = 1,370)',loc='left',fontweight='bold')
for s in ('top','right'): b.spines[s].set_visible(False)
b.grid(axis='y',ls='--',alpha=0.35); b.set_axisbelow(True)
plt.tight_layout(); plt.savefig('figs/Figure5_allele_level.png',dpi=300,bbox_inches='tight',facecolor='white')
print({v:dict(C[v]) for v in C})
