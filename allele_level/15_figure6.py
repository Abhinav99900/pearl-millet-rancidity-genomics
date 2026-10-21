import csv,math,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
plt.rcParams.update({'font.family':'DejaVu Sans'})
L={'chr1':275469374,'chr2':243054348,'chr3':300907625,'chr4':191830058,'chr5':158694342,'chr6':240570548,'chr7':154011256}
G=list(csv.DictReader(open('uni/final_gene_classes.tsv'),delimiter='\t'))
gap=0.012*2*math.pi; tot=sum(L.values()); span=2*math.pi-gap*7
start={};a=math.pi/2
for c in L: start[c]=a; a-= L[c]/tot*span+gap
def ang(c,p): return start[c]-int(p)/tot*span
fig=plt.figure(figsize=(11,11)); ax=fig.add_subplot(111,projection='polar'); ax.set_axis_off(); ax.set_ylim(0,10.6)
cols=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf']
for i,c in enumerate(L):
    th=np.linspace(start[c],start[c]-L[c]/tot*span,300)
    ax.fill_between(th,9.2,9.7,color=cols[i],alpha=0.9,lw=0)
    for k in range(0,11):
        t=start[c]-k/10*L[c]/tot*span; ax.plot([t,t],[9.7,9.85],color='black',lw=0.6)
    mid=start[c]-L[c]/tot*span/2
    ax.text(mid,10.25,f'{c}\n{L[c]/1e6:.0f} Mb',ha='center',va='center',fontsize=12,fontweight='bold')
LEG=[]
tracks=[('LR-specific','#1b7837',8.2,9.0),('HR-specific','#c51b7d',7.2,8.0),('Both varieties','#8c6bb1',6.2,7.0),('Complex locus','#fdae61',5.2,6.0)]
for k,col,r0,r1 in tracks:
    for c in L:
        th=np.linspace(start[c],start[c]-L[c]/tot*span,200); ax.fill_between(th,r0,r1,color='#f2f2f2',lw=0)
    n=0
    for g in G:
        if g['allele_class']==k:
            t=ang(g['chrom'],(int(g['start'])+int(g['end']))//2); ax.plot([t,t],[r0,r1],color=col,lw=1.1); n+=1
    LEG.append((k,n,col))
lab={'Pgl_GLEAN_10017888':'SDP1','Pgl_GLEAN_10036216':'MAGL','Pgl_GLEAN_10033116':'POX (PER52-like)','Pgl_GLEAN_10005615':'LOX3','Pgl_GLEAN_10011213':'PgTAGLip3'}
for g in G:
    if g['gene_id'] in lab:
        t=ang(g['chrom'],g['start']); r={'Pgl_GLEAN_10036216':3.7}.get(g['gene_id'],4.6)
        ax.plot([t,t],[r+0.2,5.2],color='black',lw=1); ax.plot(t,r+0.2,'o',color='black',ms=4)
        dt={'Pgl_GLEAN_10011213':-0.09,'Pgl_GLEAN_10036216':0.06}.get(g['gene_id'],0)
        ax.text(t+dt,r-0.35,lab[g['gene_id']],ha='center',va='center',fontsize=10,fontweight='bold',
                color='#1b7837' if g['allele_class']!='Complex locus' else '#b35806')
for i,(k,n,col) in enumerate(LEG):
    fig.text(0.5,0.49-i*0.032,f'■ {k} ({n})',ha='center',va='center',fontsize=13,fontweight='bold',color=col)
fig.text(0.5,0.525,'Tracks (outer → inner)',ha='center',fontsize=12,color='#333')
plt.savefig('figs/Figure6_circos_allele.png',dpi=300,bbox_inches='tight',facecolor='white')
