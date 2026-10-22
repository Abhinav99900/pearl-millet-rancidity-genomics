#!/usr/bin/env python3
"""Threshold sensitivity on the joint calls (same 9 settings as 16_threshold_sensitivity.py)."""
import os,shutil,subprocess,json,csv,collections,re,tempfile
SCEN=[('baseline',4,0.35,100,2.0),('K=3',3,0.35,100,2.0),('K=5',5,0.35,100,2.0),('AF<0.30',4,0.30,100,2.0),
      ('AF<0.40',4,0.40,100,2.0),('W=50bp',4,0.35,50,2.0),('W=200bp',4,0.35,200,2.0),('depth>1.5x',4,0.35,100,1.5),('depth>3x',4,0.35,100,3.0)]
here=os.path.dirname(os.path.abspath(__file__)); tmp=tempfile.mkdtemp(); rows=[]
for name,K,AF,W,DM in SCEN:
    d=os.path.join(tmp,re.sub(r'[^A-Za-z0-9]','_',name))
    shutil.copytree(here,d,ignore=shutil.ignore_patterns('out','figs','old','J16_sensitivity_joint.py'))
    s=open(d+'/07_gene_classes.py').read(); a="any(x>2 for x in dr.values())"; assert a in s
    open(d+'/07_gene_classes.py','w').write(s.replace(a,"any(x>%s for x in dr.values())"%DM))
    subprocess.run(['python3','J03_allele_table_joint.py',str(K),str(AF),str(W)],cwd=d,check=True,capture_output=True)
    for step in ['06_validate_reads.py','07_gene_classes.py','08_select_panel.py']:
        subprocess.run(['python3',step],cwd=d,check=True,capture_output=True)
    R=list(csv.DictReader(open(d+'/uni/coding_alleles_validated.tsv'),delimiter='\t'))
    G=list(csv.DictReader(open(d+'/uni/final_gene_classes.tsv'),delimiter='\t'))
    rk=json.load(open(d+'/uni/ranking2.json')); C=collections.Counter(x['allele_class'] for x in G)
    sp=collections.Counter(x['variety'] for x in R if x['read_check']=='specific')
    sel={r['family']:r['gene_id'] for r in rk if r.get('selected')=='YES'}
    rows.append([name,sum(x['paralog_flag']=='True' for x in R),sp['GHB-732'],sp['GHB-538'],sum(x['read_check']=='not_checked' for x in R),
                 C['LR-specific'],C['HR-specific'],C['Both varieties'],C['Complex locus'],sum(r['eligibility']=='Eligible' for r in rk),
                 sel.get('TAG lipase','none'),sel.get('MAG lipase','none'),sel.get('Class III peroxidase','none')])
    print('\t'.join(map(str,rows[-1])),flush=True)
with open(os.path.join(here,'uni','threshold_sensitivity_joint.tsv'),'w') as f:
    f.write('\t'.join(['setting','paralog_flagged','validated_732','validated_538','not_checked','LR_specific','HR_specific','both','complex','eligible','TAG_lipase','MAG_lipase','class_III_peroxidase'])+'\n')
    for r in rows: f.write('\t'.join(map(str,r))+'\n')
shutil.rmtree(tmp)
