#!/usr/bin/env python3
"""Threshold sensitivity of the allele-level classification and panel selection.

Re-runs steps 03, 06, 07 and 08 in temporary copies of this directory with one threshold changed
at a time (paralog-cluster size, allele-fraction cut-off and window; complex-locus depth multiple).
Requires the SnpEff-annotated VCFs from step 02 in vcf/. Alleles that are reliable only under an
alternative setting have no read-level data and are counted as not validated ('not_checked').
Output: uni/threshold_sensitivity.tsv
"""
import os,shutil,subprocess,json,csv,collections,re,tempfile
SCEN=[('baseline',4,0.35,100,2.0),('K=3',3,0.35,100,2.0),('K=5',5,0.35,100,2.0),('AF<0.30',4,0.30,100,2.0),
      ('AF<0.40',4,0.40,100,2.0),('W=50bp',4,0.35,50,2.0),('W=200bp',4,0.35,200,2.0),('depth>1.5x',4,0.35,100,1.5),('depth>3x',4,0.35,100,3.0)]
here=os.path.dirname(os.path.abspath(__file__)); tmp=tempfile.mkdtemp(); rows=[]
for name,K,AF,W,DM in SCEN:
    d=os.path.join(tmp,re.sub(r'[^A-Za-z0-9]','_',name))
    shutil.copytree(here,d,ignore=shutil.ignore_patterns('out','figs','*.png','16_threshold_sensitivity.py'))
    def patch(fn,pairs):
        s=open(os.path.join(d,fn)).read()
        for a,b in pairs: assert a in s,(fn,a); s=s.replace(a,b)
        open(os.path.join(d,fn),'w').write(s)
    patch('03_allele_table.py',[("r['AF']<0.35:","r['AF']<%s:"%AF),("bisect.bisect_right(ps,p+100)","bisect.bisect_right(ps,p+%d)"%W),("if j-i>=4:","if j-i>=%d:"%K)])
    patch('06_validate_reads.py',[("    s=X[r['variety']][0][key]; o=X[r['variety']][1][key]",
        "    if key not in X[r['variety']][0] or key not in X[r['variety']][1]:\n        r['read_check']='not_checked'; cnt[(r['variety'],r['impact'],'not_checked')]+=1; out.append(r); continue\n    s=X[r['variety']][0][key]; o=X[r['variety']][1][key]")])
    patch('07_gene_classes.py',[("any(x>2 for x in dr.values())","any(x>%s for x in dr.values())"%DM)])
    for step in ['03_allele_table.py','06_validate_reads.py','07_gene_classes.py','08_select_panel.py']:
        subprocess.run(['python3',step],cwd=d,check=True,capture_output=True)
    R=list(csv.DictReader(open(d+'/uni/coding_alleles_validated.tsv'),delimiter='\t'))
    G=list(csv.DictReader(open(d+'/uni/final_gene_classes.tsv'),delimiter='\t'))
    rk=json.load(open(d+'/uni/ranking2.json')); C=collections.Counter(x['allele_class'] for x in G)
    sp=collections.Counter(x['variety'] for x in R if x['read_check']=='specific')
    sel={r['family']:r['gene_id'] for r in rk if r.get('selected')=='YES'}
    rows.append([name,sum(x['paralog_flag']=='True' for x in R),sp['GHB-732'],sp['GHB-538'],sum(x['read_check']=='not_checked' for x in R),
                 C['LR-specific'],C['HR-specific'],C['Both varieties'],C['Complex locus'],sum(r['eligibility']=='Eligible' for r in rk),
                 sel.get('TAG lipase','none'),sel.get('MAG lipase','none'),sel.get('Class III peroxidase','none')])
    print(rows[-1],flush=True)
with open(os.path.join(here,'uni','threshold_sensitivity.tsv'),'w') as f:
    f.write('\t'.join(['setting','paralog_flagged','validated_732','validated_538','not_checked','LR_specific','HR_specific','both','complex','eligible','TAG_lipase','MAG_lipase','class_III_peroxidase'])+'\n')
    for r in rows: f.write('\t'.join(map(str,r))+'\n')
shutil.rmtree(tmp)
