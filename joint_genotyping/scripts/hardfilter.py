"""Hard filters exactly as run in the original pipeline (bcftools view -i, thesis Appendix Step 10):
SNP  : QD>=1.8 && FS<=80 && SOR<=4.0 && MQ>=40 && (MQRankSum>=-15 || .) && (ReadPosRankSum>=-10 || .)
INDEL: QD>=1.8 && FS<=200 && SOR<=10.0 && (ReadPosRankSum>=-20 || .)
A missing QD/FS/SOR/MQ fails (bcftools semantics); missing rank-sum tests pass.
Only pure SNP and pure INDEL records are retained (mixed, MNP-containing and '*' records drop out)."""
def info(s):
    d={}
    for x in s.split(';'):
        if '=' in x: k,v=x.split('=',1); d[k]=v
    return d
def vtype(ref,alts):
    if '*' in alts: return None
    if len(ref)==1 and all(len(a)==1 for a in alts): return 'SNP'
    if all(len(a)!=len(ref) for a in alts): return 'INDEL'
    return None
def passes(ref,alts,inf):
    t=vtype(ref,alts)
    if t is None: return None
    d=info(inf)
    def v(k):
        try: return float(d[k])
        except (KeyError,ValueError): return None
    req=lambda k,ok: v(k) is not None and ok(v(k))
    opt=lambda k,ok: v(k) is None or ok(v(k))
    if t=='SNP':
        return (req('QD',lambda x:x>=1.8) and req('FS',lambda x:x<=80) and req('SOR',lambda x:x<=4.0) and req('MQ',lambda x:x>=40)
                and opt('MQRankSum',lambda x:x>=-15) and opt('ReadPosRankSum',lambda x:x>=-10))
    return req('QD',lambda x:x>=1.8) and req('FS',lambda x:x<=200) and req('SOR',lambda x:x<=10.0) and opt('ReadPosRankSum',lambda x:x>=-20)
