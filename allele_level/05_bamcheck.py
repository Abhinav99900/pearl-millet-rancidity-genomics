import struct,zlib,sys,collections
def bgzf_block(f,coff):
    f.seek(coff); h=f.read(18)
    if len(h)<18: return None,0
    bsize=struct.unpack('<H',h[16:18])[0]
    data=f.read(bsize-17)
    return zlib.decompress(data[:-8],-15),bsize+1
class Reader:
    def __init__(s,path):
        s.f=open(path,'rb'); s.coff=0; s.buf=b''; s.uoff=0
    def seek(s,voff):
        s.coff=voff>>16; s.buf,s.bs=bgzf_block(s.f,s.coff); s.uoff=voff&0xffff
    def read(s,n):
        out=b''
        while len(out)<n:
            if s.uoff>=len(s.buf):
                s.coff+=s.bs; s.buf,s.bs=bgzf_block(s.f,s.coff); s.uoff=0
                if s.buf is None: return out
            take=min(n-len(out),len(s.buf)-s.uoff); out+=s.buf[s.uoff:s.uoff+take]; s.uoff+=take
        return out
def header(r):
    r.seek(0); assert r.read(4)==b'BAM\1'
    lt=struct.unpack('<i',r.read(4))[0]; r.read(lt)
    n=struct.unpack('<i',r.read(4))[0]; names=[]
    for i in range(n):
        l=struct.unpack('<i',r.read(4))[0]; names.append(r.read(l)[:-1].decode()); r.read(4)
    return names
def reg2bins(beg,end):
    L=[0]; end-=1
    for sh,off in ((26,1),(23,9),(20,73),(17,585),(14,4681)):
        L+=range(off+(beg>>sh),off+(end>>sh)+1)
    return L
def bai_chunks_UNUSED(bai,tid,beg,end):
    f=open(bai,'rb'); assert f.read(4)==b'BAI\1'
    nref=struct.unpack('<i',f.read(4))[0]
    for i in range(nref):
        nb=struct.unpack('<i',f.read(4))[0]; bins={}
        for j in range(nb):
            b,nc=struct.unpack('<Ii',f.read(8)); ch=[struct.unpack('<QQ',f.read(16)) for k in range(nc)]
            if i==tid: bins[b]=ch
        ni=struct.unpack('<i',f.read(4))[0]; lin=f.read(8*ni)
        if i==tid:
            lin=struct.unpack('<%dQ'%ni,lin); minoff=lin[min(beg>>14,ni-1)] if ni else 0
            ch=[c for b in reg2bins(beg,end) for c in bins.get(b,[]) if c[1]>minoff]
            return sorted(ch)
SEQ='=ACMGRSVTWYHKDBN'
def reads(bam,bai,chrom,beg,end):
    r=Reader(bam); names=header(r); tid=names.index(chrom)
    seen=set()
    for cs,ce in bai_chunks(bai,tid,beg,end):
        r.seek(cs)
        while True:
            if (r.coff<<16|r.uoff)>=ce: break
            bs=r.read(4)
            if len(bs)<4: break
            rec=r.read(struct.unpack('<i',bs)[0])
            rid,pos,lrn,mq,bin_,ncig,flag,lseq=struct.unpack('<iiBBHHHi',rec[:20])
            if rid!=tid or pos>end: break
            if flag & 0x704: continue
            p=32; qn=rec[p:p+lrn]; p+=lrn
            if (qn,flag) in seen: continue
            cig=struct.unpack('<%dI'%ncig,rec[p:p+4*ncig]); p+=4*ncig
            sq=''.join(SEQ[(rec[p+i//2]>>(4*(1-i%2)))&15] for i in range(lseq))
            yield pos,mq,flag,[(c>>4,'MIDNSHP=X'[c&15]) for c in cig],sq
            seen.add((qn,flag))

def load_bai(bai):
    f=open(bai,'rb'); assert f.read(4)==b'BAI\1'
    nref=struct.unpack('<i',f.read(4))[0]; IDX=[]
    for i in range(nref):
        nb=struct.unpack('<i',f.read(4))[0]; bins={}
        for j in range(nb):
            b,nc=struct.unpack('<Ii',f.read(8)); raw=f.read(16*nc)
            bins[b]=[struct.unpack_from('<QQ',raw,16*k) for k in range(nc)]
        ni=struct.unpack('<i',f.read(4))[0]; lin=struct.unpack('<%dQ'%ni,f.read(8*ni))
        IDX.append((bins,lin))
    return IDX
def chunks(IDX,tid,beg,end):
    bins,lin=IDX[tid]; ni=len(lin); minoff=lin[min(beg>>14,ni-1)] if ni else 0
    return sorted(c for b in reg2bins(beg,end) for c in bins.get(b,[]) if c[1]>minoff)
def fetch(r,IDX,tid,beg,end):
    seen=set()
    for cs,ce in chunks(IDX,tid,beg,end):
        r.seek(cs)
        while (r.coff<<16|r.uoff)<ce:
            bs=r.read(4)
            if len(bs)<4: break
            rec=r.read(struct.unpack('<i',bs)[0])
            rid,pos,lrn,mq,bin_,ncig,flag,lseq=struct.unpack('<iiBBHHHi',rec[:20])
            if rid!=tid or pos>end: break
            if flag & 0x704 or mq<20: continue
            p=32; qn=rec[p:p+lrn]; p+=lrn
            if (qn,flag) in seen: continue
            seen.add((qn,flag))
            cig=[(c>>4,'MIDNSHP=X'[c&15]) for c in struct.unpack('<%dI'%ncig,rec[p:p+4*ncig])]; p+=4*ncig
            sq=''.join(SEQ[(rec[p+i//2]>>(4*(1-i%2)))&15] for i in range(lseq))
            yield pos,cig,sq
W=10
def evid(reads,pos,ref,alt):
    s0=pos-1; kind='snv' if len(ref)==len(alt) else ('ins' if len(alt)>len(ref) else 'del')
    L=abs(len(alt)-len(ref)); lo,hi=s0-W,s0+len(ref)+W
    n=0;A=0;R=0;O=0;C=0
    for rpos,cig,sq in reads:
        rp=rpos;qp=0;spans_lo=None;spans_hi=None;evs=[];base=None;clip=False
        for l,o in cig:
            if o in 'M=X':
                if rp<=s0<rp+l and kind=='snv': base=sq[qp+s0-rp:qp+s0-rp+len(alt)]
                if rp<=lo<rp+l: spans_lo=True
                if rp<=hi<rp+l: spans_hi=True
                rp+=l;qp+=l
            elif o=='I':
                if lo<=rp<=hi: evs.append(('I',l))
                qp+=l
            elif o=='D':
                if rp<=hi and rp+l>=lo: evs.append(('D',l))
                rp+=l
            elif o=='S':
                if lo<=rp<=hi: clip=True
                qp+=l
        if not (spans_lo or spans_hi or evs or clip): continue
        if kind=='snv':
            if base is None and not evs: continue
            n+=1
            if base==alt: A+=1
            elif evs: O+=1
            elif base==ref: R+=1
            else: O+=1
        else:
            match=[e for e in evs if e[0]==('I' if kind=='ins' else 'D') and abs(e[1]-L)<=max(2,0.2*L)]
            if match: n+=1;A+=1
            elif clip: n+=1;C+=1
            elif evs: n+=1;O+=1
            elif spans_lo and spans_hi: n+=1;R+=1
    return kind,n,A,R,O,C
if __name__=='__main__':
    bam,bai,sites,out=sys.argv[1:5]
    r=Reader(bam); names=header(r); IDX=load_bai(bai); T={n:i for i,n in enumerate(names)}
    with open(out,'w') as fo:
        for line in open(sites):
            c,p,ref,alt=line.rstrip('\n').split('\t'); p=int(p)
            rd=list(fetch(r,IDX,T[c],p-170,p+len(ref)+170))
            k,n,A,R,O,C=evid(rd,p,ref,alt)
            fo.write(f'{c}\t{p}\t{ref}\t{alt}\t{k}\t{n}\t{A}\t{R}\t{O}\t{C}\n')
