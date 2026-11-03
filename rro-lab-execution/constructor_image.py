"""Complete block encoding of W_U, with literal copier templates and bindings.

The source binary table is fully materialized. The expanded non-erasing table
and Wang word are represented by a finite closed block image. No evaluator
entry is synthesized from runtime input. Code fetch is an explicit interface
for paid-runtime B; this file's Python arithmetic is an image reader/reference,
not a backed LAB controller. Expanded words are not asserted to be stored.
"""
from __future__ import annotations
import argparse,bisect,gzip,hashlib,json,mmap,struct,tempfile
from pathlib import Path
from itertools import product
from copier import Copier

HERE=Path(__file__).resolve().parent


def templates():
    out={}
    for w0,w1,d0,d1 in product(range(2),range(2),('L','R'),('L','R')):
        kind=8*w0+4*w1+2*int(d0=='R')+int(d1=='R')
        c=Copier((((w0,d0,1),(w1,d1,2)),((0,'R',3),(1,'R',3)),((0,'R',3),(1,'R',3))))
        boundaries={v:k for k,v in c.boundaries.items()};size=boundaries[1];ext={boundaries[1]:-1,boundaries[2]:-2}
        rows=[]
        for row in c.rows[:size]:
            rr=[]
            for w,d,n in row:
                dest=ext[n] if n in ext else n
                assert dest<0 or dest<size,(kind,n,size)
                rr.append((w,d,dest))
            rows.append(rr)
        out[kind]=rows
    assert set(out)==set(range(16))
    return out

class BinaryTable:
    def __init__(self,path):
        self.f=open(path,'rb');self.data=mmap.mmap(self.f.fileno(),0,access=mmap.ACCESS_READ)
        assert self.data[:8]==b'RROTABLE';alpha,self.n=struct.unpack_from('<II',self.data,8)
        assert alpha==2 and len(self.data)==16+12*self.n
    def row(self,q):
        assert 0<=q<self.n
        return [struct.unpack_from('<BbI',self.data,16+12*q+6*b) for b in (0,1)]
    def kind(self,q):
        x,y=self.row(q)
        for w,d,n in (x,y):assert w<2 and d in (-1,1) and n<=self.n
        return 8*x[0]+4*y[0]+2*int(x[1]==1)+int(y[1]==1)

def build_image(binary_path,directory):
    directory=Path(directory);directory.mkdir(exist_ok=True)
    source=BinaryTable(binary_path);ts=templates();counts=[0]*16;nstates=0
    bases_path=directory/'copier-bases.bin'
    with bases_path.open('wb') as f:
        for q in range(source.n):
            f.write(struct.pack('<Q',nstates));k=source.kind(q);counts[k]+=1;nstates+=len(ts[k])
        f.write(struct.pack('<Q',nstates))
    with bases_path.open('rb') as bf:
        bases=mmap.mmap(bf.fileno(),0,access=mmap.ACCESS_READ)
        base=lambda q:struct.unpack_from('<Q',bases,8*q)[0]
        target=directory/'W_U.block.gz';digest=hashlib.sha256();logical_bytes=0
        with target.open('wb') as file:
            with gzip.GzipFile(filename='',mode='wb',fileobj=file,mtime=0,compresslevel=6) as f:
                def put(b):
                    nonlocal logical_bytes
                    f.write(b);digest.update(b);logical_bytes+=len(b)
                put(b'RROWUC1\0'+struct.pack('<IQ',source.n,nstates))
                for k in range(16):
                    put(struct.pack('<I',len(ts[k])))
                    for row in ts[k]:
                        for w,d,n in row:put(struct.pack('<Bbq',w,-1 if d=='L' else 1,n))
                put(bases)
                for q in range(source.n):
                    x,y=source.row(q);k=source.kind(q)
                    # Actual resolved destination addresses, including sole halt.
                    put(struct.pack('<BQQ',k,base(x[2]),base(y[2])))
        result={'schema':1,'source_binary_states':source.n,'template_count':16,'template_states':{str(k):len(v) for k,v in ts.items()},'template_use_counts':counts,'non_erasing_states':nstates,'non_erasing_rules':2*nstates,'W_U_instructions':13*nstates+1,'W_U_halt_PC':13*nstates+1,'image_bytes':target.stat().st_size,'uncompressed_image_bytes':logical_bytes,'image_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'uncompressed_image_sha256':digest.hexdigest(),'bases_bytes':len(bases),'bases_sha256':hashlib.sha256(bases).hexdigest(),'flat_W_U_materialized':False,'code_encoding':'Complete literal-template / resolved-binding block image, not a flat instruction array. Fetch/decompression and backing remain B obligations.'}
    (directory/'image-results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');return result

class Image:
    """Read a decompressed block image; addresses and words are public code."""
    def __init__(self,path):
        self.f=open(path,'rb');self.data=mmap.mmap(self.f.fileno(),0,access=mmap.ACCESS_READ)
        assert self.data[:8]==b'RROWUC1\0';self.n,self.nm=struct.unpack_from('<IQ',self.data,8);pos=20;self.ts=[]
        for k in range(16):
            n=struct.unpack_from('<I',self.data,pos)[0];pos+=4
            self.ts.append((pos,n));pos+=20*n
        self.base_start=pos;self.binding_start=pos+8*(self.n+1)
        assert len(self.data)==self.binding_start+17*self.n
    def base(self,q):return struct.unpack_from('<Q',self.data,self.base_start+8*q)[0]
    def transition(self,state,bit):
        assert state<self.nm and bit in (0,1)
        lo,hi=0,self.n
        while lo+1<hi:
            mid=(lo+hi)//2
            if self.base(mid)<=state:lo=mid
            else:hi=mid
        k,d0,d1=struct.unpack_from('<BQQ',self.data,self.binding_start+17*lo)
        at,n=self.ts[k];offset=state-self.base(lo);assert offset<n
        w,d,dest=struct.unpack_from('<Bbq',self.data,at+20*offset+10*bit)
        dest=(d0 if dest==-1 else d1) if dest<0 else self.base(lo)+dest
        return w,'L' if d==-1 else 'R',dest
    def fetch(self,pc):
        if pc==13*self.nm:return 'M'
        if pc==13*self.nm+1:return None
        assert 0<=pc<13*self.nm
        q,j=divmod(pc,13)
        if j==0:return 'R'
        if j==1:return ('J',13*q+8)
        read=0 if j<8 else 1;w,d,dest=self.transition(q,read)
        if read==0:
            block=(('R','M','M','M','M') if w==0 else ('M','R','M','M','M')) if d=='R' else (('L','L','L','M','M') if w==0 else ('M','L','L','L','M'))
            block=block+(('J',13*dest),);return block[j-2]
        block=(('R','M','M','M') if d=='R' else ('L','L','L','M'))+(('J',13*dest),);return block[j-8]


def image_checks():
    # Every finite template row, both read cases, and all 13 Wang entries,
    # compared with the unchanged generators on independently chosen programs.
    from nwmg import compile_non_erasing
    checked=states=0
    with tempfile.TemporaryDirectory(prefix='rro-image-') as td:
        td=Path(td)
        for k in range(16):
            w0,w1=k>>3,(k>>2)&1;d0='R' if k&2 else 'L';d1='R' if k&1 else 'L'
            source=(((w0,d0,0),(w1,d1,1)),)
            raw=b'RROTABLE'+struct.pack('<II',2,1)+b''.join(struct.pack('<BbI',w,-1 if d=='L' else 1,n) for row in source for w,d,n in row)
            b=td/'source';b.write_bytes(raw);result=build_image(b,td)
            ip=td/'image'
            with gzip.open(td/'W_U.block.gz','rb') as f:ip.write_bytes(f.read())
            im=Image(ip);full=Copier(source);wp=compile_non_erasing(full.rows)
            assert im.nm==len(full.rows)
            for q,row in enumerate(full.rows):
                for bit in (0,1):assert im.transition(q,bit)==row[bit]
                for j in range(13):assert im.fetch(13*q+j)==wp[13*q+j]
            assert im.fetch(len(wp)-1)=='M' and im.fetch(len(wp)) is None
            checked+=1;states+=im.nm
    return {'images':checked,'non_erasing_states_compared':states,'transition_cases':2*states,'W_U_instructions_compared':13*states+16}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--checks-only',action='store_true');args=ap.parse_args()
    print(json.dumps(image_checks()),flush=True)
    if not args.checks_only:print(json.dumps(build_image(args.binary,args.output)),flush=True)
