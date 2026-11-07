"""Executable content port of the finite StackSimulation controller.

Apache-2.0 derivative: RBarish-UTokyo/NPCompleteness-Lean4 revision
3e603e795b44fedcf7d6dd65b7d1e30547ead925, Complexity/{StackCompiler,
StackPop,StackEncoding,StackSimulation}.lean. See compiler-provenance.md
and RBarishBinaryMachine-LICENSE.txt. Symbolic labels replace Fin encodings;
reachable labels are enumerated injectively, without changing a local rule.
The explicit equal-length binary adapter below is repo-derived ordinary code.
"""
from __future__ import annotations
from collections import deque,Counter
import struct,json,hashlib,gzip
from pathlib import Path

# Symbols: blank=0, bit false=1, bit true=2, separator=3.
# Moves: -1 left, 0 stay, +1 right. None is a terminal instruction.

def four_rule(code,q,a):
    k=q[0]
    if k=='finish':return None
    if k=='main':
        x=code[q[1]];op=x[0]
        if op=='halt':return a,1,('finish',x[1])
        if op=='goto':return a,0,('main',x[1])
        if op=='push':return a,0,('push',x[3],x[2],'seek',x[1]+1)
        return a,0,('seek',q[1],x[1]+1)
    if k=='seek':
        _,source,c=q
        if c==0:return a,0,('inspect',source)
        if a==0:return None
        return a,1,('seek',source,c-int(a==3))
    if k=='inspect':
        x=code[q[1]];op=x[0]
        if op not in ('pop','peek') or a==0:return None
        nxt=x[2 if a==3 else a+2]
        if op=='pop' and a!=3:return a,0,('delete',nxt,0)
        return a,0,('rewind',nxt,0)
    if k=='rewind':
        _,n,phase=q
        if phase==1:return a,0,('main',n)
        return (0,1,('rewind',n,1)) if a==0 else (a,-1,q)
    if k=='delete':
        _,n,phase=q
        if phase==0:return a,1,('delete',n,1)
        if phase==1:return (0,-1,('delete',n,3)) if a==0 else (a,-1,('delete',n,5+a))
        if phase==2:return a,1,('delete',n,1)
        if phase==3:return 0,0,('delete',n,4)
        if phase==4:return a,0,('delete',n,9) # success -> rewindFromBlank
        if 5<=phase<=8:return phase-5,1,('delete',n,2)
        if phase==9:return a,-1,('delete',n,10)
        if phase==10:return a,0,('delete',n,11)
        if phase==11:return (0,1,('delete',n,12)) if a==0 else (a,-1,q)
        if phase==12:return a,0,('main',n)
        raise AssertionError(q)
    if k=='push':
        _,n,b,phase,c=q
        if phase=='seek':
            if c==0:return a,0,('push',n,b,'insert',1+b)
            if a==0:return None
            return a,1,('push',n,b,'seek',c-int(a==3))
        if phase=='insert':
            if a==0:return c,0,('push',n,b,'insert-halt',0)
            return c,1,('push',n,b,'insert',a)
        if phase=='insert-halt':return a,0,('push',n,b,'rewind',0)
        if phase=='rewind':
            if c==1:return a,0,('main',n)
            return (0,1,('push',n,b,'rewind',1)) if a==0 else (a,-1,q)
        raise AssertionError(q)
    raise AssertionError(q)


def reachable_four(program,progress=None):
    code=program['code'];labels=[('main',program['start'])];index={labels[0]:0};rows=[]
    for q in labels:
        row=[]
        for symbol in range(4):
            rule=four_rule(code,q,symbol)
            if rule is None:row.append(None);continue
            w,d,n=rule
            if n not in index:index[n]=len(labels);labels.append(n)
            row.append((w,d,index[n]))
        rows.append(tuple(row))
        if progress and len(rows)%100000==0:progress(len(rows),len(labels))
    assert len(rows)==len(labels)
    return labels,rows


def binary_row(four,q,phase):
    """Eleven allocated states per source state; all rows are materialized.

    0 is pair-aligned; 1/2 retain the first bit; 3..6 write the first
    result bit; 7..10 finish the selected two-square movement. Unreachable
    phases attached to a terminal symbol still have closed preserving rules.
    """
    halt=11*len(four)
    if phase==0:return ((0,'R',11*q+1),(1,'R',11*q+2))
    if phase in (1,2):
        hi=phase-1;row=[]
        for b in (0,1):
            sym=2*hi+b;r=four[q][sym]
            row.append((b,'L',halt) if r is None else (r[0]%2,'L',11*q+3+sym))
        return tuple(row)
    sym=(phase-3)%4;r=four[q][sym]
    if r is None:return ((0,'L',halt),(1,'L',halt))
    w,d,n=r
    if phase<7:return ((w//2,'L' if d==-1 else 'R',11*q+7+sym),)*2
    return tuple((b,'R' if d==1 else 'L',11*n) for b in (0,1))

def aligned_binary(four):
    rows=[binary_row(four,q,phase) for q in range(len(four)) for phase in range(11)]
    return None,rows

class MappedFour:
    def __init__(self,path):
        import mmap
        self.file=open(path,'rb');self.data=mmap.mmap(self.file.fileno(),0,access=mmap.ACCESS_READ)
        assert self.data[:8]==b'RROTABLE'
        alphabet,self.n=struct.unpack_from('<II',self.data,8);assert alphabet==4 and len(self.data)==16+24*self.n
    def __len__(self):return self.n
    def __getitem__(self,q):
        assert 0<=q<self.n
        out=[]
        for a in range(4):
            w,d,n=struct.unpack_from('<BbI',self.data,16+24*q+6*a)
            assert n==0xffffffff or (w<4 and d in (-1,0,1) and n<self.n)
            out.append(None if n==0xffffffff else (w,d,n))
        return out

def stream_adapter(source,destination):
    four=MappedFour(source);n=len(four);digest=hashlib.sha256();size=0
    with open(destination,'wb') as f:
        def put(b):
            nonlocal size
            f.write(b);digest.update(b);size+=len(b)
        put(b'RROTABLE'+struct.pack('<II',2,11*n))
        # One source row at a time, with exhaustive checks of all 22 binary
        # transitions. No expansion dictionary or input is used.
        class OneRow:
            def __init__(self,row,n):self.row=row;self.n=n
            def __len__(self):return self.n
            def __getitem__(self,q):return self.row
        for q in range(n):
            row=four[q];view=OneRow(row,n)
            for phase in range(11):
                for w,d,dest in binary_row(view,q,phase):
                    assert w in (0,1) and d in ('L','R') and 0<=dest<=11*n
                    put(struct.pack('<BbI',w,-1 if d=='L' else 1,dest))
    return {'binary_states':11*n,'binary_rules':22*n,'binary_bytes':size,'binary_sha256':digest.hexdigest(),'four_states':n,'layout':'11 states per four-symbol state; row-major little-endian RROTABLE'}


def initial_four(bits,nregisters):
    return [3]+[1+int(b) for b in bits]+[3]*nregisters

def initial_binary(bits,nregisters):
    return ''.join(format(a,'02b') for a in initial_four(bits,nregisters))

def run_four(rows,initial,fuel):
    tape={i:a for i,a in enumerate(initial) if a};h=q=0;steps=0
    while steps<fuel:
        r=rows[q][tape.get(h,0)];steps+=1
        if r is None:return tape,h,steps
        w,d,q=r;tape[h]=w;h+=d
    raise RuntimeError(('four fuel',steps,q,h))

def run_binary(rows,initial,fuel):
    tape={i:int(a) for i,a in enumerate(initial) if int(a)};h=q=0;steps=0
    while q<len(rows) and steps<fuel:
        w,d,q=rows[q][tape.get(h,0)];tape[h]=w;h+=1 if d=='R' else -1;steps+=1
    if q!=len(rows):raise RuntimeError(('binary fuel',steps,q,h))
    return tape,h,steps


def decode_four(tape,head):
    out=[];p=head
    while tape.get(p,0) in (1,2):out.append(str(tape[p]-1));p+=1
    return ''.join(out)

def decode_binary(tape,head):
    assert head%2==0
    def read(p):return 2*tape.get(2*p,0)+tape.get(2*p+1,0)
    return decode_four({p:read(p) for p in range(head//2, max(tape,default=0)//2+2)},head//2)


def table_bytes(rows,alphabet):
    out=bytearray(b'RROTABLE'+struct.pack('<II',alphabet,len(rows)))
    for row in rows:
        for r in row:
            if r is None:out.extend(struct.pack('<BbI',0,0,0xffffffff))
            else:
                w,d,q=r;out.extend(struct.pack('<BbI',w,(-1 if d=='L' else 1) if isinstance(d,str) else d,q))
    return bytes(out)


def generic_checks():
    from constructor_ir import Assembler,primitives,run_graph
    # Every primitive stack operation and its empty/zero/one branch on all
    # represented words of length <=3, including unequal stack lengths.
    from itertools import product
    count=0;fsteps=bsteps=0
    for op in ['push0','push1','pop','peek']:
        for register in range(3):
            code=[('push',register,int(op[-1]),1),('halt',True)] if op.startswith('push') else [(op,register,1,2,3),('halt',True),('halt',True),('halt',True)]
            p={'start':0,'code':code,'registers':{str(i):i for i in range(3)}}
            labels,four=reachable_four(p);bl,binary=aligned_binary(four)
            for length in range(4):
                for word in product((0,1),repeat=length):
                    regs=[list(word),[1,0],[]]
                    initial=[]
                    for r in regs:initial.extend([3]+[1+b for b in r])
                    initial.append(3)
                    expected=[x.copy() for x in regs]
                    if op.startswith('push'):expected[register].insert(0,int(op[-1]))
                    elif op=='pop' and expected[register]:expected[register].pop(0)
                    ft,fh,fs=run_four(four,initial,10000)
                    bt,bh,bs=run_binary(binary,''.join(format(x,'02b') for x in initial),100000)
                    expected_tape=[]
                    for r in expected:expected_tape.extend([3]+[1+b for b in r])
                    expected_tape.append(3)
                    assert fh==1 and bh==2,(op,register,word,fh,bh)
                    assert all(ft.get(i,0)==x for i,x in enumerate(expected_tape))
                    assert all(2*bt.get(2*i,0)+bt.get(2*i+1,0)==x for i,x in enumerate(expected_tape))
                    assert decode_four(ft,fh)==''.join(map(str,expected[0]))==decode_binary(bt,bh)
                    count+=1;fsteps+=fs;bsteps+=bs
    return {'instances':count,'four_steps_including_halt':fsteps,'binary_steps':bsteps}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);ap.add_argument('--checks-only',action='store_true');ap.add_argument('--from-four',type=Path);args=ap.parse_args()
    print(json.dumps(generic_checks()),flush=True)
    if args.from_four:
        result=stream_adapter(args.from_four,args.output/'evaluator.binary.bin')
        old=json.loads((args.output/'compiler-results.json').read_text());old.update(result);(args.output/'compiler-results.json').write_text(json.dumps(old,indent=2,sort_keys=True)+'\n');print(json.dumps(old),flush=True)
    elif not args.checks_only:
        from evaluator_program import build
        p=build();labels,four=reachable_four(p,lambda n,m:print('four closure',n,m,flush=True));print('four complete',len(four),flush=True)
        args.output.mkdir(exist_ok=True)
        raw=table_bytes(four,4);(args.output/'evaluator.four.bin').write_bytes(raw)
        summary={'four_states':len(four),'four_bytes':len(raw),'four_sha256':hashlib.sha256(raw).hexdigest(),'phase_counts':dict(Counter(q[0] for q in labels))}
        (args.output/'compiler-results.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n');print(json.dumps(summary),flush=True)
        summary.update(stream_adapter(args.output/'evaluator.four.bin',args.output/'evaluator.binary.bin'))
        (args.output/'compiler-results.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n');print(json.dumps(summary),flush=True)
