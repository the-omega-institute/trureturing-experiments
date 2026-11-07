"""FS1 bounded hosted macro probe, not an installed native controller.

Expected invariant: original paid LAB edits, full ancestor routes, distinct
retained identities, and stack/B boundary semantics agree at every tested cut.
Falsifiers: a changed old stack at a busy growth cut, an illegal edit, merged
payloads, omitted ancestor append, source replacement or lost receipt.
No evaluator result is computed; only bounded prefixes use its fixed graph.
"""
from __future__ import annotations
import hashlib,json,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from lab import Arena,Lab
from constructor_checks import graph_bytes
from evaluator_program import build
from framing import parse_arithmetic,arithmetic_input,REPLIES
from session import prepare_source,OriginalPort,Session

HERE=Path(__file__).resolve().parent
OLD=HERE.parent
assert sys.version_info >= (3,12) and __debug__

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Frontier:
    """FIFO branch directory. All queue records and routes are retained."""
    def __init__(self,a):
        self.a=a;self.lab=Lab(a,capture_cuts=False)
        self.header=a.alloc('FS1_queue_header')
        for k in ('head','tail'):a.write(self.header,k,None)
        self.retained=[];self.data=[];self.depths=[]
        assert self.lab.grow('P',self.lab.root)
        self.enqueue(self.lab.root)
    def enqueue(self,p):
        a=self.a;q=a.alloc('FS1_queue_record');self.retained.append(q)
        a.write(q,'value',p);a.write(q,'next',None)
        tail=a.read(self.header,'tail')
        if a.test('null',tail is None):a.write(self.header,'head',q)
        else:a.write(tail,'next',q)
        a.write(self.header,'tail',q)
    def bit(self,b):
        a=self.a;lab=self.lab
        q=a.read(self.header,'head');assert q is not None
        p=a.read(q,'value');n=a.read(q,'next')
        a.write(self.header,'head',n)
        if a.test('null',n is None):a.write(self.header,'tail',None)
        self.depths.append(len(lab.address(p)))
        assert lab.grow('S',p)
        left=a.read(p,'left');right=a.read(p,'right')
        assert lab.grow('P',right) and lab.grow('S',right)
        second=a.read(right,'left');data=a.read(right,'right')
        if b:assert lab.grow('P',data)
        self.enqueue(left);self.enqueue(second);self.data.append(data)
        return data
    def value(self,p):
        tag=self.lab.tag(p);assert tag in ('alpha','beta')
        return int(tag=='beta')
    def check(self):
        self.lab.validate_routes()
        t=len(self.data);ones=sum(self.value(p) for p in self.data)
        d=sum(3*x+2 for x in self.depths)
        d+=sum(depth+2 for depth,p in zip(self.depths,self.data) if self.value(p))
        c=self.lab.persistent()
        assert (c['V'],c['H'],c['D'])==(2+5*t+ones,2*t,d)
        assert len({id(x) for x in self.data})==t
        assert len(self.retained)==2*t+1
        assert all(depth == (i+1).bit_length()-1+(i+1).bit_count()-1 for i,depth in enumerate(self.depths))
        return dict(t=t,ones=ones,depth_max=max(self.depths,default=0),**c)

class StackProbe:
    def __init__(self,p,bits):
        self.a=Arena(registers=32);self.f=Frontier(self.a);self.p=p
        self.regs=[];self.frames=[];self.pc=p['start'];self.commits=0;self.busy=None
        self.cut_count=0;self.cut_digest=hashlib.sha256()
        for _ in p['registers']:
            r=self.a.alloc('FS1_register');self.a.write(r,'head',None);self.regs.append(r)
        self.loading=True
        for b in reversed(bits):self.push(0,int(b))
        self.loading=False
        def observe(phase):
            self.cut_count+=1
            assert self.busy is not None
            assert self.words()==self.busy['old_words']
            self.cut_digest.update(phase.encode());self.cut_digest.update(str(self.a.census()).encode())
        self.f.lab.observer=observe
    def words(self):
        out=[]
        for reg in self.regs:
            word=[];p=reg.fields['head'];seen=set()
            while p is not None:
                assert p not in seen;seen.add(p)
                word.append(int(Lab.peek_tag(p.fields['value'])=='beta'));p=p.fields['next']
            out.append(word)
        return out
    def push(self,r,b):
        a=self.a;old=a.read(self.regs[r],'head');data=self.f.bit(b)
        frame=a.alloc('FS1_stack_frame');self.frames.append(frame)
        a.write(frame,'value',data);a.write(frame,'next',old)
        a.write(self.regs[r],'head',frame)
    def step(self):
        x=self.p['code'][self.pc];kind=x[0]
        if kind=='halt':return False
        self.busy=dict(old_pc=self.pc,old_words=self.words(),op=kind)
        if kind=='goto':n=x[1]
        elif kind=='push':self.push(x[1],x[2]);n=x[3]
        else:
            r=self.regs[x[1]];p=self.a.read(r,'head')
            if self.a.test('null',p is None):v=None
            else:v=self.f.value(self.a.read(p,'value'))
            n=x[2 if v is None else 3+v]
            if kind=='pop' and p is not None:self.a.write(r,'head',self.a.read(p,'next'))
        self.pc=n;self.commits+=1;self.busy=None;return True

def refstep(p,pc,regs):
    x=p['code'][pc];kind=x[0]
    if kind=='halt':return None
    if kind=='goto':return x[1]
    r=regs[x[1]]
    if kind=='push':r.insert(0,x[2]);return x[3]
    v=r[0] if r else None
    if kind=='pop' and r:r.pop(0)
    return x[2 if v is None else 3+v]

def stack_checks():
    fixed=build();raw=graph_bytes(fixed)
    assert raw==(OLD/'constructor/evaluator.stack.bin').read_bytes()
    tests=[('fixed-evaluator-prefix',fixed,(OLD/'constructor/positive-full.input.bits').read_text(),96)]
    # Empty, zero and one branches, pop and peek, repeated equal-valued slots.
    p={'registers':{'input':0},'start':0,'code':[
       ('peek',0,1,4,6),('push',0,0,2),('push',0,1,3),('pop',0,8,8,4),
       ('pop',0,8,5,8),('goto',0),('push',0,1,7),('pop',0,8,8,0),('halt',True)]}
    for bits in ('','0','1','00'):tests.append(('branch-'+bits,p,bits,20))
    rows=[]
    for name,p,bits,fuel in tests:
        s=StackProbe(p,bits);regs=[[] for _ in p['registers']];regs[0]=list(map(int,bits));pc=p['start']
        for _ in range(fuel):
            nxt=refstep(p,pc,regs)
            if nxt is None:assert not s.step();break
            assert s.step();pc=nxt;assert s.pc==pc and s.words()==regs
        counts=s.f.check()
        rows.append(dict(name=name,commits=s.commits,pc=pc,halt=p['code'][pc][0]=='halt',
                         busy_growth_cuts=s.cut_count,cut_sha256=s.cut_digest.hexdigest(),
                         work=counts,census=s.a.census(),charges=dict(s.a.charges)))
    return dict(rows=rows,fixed_graph_sha256=hashlib.sha256(raw).hexdigest(),
                installed_PC=False,grant_refusal_implementation=False,
                scope='hosted macro correspondence; fixed evaluator prefix only; no returned arithmetic result')

def b_checks():
    rows=[]
    tests=[([], '',0),(['R','M'],'00',2),(['M','M'], '0',3),(['R','L','L','R','M'], '00',8),
           ([('J',2),'M','L'],'0',8),([('J',2),'M','L'],'1',8),
           (['R',('J',0)],'11',24),(['M',('J',0)],'',24)]
    for code,bits,fuel in tests:
        a=Arena();f=Frontier(a);records=[]
        def cell(b):
            p=a.alloc('FS1_tape_record');records.append(p)
            for k,v in (('left',None),('right',None),('value',f.bit(b))):a.write(p,k,v)
            return p
        old=None;origin=None
        for b in bits or '0':
            p=cell(int(b))
            if old is not None:a.write(old,'right',p);a.write(p,'left',old)
            else:origin=p
            old=p
        head=origin;pc=h=0;tape={i:int(b) for i,b in enumerate(bits) if int(b)};commits=0
        first_events=tuple(f.lab.events);geometry_witness=False
        for _ in range(fuel):
            if pc==len(code):break
            op=code[pc];n=pc+1
            if op in ('L','R'):
                k,back,delta=('left','right',-1) if op=='L' else ('right','left',1)
                q=a.read(head,k)
                if a.test('null',q is None):q=cell(0);a.write(q,back,head);a.write(head,k,q)
                head=q;h+=delta
            elif op=='M':
                v=a.read(head,'value')
                if f.value(v)==0:assert f.lab.grow('P',v)
                tape[h]=1
            else:n=op[1] if f.value(a.read(head,'value')) else n
            pc=n;commits+=1
            assert f.value(a.read(head,'value'))==tape.get(h,0)
            mapped={0:origin};q=origin
            for delta,key,back in ((1,'right','left'),(-1,'left','right')):
                q=origin;i=0
                while q.fields[key] is not None:
                    nxt=q.fields[key];assert nxt.fields[back] is q
                    i+=delta;assert i not in mapped;mapped[i]=nxt;q=nxt
            assert len(mapped)==len(records) and mapped[h] is head
            assert {i:1 for i,r in mapped.items() if f.value(r.fields['value'])}=={i:v for i,v in tape.items() if v}
            if code==['R','M'] and commits==1:
                assert tuple(f.lab.events)==first_events and head is not origin and pc==1
                geometry_witness=True
        # Exact permanent tree formula includes later M marks in their slot.
        work=f.check();before=tuple(f.lab.events);obj=f.data[0]
        if f.value(obj):
            assert not f.lab.grow('P',obj);assert tuple(f.lab.events)==before
        rows.append(dict(code=code,input=bits,commits=commits,halt=pc==len(code),head=h,
                         work=work,charges=dict(a.charges),same_bare_tree_witness=geometry_witness))
    return rows

def source_checks():
    rows=[]
    for name in ('positive-full','positive-partial','null-nine-distinct-left'):
        raw=(OLD/'constructor'/(name+'.input.bits')).read_text();n,theta,cache,cyl=parse_arithmetic(raw)
        bits=[dict(row).get(h,0) for row in cache for h in range(15)]
        source=prepare_source(bits);root=source.root;events=tuple(source.events)
        port=OriginalPort(source);s=Session(n,theta,port)
        calls=['','LL','LRR','LLR','LL','LR','LLL','LRL']
        calls+=['R'*b+'LL' for b in range(15*n) if (b%15,bits[b]) in cache[b//15]]
        calls+=['LL','R'*15+'LL','R'*45+'LL','L'*130,'R'*(15*n+10)]
        for w in calls:s.read(w)
        assert s.cache_tuple()==cache
        assert arithmetic_input(n,theta,s.cache_tuple(),cyl)==raw
        assert set(x[1] for x in port.calls)==set(REPLIES)
        assert len(s.records)==len(calls)
        assert s.records[1][0]==s.records[4][0]=='LL'
        assert s.records[1][2].header is not s.records[4][2].header
        assert source.root is root and tuple(source.events)==events
        before=len(port.calls);s.read('LL')
        assert len(port.calls)==before+1 and s.port is port and port.root is root
        rows.append(dict(case=name,N=n,accepted_frame_sha256=hashlib.sha256(raw.encode()).hexdigest(),
                         calls=len(port.calls),operand_length_sum=sum(len(x[0]) for x in port.calls),
                         replies=sorted(set(x[1] for x in port.calls)),occurrence_receipts=len(s.records),
                         source_census=source.a.census(),observer_census=s.a.census(),
                         scope='unchanged existing hosted Session only; no evaluator/output/continuation integration into FS1'))
    return rows

def main():
    # Measurements/results only; no execution transcript or per-cut snapshots.
    result=dict(schema=1,candidate='FS1-fifo-directory-persistent-stacks',
                candidate_complete=False,old_C_complete=False,
                stack=stack_checks(),B=b_checks(),source=source_checks(),
                source_hashes={p.name:digest(p) for p in [Path(__file__),OLD/'lab.py',OLD/'session.py',
                    OLD/'framing.py',OLD/'constructor_ir.py',OLD/'evaluator_program.py']})
    out=HERE/'frontier-stack-results.json';out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(result=out.name,sha256=digest(out),candidate_complete=False,exit=0)))
if __name__=='__main__':main()
