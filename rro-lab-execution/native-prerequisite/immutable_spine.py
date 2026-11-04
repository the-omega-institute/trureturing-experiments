"""IS2: changed immutable-spine route, bounded author discrimination only.

Invariant: copied heads alias immutable spines, published words change only at
literal register writes, token leaves are never marked after P/S initialization.
Falsifiers: an old published tail changes, an unrepresented selector dispatches,
a guard/refusal consumes a pool cell, a source input bit is replaced, or an
output bit differs from the actual computed register-zero spine.
No full UJB consumer run or host resource sufficiency is claimed.
"""
from __future__ import annotations
from collections import Counter
from contextlib import contextmanager
import hashlib
import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE))
import constructor_ir as ci
import evaluator_program as ep
from field_graph import Graph
from native_install import Block, Source, Installer, BitKernel
from framing import reply, parse_arithmetic, arithmetic_input, REPLIES
from session import Session, OriginalPort, prepare_source

assert sys.version_info >= (3, 12) and __debug__
IDENTITY = 'IS2-immutable-spines-direct-native/v1'
GRAPH_LABELS = {}  # Source decoder for independent checks; never runtime dispatch.

class SharingGraph(ci.Graph):
    def clear(self, r): self.emit('clear', self.reg(r))
    def copy(self, s, d):
        if s != d: self.emit('copy', self.reg(d), self.reg(s))
    def finish(self):
        out = []
        for i, x in enumerate(self.code):
            k = x[0]
            if k == 'halt': out.append(x)
            elif k == 'goto': out.append((k, self.labels[x[1]]))
            elif k in ('push', 'copy', 'clear'): out.append((*x, i+1))
            else: out.append((k, x[1], *(self.labels[a] for a in x[2:])))
        p = dict(registers=self.registers, labels=self.labels, start=0, code=out)
        assert all(0 <= q < len(out) for x in out for q in successors(x))
        return p

@contextmanager
def shared_assembler():
    original = ci.Graph
    ci.Graph = SharingGraph
    try: yield
    finally: ci.Graph = original

def successors(x):
    k = x[0]
    if k == 'halt': return ()
    if k == 'goto': return (x[1],)
    if k in ('clear', 'copy', 'push'): return (x[-1],)
    return x[2:]

def prune(p):
    seen, todo = set(), [p['start']]
    while todo:
        q = todo.pop()
        if q not in seen: seen.add(q); todo.extend(successors(p['code'][q]))
    order = sorted(seen); ids = {q:i for i,q in enumerate(order)}; out=[]
    for q in order:
        x=p['code'][q];k=x[0]
        if k=='halt': out.append(x)
        elif k=='goto': out.append((k,ids[x[1]]))
        elif k in ('clear','copy','push'): out.append((*x[:-1],ids[x[-1]]))
        else: out.append((k,x[1],*(ids[z] for z in x[2:])))
    return dict(registers=p['registers'],start=ids[p['start']],code=out)

def canonical(p):
    return json.dumps(p,sort_keys=True,separators=(',',':')).encode()

def field_graph(p):
    """Every selector is literal wiring. No word ALU, pages or B compiler."""
    g=Graph();g.ref_names.update('stack.'+str(i) for i in p['registers'].values());g.label('start');g.j('seed')
    def assign(d,s):g.cp(d,s)
    def allocate(d):g.allocate(d)
    def reg(i):return 'stack.'+str(i)
    def emitbit(b):g.put('print.bit',b);g.emit('emit','print.bit')
    def branch_token(token,z,o):
        check=g.fresh('check.zero');g.equal(token,'token.one',o,check)
        g.label(check);g.equal(token,'token.zero',z,'fault')
    # Known actual alpha; all work changes are the literal P then S templates.
    g.label('seed');allocate('work.root')
    for i in range(4):g.set('work.root',i,'nil')
    g.get('guard.left','work.root',1);g.null('guard.left','seed.p.right','fault')
    g.label('seed.p.right');g.get('guard.right','work.root',2)
    g.null('guard.right','seed.p.reserve','fault')
    g.label('seed.p.reserve');allocate('token.one')
    g.set('token.one',0,'work.root');g.set('work.root',1,'token.one')
    # S guard: beta, with its unary marker still a physical alpha.
    g.get('guard.left','work.root',1);g.null('guard.left','fault')
    g.get('guard.right','work.root',2);g.null('guard.right','seed.s.reserve','fault')
    g.label('seed.s.reserve');allocate('token.marker');allocate('token.zero');allocate('route.header')
    g.set('token.marker',0,'token.one');g.set('token.zero',0,'work.root')
    g.set('token.one',1,'token.marker');g.set('work.root',2,'token.zero')
    g.set('work.root',3,'route.header')
    # No ancestors at either root event, hence no omitted route append.
    g.label('load');assign('load.cursor','received');assign('load.head','nil');assign('load.tail','nil')
    g.label('load.loop');g.null('load.cursor','load.end')
    g.read('load.bit','load.cursor');g.test('load.bit','load.one')
    assign('load.value','token.zero');g.j('load.append')
    g.label('load.one');assign('load.value','token.one')
    g.label('load.append');allocate('load.new');g.set('load.new',0,'nil');g.set('load.new',1,'load.value')
    g.null('load.tail','load.first');g.set('load.tail',0,'load.new');g.j('load.linked')
    g.label('load.first');assign('load.head','load.new')
    g.label('load.linked');assign('load.tail','load.new');g.get('load.cursor','load.cursor',1);g.j('load.loop')
    g.label('load.end');assign(reg(0),'load.head');g.j('q.'+str(p['start']))
    for i,x in enumerate(p['code']):
        g.label('q.'+str(i));k=x[0]
        if k in ('push','pop','copy','clear'):assign('busy.old',reg(x[1]))
        if k=='halt':g.put('result.ok',int(x[1]));g.j('print' if x[1] else 'reject')
        elif k=='goto':g.j('q.'+str(x[1]))
        elif k=='clear':assign(reg(x[1]),'nil');g.j('q.'+str(x[2]))
        elif k=='copy':assign(reg(x[1]),reg(x[2]));g.j('q.'+str(x[3]))
        elif k=='push':
            allocate('push.new');g.set('push.new',0,reg(x[1]));g.set('push.new',1,'token.one' if x[2] else 'token.zero')
            assign(reg(x[1]),'push.new');g.j('q.'+str(x[3]))
        else:
            assign('inspect.old',reg(x[1]));g.null('inspect.old','q.'+str(x[2]))
            g.get('inspect.value','inspect.old',1)
            if k=='pop':g.get(reg(x[1]),'inspect.old',0)
            branch_token('inspect.value','q.'+str(x[3]),'q.'+str(x[4]))
    # Print precisely reply(1, computed register zero), without a reference answer.
    g.label('print');emitbit(0);emitbit(1);assign('print.cursor',reg(0))
    g.label('print.length');g.null('print.cursor','print.delimiter');emitbit(1)
    g.get('print.cursor','print.cursor',0);g.j('print.length')
    g.label('print.delimiter');emitbit(0);assign('print.cursor',reg(0))
    g.label('print.payload');g.null('print.cursor','print.end')
    g.get('print.value','print.cursor',1);branch_token('print.value','print.zero','print.one')
    g.label('print.zero');emitbit(0);g.j('print.next')
    g.label('print.one');emitbit(1)
    g.label('print.next');g.get('print.cursor','print.cursor',0);g.j('print.payload')
    g.label('print.end');g.emit('end');g.emit('ack');g.emit('return');g.j('done')
    g.label('reject');emitbit(0);emitbit(0);emitbit(0);g.emit('end');g.emit('ack');g.emit('return');g.j('done')
    g.label('done');g.emit('wait');g.j('done')
    g.label('fault');g.emit('fault','invalid immutable spine/token')
    g.label('exhausted');g.emit('fault','no supplied cell')
    program=g.seal();GRAPH_LABELS[id(program)]=dict(g.labels);return program

# Literal native operations already in ML-P. The graph below is specialized to
# fixed field selectors; it does not execute field_graph.Engine or its word ALU.
OPS=('rr','rw','br','bw','address','rfield','wfield','bfield','vfield',
     'literal','nil','save','null','equal','test','goto','send','end','ack','return','fault','wait')

class Direct:
    def __init__(self,p,bits,capacity,kernel):
        assert capacity>=0
        self.p=p;self.bits_input=bits;self.output=[];self.events=Counter();self.steps=0
        self.source_blocks=[];self.native=[];self.boundary={};self.sources=[]
        def block(refs=4,width=8):
            b=Block([None]*refs,bytearray((width+7)//8),width);b.grant=1;self.source_blocks.append(b);return b
        source_cpu=block()
        self.source_regs={k:block(1,0) for k in p.refs}
        self.source_bits={k:block(0,1) for k in p.bits}
        pool=[block() for _ in range(capacity)]
        for a,b in zip(pool,pool[1:]):a.refs[0]=b
        self.source_regs['pool'].refs[0]=pool[0] if pool else None
        source_input=[block() for _ in bits]
        for b,ch in zip(source_input,bits):b.put(0,int(ch))
        for a,b in zip(source_input,source_input[1:]):a.refs[1]=b
        self.source_regs['received'].refs[0]=source_input[0] if source_input else None
        self.source_pool=pool;self.source_input=source_input
        def node(op,a=None,i=0,z=None,o=None):
            b=block(4,21);b.refs[:]=[a,None,z,z if o is None else o]
            b.bits=OPS.index(op)|(i<<5);self.native.append(b);return b
        def seq(actions,end):
            for op,a,i in reversed(actions):end=node(op,a,i,end)
            return end
        def rr(name):return ('rr',self.source_regs[name],0)
        def rw(name):return ('rw',self.source_regs[name],0)
        def br(name):return ('br',self.source_bits[name],0)
        def bw(name):return ('bw',self.source_bits[name],0)
        def addr(name):return [rr(name),('address',None,0)]
        entries=[node('goto') for _ in p.nodes]
        field_ids={n:i for i,n in enumerate(p.nodes)}
        for i,n in enumerate(p.nodes):
            op=n.op;a=n.args;z=entries[field_ids[n.zero]];o=entries[field_ids[n.one]]
            actions=[];end=z
            if op=='copy':actions=[rr(a[1]),rw(a[0])]
            elif op=='put':actions=[('literal',None,a[1]),bw(a[0])]
            elif op=='null':end=node('null',z=z,o=o);actions=[rr(a[0])]
            elif op=='equal':end=node('equal',z=z,o=o);actions=[rr(a[0]),('save',None,0),rr(a[1])]
            elif op=='test':end=node('test',z=z,o=o);actions=[br(a[0])]
            elif op=='get':actions=addr(a[1])+[('rfield',None,a[2]),rw(a[0])]
            elif op=='set':actions=addr(a[0])+[rr(a[2]),('wfield',None,a[1])]
            elif op=='read':actions=addr(a[1])+[('bfield',None,a[2]),bw(a[0])]
            elif op=='write':actions=addr(a[0])+[br(a[2]),('vfield',None,a[1])]
            elif op=='emit':actions=[br(a[0]),('send',None,0)]
            elif op in ('end','ack','return','wait','fault'):end=node(op,z=z,o=o)
            elif op!='goto':raise AssertionError(op)
            first=seq(actions,end);entries[i].refs[2:]=[first,first]
            self.boundary[entries[i]]=i
        self.source_entry=entries[0]
        self.source_boundaries=dict(self.boundary)
        self.sources=[Source(b) for b in self.source_blocks];mapping=dict(zip(self.source_blocks,self.sources))
        for s in self.sources:s.block.refs[:]=[mapping[v] if v is not None else None for v in s.block.refs]
        for a,b in zip(self.sources,self.sources[1:]):a.next=b
        self.inst=Installer(self.sources,kernel,mapping[self.source_entry]);self.inst.run();self.inst.verify()
        target=lambda b:mapping[b].forward
        self.regs={k:target(v) for k,v in self.source_regs.items()};self.bitregs={k:target(v) for k,v in self.source_bits.items()}
        self.boundary={target(b):i for b,i in self.source_boundaries.items()}
        self.pool=[target(b) for b in pool];self.received=[target(b) for b in source_input]
        self.cpu=target(source_cpu);self.pc=self.inst.execution_pc;self.events['handoff_pc_write']=1;self.width=self.inst.width
        self.labels={v:k for k,v in GRAPH_LABELS[id(p)].items()};self.frozen={b:tuple(b.refs) for b in self.received};self.cut_checks=0
        self.allocated=[];self.snapshots={};self.last_pool=self.regs['pool'].refs[0]
    @property
    def rr(self):return self.cpu.refs[0]
    @rr.setter
    def rr(self,v):self.cpu.refs[0]=v
    @property
    def saved(self):return self.cpu.refs[1]
    @saved.setter
    def saved(self,v):self.cpu.refs[1]=v
    @property
    def address(self):return self.cpu.refs[2]
    @address.setter
    def address(self,v):self.cpu.refs[2]=v
    @property
    def pc(self):return self.cpu.refs[3]
    @pc.setter
    def pc(self,v):self.cpu.refs[3]=v
    @property
    def bit(self):return self.cpu.bit(0)
    @bit.setter
    def bit(self,v):self.cpu.put(0,v)
    @property
    def done(self):return bool(self.cpu.bit(1))
    @done.setter
    def done(self,v):self.cpu.put(1,int(v))
    @property
    def refused(self):return bool(self.cpu.bit(2))
    @refused.setter
    def refused(self,v):self.cpu.put(2,int(v))
    @property
    def ended(self):return bool(self.cpu.bit(3))
    @ended.setter
    def ended(self,v):self.cpu.put(3,int(v))
    @property
    def acked(self):return bool(self.cpu.bit(4))
    @acked.setter
    def acked(self,v):self.cpu.put(4,int(v))
    def word(self,name):
        p=self.regs[name].refs[0];out=[];seen=set();zero=self.regs['token.zero'].refs[0];one=self.regs['token.one'].refs[0]
        while p is not None:
            assert p not in seen;seen.add(p);value=p.refs[1]
            assert value is zero or value is one
            out.append(int(value is one));p=p.refs[0]
        return out
    def tick(self):
        assert not self.done and not self.refused
        n=self.pc;op=OPS[n.bits&31];i=n.bits>>5;a=n.refs[0];branch=0
        # Current rights are tested for every material operand; refusal leaves
        # this same native PC and all installed fields/latches unchanged.
        required=[n,self.cpu]
        if a is not None:required.append(a)
        if op in ('rfield','wfield','bfield','vfield'):required.append(self.address)
        if any(x is None or not x.grant for x in required):self.refused=True;self.events['grant_refusals']+=1;return
        self.steps+=1;self.events[op]+=1;self.events['grant_tests']+=len(required)
        self.events['native_label_bits']+=21;self.events['native_operand_reference_reads']+=int(a is not None)
        if op in ('wfield','vfield') and self.address in self.frozen:
            raise AssertionError('write to an immutable published frame/token')
        if op=='rr':self.rr=a.refs[0]
        elif op=='rw':a.refs[0]=self.rr
        elif op=='br':self.bit=a.bit(0)
        elif op=='bw':a.put(0,self.bit)
        elif op=='address':self.address=self.rr
        elif op=='rfield':self.rr=self.address.refs[i]
        elif op=='wfield':
            if self.address is self.regs['work.root'].refs[0]:
                if i==1 and self.address.refs[1] is None and self.rr is not None:self.events['work_P']+=1
                if i==3 and self.address.refs[3] is None and self.rr is not None:self.events['work_S']+=1
            self.address.refs[i]=self.rr
        elif op=='bfield':self.bit=self.address.bit(i)
        elif op=='vfield':self.address.put(i,self.bit)
        elif op=='literal':self.bit=i
        elif op=='nil':self.rr=None
        elif op=='save':self.saved=self.rr
        elif op=='null':branch=int(self.rr is None)
        elif op=='equal':branch=int(self.rr is self.saved)
        elif op=='test':branch=self.bit
        elif op=='send':self.output.append(str(self.bit))
        elif op=='end':self.ended=True;self.acked=False
        elif op=='ack':assert self.ended;self.acked=True
        elif op=='return':assert self.acked
        elif op=='wait':self.done=True
        elif op=='fault':self.refused=True;return
        elif op!='goto':raise AssertionError(op)
        self.events['controller_selection']+=1;self.events['successor_read']+=1;self.events['control_write']+=1
        self.pc=n.refs[3 if branch else 2]
        self.cut_checks+=1
        assert ''.join(str(b.bit(0)) for b in self.received)==self.bits_input
        if self.pc in self.boundary:
            idx=self.boundary[self.pc];row=self.p.nodes[idx]
            if row.op=='copy' and row.args[0]=='pool':
                newpool=self.regs['pool'].refs[0]
                if newpool is not self.last_pool:self.allocated.append(self.last_pool);self.last_pool=newpool
            label=self.labels.get(idx,'')
            if label.startswith('q.') or label=='print':
                for name,slot in self.regs.items():
                    if name.startswith('stack.'):
                        cell=slot.refs[0]
                        while cell is not None:
                            if cell in self.frozen:break
                            self.frozen[cell]=tuple(cell.refs);cell=cell.refs[0]
                for name in ('work.root','token.one','token.zero','token.marker','route.header'):
                    cell=self.regs[name].refs[0];self.frozen[cell]=tuple(cell.refs)
    def run(self,limit=1_000_000):
        while not self.done and not self.refused and self.steps<limit:self.tick()
        assert self.done or self.refused or self.steps==limit
        return ''.join(self.output)
    def census(self):return self.inst.census()

def small_program(name):
    with shared_assembler():
        a=ci.Assembler();ci.primitives(a)
        def main(a):
            if name=='alias':
                a.g.copy('input','other');a.g.pop('input','bad','one','one');a.g.mark('one');a.g.push('input',1)
                a.g.copy('other','check');a.g.copy('input','other');a.g.clear('input');a.g.copy('check','input')
            else:
                a.call('popword',['arg.0','input'],'input');a.call('popword',['arg.1','input'],'input')
                a.call(name,['computed'], 'arg.0','arg.1');a.g.copy('computed','input')
        return prune(a.compile(main))

def check_tiny(kernel):
    rows=[]
    for name,bits,expected in [('alias','0101','0101'),('add','','101'),('mul','','011'),('add-wide','','')]:
        # Encoded actual arguments: 2 and 3, built as input words, not answers.
        if name!='alias':
            args=[(1<<130)+5,3] if name=='add-wide' else [2,3];bits=''.join(''.join('1'+b for b in bin(x)[:1:-1])+'0' for x in args)
            value=sum(args) if name.startswith('add') else args[0]*args[1];expected=bin(value)[:1:-1]
        p=small_program('add' if name=='add-wide' else name);graph=field_graph(p);e=Direct(graph,bits,2048,kernel)
        output=e.run();assert e.done and not e.refused,(name,e.steps)
        assert output==reply(1,expected),(name,output,expected)
        assert e.word('stack.0')==list(map(int,expected))
        assert e.events['work_P']==e.events['work_S']==1
        if name=='alias':
            assert e.regs['stack.0'].refs[0] is e.regs['stack.'+str(p['registers']['check'])].refs[0]
            assert e.word('stack.'+str(p['registers']['other']))==[1,1,0,1]
        # Independently reconstruct IS.4-5 from the actual full layouts.
        k,l,f,r,b=2048,len(bits),len(e.native),len(e.regs),len(e.bitregs)
        n=k+l+f+r+b+1;rr=4*(k+l+f+1)+r;bb=8*(k+l+1)+21*f+b
        c=e.census();seed=e.inst.seed_census()
        assert n==len(e.sources)
        assert (c['objects'],c['refs'],c['bits'])==(3*n+1+seed['objects'],2*rr+5*n+15+seed['refs'],2*bb+2+seed['bits'])
        assert c['total_bits']==c['objects']+c['refs']*c['reference_width']+c['bits']
        # A copied tail retains identity; no token was written after loading.
        zero=e.regs['token.zero'].refs[0];one=e.regs['token.one'].refs[0]
        assert zero.refs[1] is None and zero.refs[2] is None
        assert one.refs[1] is not None and one.refs[2] is None
        rows.append(dict(name=name,input_sha256=hashlib.sha256(bits.encode()).hexdigest(),
            output=output,steps=e.steps,field_vertices=len(graph.nodes),native_vertices=len(e.native),
            counts=dict(e.events),census=e.census(),seed_census=e.inst.seed_census(),source_blocks=len(e.sources),reference_registers=len(e.regs),bit_registers=len(e.bitregs),immutable_tokens=True,all_native_cuts_checked=e.cut_checks))
    terminals=[]
    for decision in (False,True):
        p=dict(registers={'input':0},start=0,code=[('halt',decision)])
        e=Direct(field_graph(p),'',16,kernel);out=e.run()
        assert e.done and not e.refused and out==reply(int(decision),'')
        assert e.bitregs['result.ok'].bit(0)==int(decision)
        terminals.append(dict(decision=decision,output=out,steps=e.steps))
    # Exhaustion and denied current grant are observable refusals, never a reply.
    p=small_program('alias');e=Direct(field_graph(p),'0',0,kernel);e.run();assert e.refused and not e.output
    r=Direct(field_graph(p),'0',16,kernel);old=r.pc;old.grant=0;r.tick()
    assert r.refused and r.pc is old and not r.output and r.steps==0
    return dict(rows=rows,terminal_cases=terminals,empty_pool_refusal=True,grant_refusal_same_pc=True)

def source_cases():
    rows=[]
    for name in ('positive-full','positive-partial','null-nine-distinct-left'):
        path=BASE/'constructor'/(name+'.input.bits');raw=path.read_text();n,theta,cache,cyl=parse_arithmetic(raw)
        source=prepare_source([dict(row).get(h,0) for row in cache for h in range(15)])
        port=OriginalPort(source);s=Session(n,theta,port)
        source_bits=[dict(row).get(h,0) for row in cache for h in range(15)]
        words=['','LL','LRR','LLR','LL','LR','LLL','LRL']
        words+=['R'*b+'LL' for b in range(15*n) if (b%15,source_bits[b]) in cache[b//15]]
        words+=['LL','R'*15+'LL','R'*45+'LL','L'*130,'R'*(15*n+10),'LL']
        for word in words:s.read(word)
        acquired=arithmetic_input(n,theta,s.cache_tuple(),cyl);assert acquired==raw
        assert len(s.records)==len(words) and len({id(x[2].header) for x in s.records})==len(words)
        assert sorted({v for _,v,_ in s.records})==sorted(REPLIES)
        rows.append(dict(case=name,N=n,input_bits=len(raw),input_sha256=hashlib.sha256(raw.encode()).hexdigest(),
            calls=len(words),operand_length_sum=sum(map(len,words)),receipts=len(s.records),
            source_census=source.a.census(),observer_census=s.a.census(),
            scope='unchanged hosted acquisition only; not candidate evaluation, return or native continuation'))
    return rows

def main():
    old=ep.build()
    from constructor_checks import graph_bytes
    assert graph_bytes(old)==(BASE/'constructor/evaluator.stack.bin').read_bytes()
    with shared_assembler():new=ep.build()
    pruned=prune(new)
    counts=Counter(x[0] for x in pruned['code'])
    # No full evaluator native installation or arithmetic execution here.
    graph=field_graph(pruned)
    field_counts=Counter(n.op for n in graph.nodes)
    native_bound=len(graph.nodes)+sum({'copy':2,'put':2,'null':2,'equal':4,'test':2,
        'get':4,'set':4,'read':4,'write':4,'emit':2,'goto':0,'end':1,'ack':1,'return':1,'wait':1,'fault':1}[n.op] for n in graph.nodes)
    with BitKernelContext() as kernel:tiny=check_tiny(kernel)
    paths=['constructor_ir.py','constructor_checks.py','evaluator_program.py','field_graph.py','native_install.py','framing.py','session.py',
        'constructor/evaluator.stack.bin','constructor/constructor-results.json','native-prerequisite/immutable_spine.py']
    result=dict(candidate=IDENTITY,complete=False,full_UJB_evaluator_runs=0,native_full_evaluator_installations=0,
        original_stack_states=len(old['code']),new_stack_states=len(new['code']),reachable_stack_states=len(pruned['code']),
        registers=len(pruned['registers']),stack_counts=dict(counts),field_vertices=len(graph.nodes),
        field_counts=dict(field_counts),direct_native_vertices=native_bound,reference_registers=len(graph.refs),bit_registers=len(graph.bits),
        owner_layout=dict(n='K+L+F+R+B+1',r='4*(K+L+F+1)+R',b='8*(K+L+1)+21*F+B',
            retained_objects='3*n+1+seed.O',retained_reference_slots='2*r+5*n+15+seed.R',
            retained_bits='2*b+2+seed.B',peak_scope='standalone native owner only; constant after installation; all source shadows and unused pool cells retained',
            sufficient_prefix_capacity='K >= 5+L+T+1 for T commits and one unfinished next source instruction',
            allocation_exact='5+L+number of executed source push instructions'),
        fixed_candidate_sha256=hashlib.sha256(canonical(pruned)).hexdigest(),
        tiny=tiny,original_chronologies=source_cases(),
        source_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in paths},
        remaining='No noncircular complete-session push/step bound, concrete host H/W upper bounds or supplied reservation. Full native service splice, post-evaluation original Read and independent experimental/reconciliation units unexecuted.',
        owned_jobs_joined=True)
    dest=Path(__file__).with_name('immutable-spine-results.json');dest.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('tiny','source_hashes','original_chronologies')},sort_keys=True))

class BitKernelContext:
    def __enter__(self):self.k=BitKernel();return self.k
    def __exit__(self,*args):self.k.close()

if __name__=='__main__':main()
