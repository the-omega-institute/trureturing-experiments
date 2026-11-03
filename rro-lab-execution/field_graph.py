"""Direct finite field graphs for the conditional ML-I boundary.

This is a new executing path, not a relabelling of lab.Machine's host state.
All data-dependent persistent operands are bits or actual Cell references.
The finite graph itself is the control location; dictionaries/labels are linker
inputs only. A reserved finite pool has one common fixed reference width.
The interpreter is a field-graph model, not a physical-device implementation.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
import hashlib
import json


@dataclass(eq=False, slots=True)
class Cell:
    refs: list
    bits: int = 0


@dataclass(eq=False, slots=True)
class Node:
    op: str
    args: tuple
    zero: object = None
    one: object = None


# Direct graph actions; indices in args are fixed wiring, not computed offsets.
OPS = ('put','test','copy','null','equal','get','set','read','write',
       'bind','goto','indirect','wait','fault','page','emit','end','ack','edge','return')


class Graph:
    def __init__(self):
        self.rows=[]
        self.labels={}
        self.serial=0
        self.ref_names=set()
        self.bit_names=set()

    def fresh(self, prefix='k'):
        self.serial += 1
        return prefix+str(self.serial)

    def label(self, name):
        assert name not in self.labels, name
        self.labels[name]=len(self.rows)

    def emit(self, op, *args, zero=None, one=None):
        assert op in OPS
        self.rows.append([op, tuple(args), zero, one])

    def j(self, target):self.emit('goto',zero=target)
    def put(self, dst, value):
        assert value in (0,1)
        self.bit_names.add(dst);self.emit('put',dst,value)
    def test(self, bit, one, zero=None):
        self.bit_names.add(bit);self.emit('test',bit,zero=zero,one='@next' if one is None else one)
    def cp(self, dst, src):
        self.ref_names.update((dst,src));self.emit('copy',dst,src)
    def null(self, src, yes, no=None):
        self.ref_names.add(src);self.emit('null',src,zero=no,one=yes)
    def equal(self, a, b, yes, no=None):
        self.ref_names.update((a,b));self.emit('equal',a,b,zero=no,one=yes)
    def get(self, dst, src, field):
        assert 0<=field<4
        self.ref_names.update((dst,src));self.emit('get',dst,src,field)
    def set(self, dst, field, src):
        assert 0<=field<4
        self.ref_names.update((dst,src));self.emit('set',dst,field,src)
    def read(self, dst, src, field=0):
        assert 0<=field<8
        self.bit_names.add(dst);self.ref_names.add(src);self.emit('read',dst,src,field)
    def write(self, dst, field, src):
        assert 0<=field<8
        self.ref_names.add(dst);self.bit_names.add(src);self.emit('write',dst,field,src)
    def bind(self, dst, label):
        self.ref_names.add(dst);self.emit('bind',dst,label)
    def ret(self, src):
        self.ref_names.add(src);self.emit('indirect',src)
    def bitcopy(self,dst,src):
        one,done=self.fresh(),self.fresh()
        self.test(src,one);self.put(dst,0);self.j(done)
        self.label(one);self.put(dst,1);self.label(done)
    def boolean(self,dst,inputs,truth):
        done=self.fresh()
        def visit(i,word):
            if i==len(inputs):
                self.put(dst,int(truth(*word)));self.j(done);return
            one=self.fresh();self.test(inputs[i],one)
            visit(i+1,word+(0,));self.label(one);visit(i+1,word+(1,))
        visit(0,());self.label(done)
    def word(self,name):
        return tuple(f'{name}.{i}' for i in range(64))
    def constant(self,name,n):
        assert 0<=n<1<<64
        for i,b in enumerate(self.word(name)):self.put(b,(n>>i)&1)
    def wcopy(self,dst,src):
        if dst!=src:
            for a,b in zip(self.word(dst),self.word(src)):self.bitcopy(a,b)
    def alu(self,op,dst,a,b=None,imm=0):
        # Snapshotting also covers every overlapping source/destination case.
        self.wcopy('alu.a',a)
        if b is not None:self.wcopy('alu.b',b)
        aa,bb,out=self.word('alu.a'),self.word('alu.b'),self.word(dst)
        if op in ('copy','shl','shr'):
            for i,target in enumerate(out):
                at=i if op=='copy' else i-imm if op=='shl' else i+imm
                if 0<=at<64:self.bitcopy(target,aa[at])
                else:self.put(target,0)
        elif op in ('add','sub'):
            self.put('carry',0)
            for i in range(64):
                self.boolean('sum',(aa[i],bb[i],'carry'),lambda x,y,c:x^y^c)
                if op=='add':self.boolean('nextcarry',(aa[i],bb[i],'carry'),lambda x,y,c:(x&y)|(x&c)|(y&c))
                else:self.boolean('nextcarry',(aa[i],bb[i],'carry'),lambda x,y,c:((1-x)&(y|c))|(y&c))
                self.bitcopy(out[i],'sum');self.bitcopy('carry','nextcarry')
            ok=self.fresh();self.test('carry','fault',ok);self.label(ok)
        elif op in ('eq','lt'):
            self.put('predicate',int(op=='eq'));self.put('decided',0)
            for i in reversed(range(64)):
                if op=='eq':self.boolean('predicate',('predicate',aa[i],bb[i]),lambda p,x,y:p&(1^(x^y)))
                else:
                    done=self.fresh();different=self.fresh()
                    self.test('decided',done)
                    self.boolean('different',(aa[i],bb[i]),lambda x,y:x^y)
                    self.test('different',different);self.j(done)
                    self.label(different);self.bitcopy('predicate',bb[i]);self.put('decided',1)
                    self.label(done)
            self.bitcopy(out[0],'predicate')
            for i in range(1,64):self.put(out[i],0)
        else:raise AssertionError(op)

    def allocate(self,dst):
        # Extract a real reserved grant. No integer is cast into a reference.
        self.null('pool','exhausted')
        self.cp(dst,'pool');self.get('pool',dst,0)
        self.set(dst,0,'nil');self.put('one',1);self.write(dst,7,'one')

    def append(self,header,bit,prefix='append'):
        c,t=self.fresh(prefix),self.fresh(prefix)
        empty,done=self.fresh(),self.fresh()
        self.allocate(c);self.write(c,0,bit)
        self.get(t,header,2);self.null(t,empty)
        self.set(t,1,c);self.j(done)
        self.label(empty);self.set(header,1,c)
        self.label(done);self.set(header,2,c)

    def seal(self):
        assert self.rows
        nodes=[Node(op,args) for op,args,_,_ in self.rows]
        for i,(node,row) in enumerate(zip(nodes,self.rows)):
            node.zero=nodes[self.labels[row[2]]] if row[2] is not None else nodes[min(i+1,len(nodes)-1)]
            node.one=nodes[min(i+1,len(nodes)-1)] if row[3]=='@next' else nodes[self.labels[row[3]]] if row[3] is not None else node.zero
            if node.op=='bind':node.args=(node.args[0],nodes[self.labels[node.args[1]]])
        return Program(tuple(nodes),{k:nodes[v] for k,v in self.labels.items()},
                       tuple(sorted(self.ref_names)),tuple(sorted(self.bit_names)))


@dataclass(slots=True)
class Program:
    nodes: tuple
    entries: dict
    refs: tuple
    bits: tuple

    def fingerprint(self):
        ids={v:i for i,v in enumerate(self.nodes)}
        digest=hashlib.sha256()
        for n in self.nodes:
            row=(n.op,tuple(ids[x] if isinstance(x,Node) else x for x in n.args),ids[n.zero],ids[n.one])
            digest.update(json.dumps(row,separators=(',',':')).encode());digest.update(b'\n')
        digest.update(json.dumps((self.refs,self.bits),separators=(',',':')).encode())
        return digest.hexdigest()


class Engine:
    """An explicitly bounded finite supply, with arbitrary finite capacity.

    The supplied native finite transition graph is the control location. The
    trace's field access templates include its code selector and successor;
    no numeric PC, Python loop counter, audit output or reference oracle chooses
    a graph successor. tick() is a graph action, not one physical clock tick.
    The model's atomic-action boundary remains stated rather than hidden.
    """
    def __init__(self,program,capacity,entry='start',bindings=None,store=None,capture=False):
        assert capacity>0
        self.program=program
        self.store=store
        self.capture=capture
        self._code_bits=None
        self.registers={name:None for name in program.refs}
        self.bits={name:0 for name in program.bits}
        self.cells=[Cell([None]*4) for _ in range(capacity)]
        self.width=(capacity+len(program.nodes)+len(program.refs)+len(program.bits)+2+(store.census()['objects'] if store else 0)).bit_length()
        self.pc=program.entries[entry]
        self.waiting=False
        self.output=[]  # Physical environment sink; never read by execution.
        self.charges=Counter()
        self.cuts=[]
        self.events=[]
        self.ticks=0
        # These are host-initialization coordinates only. They are NOT a
        # resumable paid installer or an executing initialization program.
        self.install_trace=[]
        for i,c in enumerate(self.cells):
            c.refs[0]=self.cells[i+1] if i+1<len(self.cells) else None
            for j in range(4):self.install_trace.append(('reference',i,j))
            for j in range(8):self.install_trace.append(('bit',i,j))
        self.registers['pool']=self.cells[0]
        if bindings:
            for name,value in bindings.items():self.registers[name]=value
        self.charges['installation_reference_writes']=4*capacity+len(self.registers)+2*len(program.nodes)
        self.charges['installation_bit_writes']=8*capacity+len(self.bits)+self.code_bits()
        self.charges['grant_deliveries']=capacity+len(program.nodes)+2
        self.charges['installed_pool_objects']=capacity
        self.last=None
        self.check_census()

    def code_bits(self):
        # Finite graph layout: 5 opcode bits, fixed selectors in a canonical
        # literal encoding. Text is linker syntax, never a runtime string parse.
        if self._code_bits is None:
            self._code_bits=sum(5+8*len(json.dumps(tuple('label' if isinstance(a,Node) else a for a in n.args),separators=(',',':')).encode()) for n in self.program.nodes)
        return self._code_bits

    def rr(self,name):
        self.charges['reference_register_reads']+=1
        return self.registers[name]
    def rw(self,name,value):
        self.registers[name]=value;self.charges['reference_register_writes']+=1
    def br(self,name):
        self.charges['bit_register_reads']+=1
        return self.bits[name]
    def bw(self,name,value):
        assert value in (0,1)
        self.bits[name]=value;self.charges['bit_register_writes']+=1
    def get(self,cell,i):
        assert isinstance(cell,Cell)
        self.charges['reference_field_reads']+=1;return cell.refs[i]
    def set(self,cell,i,value):
        assert isinstance(cell,Cell)
        self.charges['reference_field_writes']+=1;cell.refs[i]=value
    def read(self,cell,i):
        assert isinstance(cell,Cell)
        self.charges['bit_field_reads']+=1;return (cell.bits>>i)&1
    def write(self,cell,i,bit):
        assert isinstance(cell,Cell) and bit in (0,1)
        cell.bits=(cell.bits&~(1<<i))|(bit<<i);self.charges['bit_field_writes']+=1

    def tick(self):
        assert not self.waiting,'waiting state requires an external resume signal'
        n=self.pc;op,args=n.op,n.args
        self.charges['controller_selections']+=1
        self.charges['code_opcode_bit_reads']+=5
        self.charges['code_selector_bit_reads']+=8*len(json.dumps(tuple('label' if isinstance(a,Node) else a for a in args),separators=(',',':')).encode())
        branch=False;destination=None
        if op=='put':self.bw(args[0],args[1])
        elif op=='test':branch=bool(self.br(args[0]));self.charges['bit_tests']+=1
        elif op=='copy':self.rw(args[0],self.rr(args[1]))
        elif op=='null':branch=self.rr(args[0]) is None;self.charges['null_tests']+=1
        elif op=='equal':branch=self.rr(args[0]) is self.rr(args[1]);self.charges['reference_equality_tests']+=1
        elif op=='get':self.rw(args[0],self.get(self.rr(args[1]),args[2]))
        elif op=='set':self.set(self.rr(args[0]),args[1],self.rr(args[2]))
        elif op=='read':self.bw(args[0],self.read(self.rr(args[1]),args[2]))
        elif op=='write':self.write(self.rr(args[0]),args[1],self.br(args[2]))
        elif op=='bind':self.rw(args[0],args[1]);self.charges['code_reference_reads']+=1
        elif op=='goto':pass
        elif op=='indirect':destination=self.rr(args[0]);assert isinstance(destination,Node)
        elif op=='edge':
            fork=self.rr(args[1]);self.charges['directory_reference_reads']+=1
            self.rw(args[0],fork.left if args[2]==0 else fork.right)
        elif op=='page':
            page=self.rr(args[1]);i=args[2]
            self.bw(args[0],(page.data[i//8]>>(i%8))&1);self.charges['page_bit_reads']+=1
        elif op=='emit':self.output.append(self.br(args[0]));self.charges['output_bits']+=1
        elif op=='end':self.charges['output_ends']+=1
        elif op=='ack':self.charges['acknowledgements']+=1
        elif op=='return':self.charges['returns']+=1
        elif op=='wait':self.waiting=True
        elif op=='fault':raise RuntimeError(args[0] if args else 'graph fault')
        else:raise AssertionError(op)
        self.charges['code_reference_reads']+=1
        self.charges['control_reference_writes']+=1
        self.pc=destination if destination is not None else n.one if branch else n.zero
        self.ticks+=1;self.last=n
        if self.capture:self.events.append(op)

    def resume(self):
        assert self.waiting
        self.waiting=False;self.charges['resume_signals']+=1

    def run(self,bound=10000000):
        # This bound belongs to the experiment scheduler, never to the machine.
        for _ in range(bound):
            if self.waiting:return
            self.tick()
        raise RuntimeError('external author-check budget reached')

    def census(self):
        # Counts actual retained graph, pool, register and control locations.
        refs=4*len(self.cells)+2*len(self.program.nodes)+len(self.registers)+1
        refs+=sum(n.op=='bind' for n in self.program.nodes)
        bits=8*len(self.cells)+len(self.bits)+self.code_bits()+1
        objects=len(self.cells)+len(self.program.nodes)+2
        if self.store:
            rom=self.store.census();objects+=rom['objects'];refs+=rom['reference_fields'];bits+=rom['bit_fields']
        return dict(objects=objects,reference_slots=refs,bit_fields=bits,
                    grants=objects,reference_width=self.width,
                    bits=bits+refs*self.width+objects,
                    reserved_pool_cells=len(self.cells),
                    occupied_pool_cells=sum((c.bits>>7)&1 for c in self.cells))

    def check_census(self):
        assert all(c.bits in range(256) and len(c.refs)==4 for c in self.cells)
        assert all(v in (0,1) for v in self.bits.values())
        c=self.census();assert c['objects']<(1<<self.width)
        return c

    def observe_word(self,name):
        return sum(self.bits[f'{name}.{i}']<<i for i in range(64))


def growth(g):
    # Cell layout: parent, left, right, routing-header. Header/list layout uses
    # left=head/next, right=tail. Bit 7 is the actual allocated/grant flag.
    g.label('growP')
    g.get('growth.left','target',1);g.null('growth.left','P.enabled');g.j('fault')
    g.label('P.enabled')
    g.get('growth.right','target',2);g.null('growth.right','P.reserve');g.j('fault')
    g.label('P.reserve')
    g.allocate('fresh.a');g.set('fresh.a',0,'target');g.set('target',1,'fresh.a');g.j('maintain')
    g.label('growS')
    g.get('growth.left','target',1);g.null('growth.left','fault')
    g.get('growth.right','target',2);g.null('growth.right','S.reserve');g.j('fault')
    g.label('S.reserve')
    g.allocate('fresh.a');g.allocate('fresh.b');g.allocate('fresh.header')
    g.set('fresh.a',0,'growth.left');g.set('fresh.b',0,'target')
    g.set('growth.left',1,'fresh.a');g.set('target',2,'fresh.b');g.set('target',3,'fresh.header')
    g.label('maintain')
    g.cp('growth.child','target')
    g.label('maintain.loop')
    g.get('growth.parent','growth.child',0);g.null('growth.parent','maintain.done')
    g.get('growth.left','growth.parent',1)
    g.equal('growth.child','growth.left','maintain.left')
    g.get('growth.right','growth.parent',2);g.equal('growth.child','growth.right','maintain.right');g.j('fault')
    g.label('maintain.left');g.put('growth.side',0);g.j('maintain.append')
    g.label('maintain.right');g.put('growth.side',1)
    g.label('maintain.append')
    g.get('growth.header','growth.parent',3);g.null('growth.header','fault')
    g.append('growth.header','growth.side','route')
    g.cp('growth.child','growth.parent');g.j('maintain.loop')
    g.label('maintain.done');g.ret('growth.return')


def call_growth(g,kind,target):
    ret=g.fresh('growth.return')
    g.cp('target',target);g.bind('growth.return',ret);g.j('grow'+kind);g.label(ret)


def loader(g,ready='ready'):
    g.label('load')
    g.allocate('root');call_growth(g,'P','root');call_growth(g,'S','root')
    g.get('load.frontier','root',2);call_growth(g,'P','load.frontier')
    g.label('load.loop')
    call_growth(g,'S','load.frontier')
    g.null('input','load.finish')
    g.read('load.bit','input');g.get('input','input',1)
    zero=g.fresh();g.test('load.bit','load.mark',zero)
    g.label(zero);g.j('load.next')
    g.label('load.mark');g.get('load.bitref','load.frontier',2);call_growth(g,'P','load.bitref')
    g.label('load.next');g.null('input','load.finish')
    g.get('load.frontier','load.frontier',1);g.j('load.loop')
    g.label('load.finish');g.get('head','root',2);g.put('arm',0);g.j(ready)


def driver(g,image=False,fetch='fetch',halt_continuation=None):
    g.label('ready');g.emit('wait')
    g.label('step')
    g.cp('old.head','head');g.cp('new.head','head');g.bitcopy('new.arm','arm')
    if image:
        g.constant('pc.one',1);g.alu('add','nextpc','bpc','pc.one')
        g.wcopy('input_pc','bpc');g.bind('fetch.return','fetched');g.j(fetch)
    else:
        g.null('bpc','halted');g.cp('old.bpc','bpc')
        for i in range(3):g.read(f'out_op.{i}','bpc',i)
        g.get('new.bpc','bpc',1)
    g.label('fetched')
    # Decode the actual 3-bit instruction, preserving all pending old/new state.
    for value,label in ((0,'halted'),(1,'moveL'),(2,'moveR'),(3,'mark'),(4,'jump')):
        mismatch=g.fresh()
        for i in range(3):
            if (value>>i)&1:g.test(f'out_op.{i}',None,mismatch)
            else:g.test(f'out_op.{i}',mismatch)
        g.j(label);g.label(mismatch)
    g.j('fault')
    g.label('mark');g.get('bitref','head',2);g.get('bitmark','bitref',1)
    g.null('bitmark','mark.zero');g.j('commit')
    g.label('mark.zero');call_growth(g,'P','bitref');g.j('commit')
    g.label('jump');g.get('bitref','head',2);g.get('bitmark','bitref',1)
    g.null('bitmark','commit')
    if image:g.wcopy('nextpc','out_dest')
    else:g.get('new.bpc','bpc',2)
    g.j('commit')
    g.label('moveL');g.test('arm','outward','inward')
    g.label('moveR');g.test('arm','inward','outward')
    g.label('outward');g.get('new.head','head',1);g.j('ensure')
    g.label('inward');g.get('new.head','head',0);g.equal('new.head','root','cross');g.j('ensure')
    g.label('cross');g.test('arm','cross.plus')
    g.put('new.arm',1);g.get('new.head','root',1);g.j('ensure')
    g.label('cross.plus');g.put('new.arm',0);g.get('new.head','root',2)
    g.label('ensure');g.get('move.right','new.head',2);g.null('move.right','extend');g.j('commit')
    g.label('extend');call_growth(g,'S','new.head')
    g.label('commit')
    g.cp('head','new.head');g.bitcopy('arm','new.arm')
    if image:g.wcopy('bpc','nextpc')
    else:g.cp('bpc','new.bpc')
    g.j('ready')
    g.label('halted');g.emit('wait');g.j(halt_continuation or 'halted')


def b_graph():
    g=Graph();g.label('start');g.j('load')
    loader(g);driver(g);growth(g)
    g.label('fault');g.emit('fault','illegal field graph state')
    g.label('exhausted');g.emit('fault','finite supplied pool exhausted')
    return g.seal()


def install_b(program,bits,capacity=4096,graph=None,engine_type=Engine):
    # Actual data installation (outside the execution graph, explicitly not yet
    # a paid chronological installer). It never supplies fetched answers.
    p=graph or b_graph();e=engine_type(p,capacity)
    pool=e.registers['pool']
    def take():
        nonlocal pool
        assert pool is not None
        c,pool=pool,pool.refs[0];c.refs[0]=None;c.bits|=128
        return c
    cells=[take() for _ in program]
    for i,(c,ins) in enumerate(zip(cells,program)):
        value=4 if isinstance(ins,tuple) else {'L':1,'R':2,'M':3}[ins]
        c.bits|=value
        c.refs[1]=cells[i+1] if i+1<len(cells) else None
        if value==4:
            assert 0<=ins[1]<len(cells);c.refs[2]=cells[ins[1]]
    first=last=None
    for b in bits:
        assert b in '01';c=take();c.bits|=int(b)
        if last:last.refs[1]=c
        else:first=c
        last=c
    if 'bpc' in p.refs:e.registers['bpc']=cells[0] if cells else None
    e.registers['input']=first;e.registers['pool']=pool
    return e,tuple(cells)


def observe_b(e,code):
    # Independent decoder only; never consumed by graph execution.
    root=e.registers['root'];values={}
    for arm,key in ((0,2),(1,1)):
        c=root.refs[key];j=0
        while c.refs[2] is not None:
            values[-j-1 if arm else j]=int(c.refs[2].refs[1] is not None)
            c=c.refs[1];j+=1
    h=e.registers['head'];path=[]
    while h is not root:
        parent=h.refs[0];path.append(int(parent.refs[2] is h));h=parent
    coordinate=-len(path) if e.bits['arm'] else len(path)-1
    pc=len(code) if e.registers['bpc'] is None else code.index(e.registers['bpc'])
    return coordinate,pc,sorted(k for k,v in values.items() if v)
