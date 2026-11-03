"""Installed ML-I primitive graph consuming the unchanged encoded field code.

Native instructions are exactly single field accesses, assignments, tests and
supplied local I/O/control primitives. All variable field/opcode choices are
executed bit-test trees. Native graph wiring is the base ML-I finite transition
graph; this model does not recursively simulate its Python implementation.
"""
from collections import Counter
from dataclasses import dataclass
from field_graph import Cell, OPS
from primitive_engine import PrimitiveEngine, Code, Slot
from paid_code import Page, Fork
from native_install import Block, Source, Installer, BitKernel

NOPS=('rr','rw','br','bw','address','rfield','wfield','bfield','vfield',
      'literal','nil','save','null','equal','test','shape','goto','send',
      'end','ack','return','fault','commit','codebit','coderef')
SHAPES=((4,8),(4,21),(1,0),(0,1),(0,32768),(2,0))
REFS=('pc','old.pc','a','b','zero','one','rx','ry','destination')
BITS=tuple(f'op.{i}' for i in range(5))+tuple(f'imm.{i}' for i in range(16))+('predicate','bx','occupied','waiting','ended','acknowledged')

@dataclass(eq=False,slots=True)
class Native:
    op: str
    a: object=None
    imm: int=0
    zero: object=None
    one: object=None

class Core:
    def __init__(self,r,b):
        self.nodes=[];self.r=r;self.b=b
    def node(self,op,a=None,imm=0,zero=None,one=None):
        n=Native(op,a,imm,zero,zero if one is None else one);self.nodes.append(n);return n
    def seq(self,actions,end):
        for op,a,imm in reversed(actions):end=self.node(op,a,imm,end)
        return end
    def rr(self,name):return ('rr',self.r[name],0)
    def rw(self,name):return ('rw',self.r[name],0)
    def br(self,name):return ('br',self.b[name],0)
    def bw(self,name):return ('bw',self.b[name],0)
    def addr(self,name):return [self.rr(name),('address',None,0)]
    def deref(self,name,dst):
        return self.addr(name)+[('rfield',None,0),self.rw(dst)]
    def shape(self,name,kind,yes,no):
        return self.seq([self.rr(name)],self.node('shape',imm=kind,zero=no,one=yes))
    def test(self,name,yes,no):
        return self.seq([self.br(name)],self.node('test',zero=no,one=yes))
    def const(self,name,v,end):return self.seq([('literal',None,v),self.bw(name)],end)
    def choice(self,bits,leaves,fault):
        # High-to-low complete prefix tree; all absent leaves reject.
        def tree(at,vals):
            if at==len(bits):return vals.get(0,fault)
            half=1<<(len(bits)-at-1)
            z={k:v for k,v in vals.items() if k<half}
            o={k-half:v for k,v in vals.items() if k>=half}
            if not z and not o:return fault
            return self.test(bits[at],tree(at+1,o),tree(at+1,z))
        return tree(0,leaves)
    def build(self):
        fault=self.node('fault')
        fetch=self.node('goto')
        gate=self.node('goto')
        gate.zero=gate.one=self.test('waiting',gate,fetch)
        commit=self.node('commit',zero=gate)
        writepc=lambda src:self.seq([self.rr(src),self.rw('pc')],commit)
        normal=self.test('predicate',writepc('one'),writepc('zero'))
        commit_select=self.seq([self.rr('destination')],self.node('null',zero=writepc('destination'),one=normal))
        finish=commit_select
        def dstref(actions):
            return self.shape('a',2,self.seq(self.addr('a')+actions+[('wfield',None,0)],finish),fault)
        def dstbit(actions):
            return self.shape('a',3,self.seq(self.addr('a')+actions+[('vfield',None,0)],finish),fault)
        def decision(make):
            return make(self.const('predicate',1,finish),self.const('predicate',0,finish))
        cp=dstref([self.rr('ry')])
        ops={}
        ops['copy']=self.shape('b',2,self.seq(self.deref('b','ry'),cp),fault)
        ops['bind']=self.seq([self.rr('b'),self.rw('ry')],cp)
        ops['put']=self.choice(tuple(f'imm.{i}' for i in reversed(range(16))),
                               {v:dstbit([('literal',None,v)]) for v in (0,1)},fault)
        ops['test']=self.shape('a',3,self.seq(self.addr('a')+[('bfield',None,0)],
                          decision(lambda y,n:self.node('test',zero=n,one=y))),fault)
        ops['null']=self.shape('a',2,self.seq(self.deref('a','rx')+[self.rr('rx')],
                          decision(lambda y,n:self.node('null',zero=n,one=y))),fault)
        ops['equal']=self.shape('a',2,self.shape('b',2,
                         self.seq(self.deref('a','rx')+self.deref('b','ry')+
                                  [self.rr('rx'),('save',None,0),self.rr('ry')],
                                  decision(lambda y,n:self.node('equal',zero=n,one=y))),fault),fault)
        for op,count,shape in (('get',4,0),('edge',2,5),('page',32768,4),('read',8,0),('set',4,0),('write',8,0)):
            leaves={}
            for i in range(count):
                if op in ('get','edge'):
                    leaves[i]=self.seq(self.addr('ry')+[('rfield',None,i),self.rw('ry')],cp)
                elif op in ('read','page'):
                    leaves[i]=self.seq(self.addr('ry')+[('bfield',None,i),self.bw('bx')],dstbit([self.br('bx')]))
                elif op=='set':
                    leaves[i]=self.seq(self.addr('rx')+[self.rr('ry'),('wfield',None,i)],finish)
                else:
                    leaves[i]=self.seq(self.addr('rx')+self.addr('b')+[('bfield',None,0),self.bw('bx')]+self.addr('rx')+[self.br('bx'),('vfield',None,i)],finish)
            select=self.choice(tuple(f'imm.{i}' for i in reversed(range(16))),leaves,fault)
            if op in ('get','edge','read','page'):
                ops[op]=self.shape('b',2,self.seq(self.deref('b','ry'),self.shape('ry',shape,select,fault)),fault)
            elif op=='set':
                ops[op]=self.shape('a',2,self.shape('b',2,self.seq(self.deref('a','rx')+self.deref('b','ry'),self.shape('rx',shape,select,fault)),fault),fault)
            else:
                ops[op]=self.shape('a',2,self.shape('b',3,self.seq(self.deref('a','rx'),self.shape('rx',shape,select,fault)),fault),fault)
        ops['indirect']=self.shape('a',2,self.seq(self.deref('a','rx'),self.shape('rx',1,
                                self.seq([self.rr('rx'),self.rw('destination')],finish),fault)),fault)
        ops['goto']=finish
        ops['wait']=self.const('waiting',1,finish)
        ops['fault']=fault
        ops['emit']=self.shape('a',3,self.seq(self.addr('a')+[('bfield',None,0),self.bw('bx')],
                          self.const('occupied',1,self.seq([self.br('bx'),('send',None,0)],self.const('occupied',0,finish)))),fault)
        ops['end']=self.seq([('end',None,0)],self.const('ended',1,self.const('acknowledged',0,finish)))
        ops['ack']=self.test('ended',self.seq([('ack',None,0)],self.const('acknowledged',1,finish)),fault)
        ops['return']=self.test('acknowledged',self.seq([('return',None,0)],finish),fault)
        dispatch=self.choice(tuple(f'op.{i}' for i in reversed(range(5))),{OPS.index(k):v for k,v in ops.items()},fault)
        dispatch=self.const('predicate',0,self.seq([('nil',None,0),self.rw('destination')],dispatch))
        fetch_end=dispatch
        for i,name in reversed(tuple(enumerate(('a','b','zero','one')))):
            fetch_end=self.seq(self.addr('pc')+[('coderef',None,i),self.rw(name)],fetch_end)
        for i in reversed(range(16)):
            fetch_end=self.seq(self.addr('pc')+[('codebit',None,5+i),self.bw(f'imm.{i}')],fetch_end)
        for i in reversed(range(5)):
            fetch_end=self.seq(self.addr('pc')+[('codebit',None,i),self.bw(f'op.{i}')],fetch_end)
        fetch.zero=fetch.one=self.seq([self.rr('pc'),self.rw('old.pc')],self.shape('pc',1,fetch_end,fault))
        self.entry=fetch
        return self

class SlotView:
    def __init__(self,slots,bit=False):self.slots=slots;self.bit=bit
    def __getitem__(self,name):
        o=self.slots[name];return o.bit(0) if self.bit else o.refs[0]
    def __setitem__(self,name,value):
        if self.bit:self.slots[name].put(0,value)
        else:self.slots[name].refs[0]=value
    def values(self):return (self[k] for k in self.slots)
    def __len__(self):return len(self.slots)

class NativeEngine(PrimitiveEngine):
    """Preparation is a retained literal source; first use pays its installation.

    install_b/install_input remain source preparers. No prepared target is ever
    executed: after seal, native operations read only relocated target fields.
    """
    def __init__(self,*args,kernel=None,**kwargs):
        super().__init__(*args,**kwargs)
        self.sealed=False;self.kernel=kernel.independent_header() if kernel is not None else None;self.owns_kernel=kernel is None
        self.native_steps=0;self.native_events=Counter();self.read_ref=None;self.saved_ref=None
        self.address=None;self.read_bit=0
        self.installation=None
    native_operations=NOPS
    def make_core(self,r,b,registers):
        return Core(r,b).build()
    def seal(self,checkpoint=None):
        assert not self.sealed
        # Public source graph includes every reserved cell, code, input object,
        # scalar bit/register and complete immutable image; nothing is omitted.
        source_r={k:Block([v.value],bytearray(),0) for k,v in self.registers.slots.items()}
        source_b={k:Block([],bytearray([v.value]),1) for k,v in self.bits.slots.items()}
        cpu_r={k:Block([getattr(self.cpu,k,None)],bytearray(),0) for k in REFS}
        cpu_b={k:Block([],bytearray([0]),1) for k in BITS}
        core=self.make_core(cpu_r,cpu_b,source_r)
        originals=list(self.cells)+list(self.code)+list(self.registers.slots.values())+list(self.bits.slots.values())
        if self.store:originals+=list(self.store.pages)+list(self.store.forks)
        blocks={}
        bit_slots=set(self.bits.slots.values())
        for x in originals:
            if isinstance(x,Cell):b=Block(list(x.refs),bytearray([x.bits]),8)
            elif isinstance(x,Code):b=Block([x.a,x.b,x.zero,x.one],bytearray((x.opcode|(x.immediate<<5)).to_bytes(3,'little')),21)
            elif isinstance(x,Slot):b=source_b[next(k for k,v in self.bits.slots.items() if v is x)] if x in bit_slots else source_r[next(k for k,v in self.registers.slots.items() if v is x)]
            elif isinstance(x,Page):b=Block([],bytearray(x.data),32768)
            elif isinstance(x,Fork):b=Block([x.left,x.right],bytearray(),0)
            else:raise AssertionError(type(x))
            blocks[x]=b
        for n in core.nodes:
            blocks[n]=Block([n.a,None,n.zero,n.one],bytearray((self.native_operations.index(n.op)|(n.imm<<5)).to_bytes(3,'little')),21)
        allblocks=list(dict.fromkeys(list(blocks.values())+list(cpu_r.values())+list(cpu_b.values())))
        
        for b in allblocks:b.grant=1
        sources=[Source(b) for b in allblocks];byblock=dict(zip(allblocks,sources))
        for s in sources:
            s.block.refs[:]=[byblock[v if isinstance(v,Block) else blocks[v]] if v is not None else None for v in s.block.refs]
        for a,b in zip(sources,sources[1:]):a.next=b
        if self.kernel is None:self.kernel=BitKernel()
        self.installation=Installer(sources,self.kernel,byblock[blocks[core.entry]])
        self.installation.run(checkpoint);self.installation.verify()
        target=lambda b:byblock[b].forward
        self.forward_observer={x:target(blocks[x]) for x in self.cells}
        self.registers=SlotView({k:target(v) for k,v in source_r.items()})
        self.bits=SlotView({k:target(v) for k,v in source_b.items()},True)
        self.cr={k:target(v) for k,v in cpu_r.items()};self.cb={k:target(v) for k,v in cpu_b.items()}
        self.origin={target(blocks[x]):x.origin for x in self.code}
        self.cells=[target(blocks[x]) for x in self.cells]
        self.code=[target(blocks[x]) for x in self.code]
        self.entries={k:target(blocks[v]) for k,v in self.entries.items()}
        self.cpu=None;self.store=None
        self.width=self.installation.width
        self.sealed=True
        # Observer maps, source graph and installer are retained, not runtime inputs.
        self.charges=Counter();self.ticks=0;self.microsteps=0
        return self.installation.census()
    @property
    def native_pc(self):return self.installation.execution_pc
    @native_pc.setter
    def native_pc(self,value):self.installation.execution_pc=value
    @property
    def waiting(self):return bool(self.cb['waiting'].bit(0)) if getattr(self,'sealed',False) else False
    @property
    def pc(self):return self.origin[self.cr['pc'].refs[0]] if self.sealed else self.cpu.pc.origin
    @pc.setter
    def pc(self,n):raise AssertionError('native control may not be replaced by an observer')
    def microtick(self):
        if not self.sealed:self.seal()
        n=self.native_pc;opcode=n.bits&31;imm=n.bits>>5;op=self.native_operations[opcode]
        a=n.refs[0];branch=0
        self.native_steps+=1;self.microsteps+=1;self.native_events[op]+=1
        event=self._event
        if op=='rr':self.read_ref=a.refs[0];event('reference_register_reads');event('read_reference_latch_writes')
        elif op=='rw':a.refs[0]=self.read_ref;event('reference_register_writes')
        elif op=='br':self.read_bit=a.bit(0);event('bit_register_reads');event('read_bit_latch_writes')
        elif op=='bw':a.put(0,self.read_bit);event('bit_register_writes')
        elif op=='address':self.address=self.read_ref;event('address_reference_writes')
        elif op=='save':self.saved_ref=self.read_ref;event('saved_reference_writes')
        elif op in ('rfield','coderef'):
            self.read_ref=self.address.refs[imm];event('reference_field_reads');event('read_reference_latch_writes')
            if op=='coderef':event('code_reference_reads')
        elif op=='wfield':self.address.refs[imm]=self.read_ref;event('reference_field_writes')
        elif op in ('bfield','codebit'):
            self.read_bit=self.address.bit(imm);event('bit_field_reads');event('read_bit_latch_writes')
            if op=='codebit':event('code_bit_reads')
        elif op=='vfield':self.address.put(imm,self.read_bit);event('bit_field_writes')
        elif op=='literal':self.read_bit=imm;event('literal_bit_assignments')
        elif op=='nil':self.read_ref=None;event('literal_null_assignments')
        elif op=='null':branch=int(self.read_ref is None);event('null_tests')
        elif op=='equal':branch=int(self.read_ref is self.saved_ref);event('reference_equality_tests')
        elif op=='test':branch=self.read_bit;event('bit_tests')
        elif op=='shape':branch=int(self.read_ref is not None and (len(self.read_ref.refs),self.read_ref.width)==SHAPES[imm]);event('arity_tests')
        elif op=='send':self.output.append(self.read_bit);event('output_bits')
        elif op in ('end','ack','return'):event({'end':'output_ends','ack':'acknowledgements','return':'returns'}[op])
        elif op=='commit':self.last=self.origin[self.cr['old.pc'].refs[0]];self.ticks+=1;event('field_commits')
        elif op=='goto':pass
        elif op=='fault':
            self._event('native_controller_selections');self._event('native_refusals')
            raise RuntimeError('native field graph refused illegal state or exhausted supply')
        else:raise AssertionError(op)
        self.charges['native_controller_selections']+=1
        self.charges['native_successor_reference_reads']+=1
        self.charges['native_control_reference_writes']+=1
        for kind in ('native_controller_selections','native_successor_reference_reads','native_control_reference_writes'):
            if self.capture:self.events.append(kind)
            if self.observer:self.observer(kind)
        self.native_pc=n.refs[3 if branch else 2]
    def tick(self):
        if not self.sealed:self.seal()
        assert not self.waiting
        stop=self.ticks+1
        while self.ticks<stop:self.microtick()
    def run(self,bound=10000000):
        if not self.sealed:self.seal()
        return super().run(bound)
    def resume(self):
        assert self.sealed and self.waiting
        self.cb['waiting'].put(0,0);self._event('resume_signals');self._event('bit_register_writes')
    def bw(self,name,value):
        if not self.sealed:return super().bw(name,value)
        assert value in (0,1)
        self.bits.slots[name].put(0,value);self._event('input_bits');self._event('bit_register_writes')
    def census(self):
        if not getattr(self,'sealed',False):return super().census()
        c=self.installation.census()
        return c
    def check_census(self):
        if not getattr(self,'sealed',False):return super().check_census()
        assert all(c.width==8 and len(c.refs)==4 and c.grant for c in self.cells)
        assert all(v.bit(0) in (0,1) for v in self.cb.values())
        assert all(c.width==21 and len(c.refs)==4 and c.grant for c in self.code)
        return self.census()
    def close(self):
        if self.owns_kernel and self.kernel:self.kernel.close()
