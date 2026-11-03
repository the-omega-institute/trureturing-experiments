"""Chronological relocation of a supplied literal public graph, including pages.

Source material is retained, uncompressed ML-I input. Inflation is public-source
preparation, never an installation transition. Copied target buffers are consumed
by native_engine directly; the initializer never substitutes a bulk clone.
"""
from __future__ import annotations
from collections import Counter
import ctypes as C
import os
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent

@dataclass(eq=False, slots=True)
class Block:
    refs: list
    data: bytearray
    width: int
    next: object=None
    grant: int=0
    @property
    def bits(self):return int.from_bytes(self.data,'little')
    @bits.setter
    def bits(self,n):self.data[:]=n.to_bytes(len(self.data),'little')
    def bit(self,i):return (self.data[i//8]>>(i%8))&1
    def put(self,i,b):
        assert b in (0,1)
        at,mask=divmod(i,8);mask=1<<mask
        self.data[at]=(self.data[at]&~mask)|(b*mask)

@dataclass(eq=False, slots=True)
class Source:
    block: Block
    next: object=None
    forward: object=None
    grant: int=1

class Wire(C.Structure):pass
Wire._fields_=[('next',C.POINTER(Wire)),('byte',C.c_uint16),('mask',C.c_uint8),('action',C.c_uint8),('grant',C.c_uint8)]
class Init(C.Structure):
    _fields_=[('source',C.POINTER(C.c_uint8)),('target',C.POINTER(C.c_uint8)),
              ('cursor',C.POINTER(Wire)),
              ('latch',C.c_uint8),('done',C.c_uint8),('grant',C.c_uint8),
              ('counts',C.c_uint64*4)]

class BitKernel:
    """Fixed native field locations and a literal actual-reference cursor chain."""
    def __init__(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='rro-native-bits-')
        lib=Path(self.tmp.name)/'bits.so'
        descriptors=tuple(int(x) for x in os.environ.get('WORKTREE_SCOPE_FDS','').split(',') if x.isdigit())
        subprocess.run(['cc','-std=c11','-O3','-shared','-fPIC',str(HERE/'native_bits.c'),'-o',str(lib)],check=True,pass_fds=descriptors)
        self.lib=C.CDLL(str(lib));self.lib.initialize.argtypes=[C.POINTER(Init),C.c_uint64]
        self.lib.initialize.restype=C.c_uint64
        self.state=Init();self.state.grant=1
        self.end=Wire();self.end.action=3;self.end.grant=1
        self.wires={}
        for width in (1,3,8,21,32768):
            wires=(Wire*(3*width))()
            for i,w in enumerate(wires):
                w.byte=(i//3)//8;w.mask=1<<((i//3)%8);w.action=i%3;w.grant=1
                w.next=C.pointer(wires[i+1]) if i+1<len(wires) else C.pointer(self.end)
            self.wires[width]=wires
    def independent_header(self):
        # Independent ML-I instances share immutable wiring, with their own
        # supplied and counted initial state header. No bank binding is copied.
        k=object.__new__(BitKernel)
        k.tmp=None;k.lib=self.lib;k.end=self.end;k.wires=self.wires;k.state=Init();k.state.grant=1
        return k
    def run(self,limit=(1<<63)):return self.lib.initialize(C.byref(self.state),limit)
    def census(self,w):
        assert self.state.grant==self.end.grant==1
        assert all(w.grant==1 for a in self.wires.values() for w in a)
        wires=1+sum(len(x) for x in self.wires.values())
        return dict(objects=wires+1,refs=wires+3,bits=22*wires+2,
                    total=22*wires+2+w*(wires+3)+wires+1)
    def close(self):
        if self.tmp is not None:self.tmp.cleanup()


@dataclass(eq=False,slots=True)
class InstallNode:
    action: str
    literal: int=0
    zero: object=None
    one: object=None
    binding: object=None
    grant: int=1

LAYOUTS=((4,8),(4,21),(1,0),(0,1),(0,32768),(2,0))

class Installer:
    """Literal native initializer graph, followed by literal relocation graph.

    Each reference position is a different installed native control vertex. The
    only unbounded cursor follows Source.next. The C bit graph has three actual
    vertices per field and actual successor links. No offset/phase arithmetic is
    a transition. This initial ML-I controller is supplied finite primitive
    control, counted from the first cut, with its original current-rights premise.
    """
    def __init__(self,sources,kernel,entry=None):
        assert sources
        self.sources=sources  # independent inspection, never consulted by tick
        self.head=sources[0];self.cursor=self.head
        self.entry_descriptor=self.head if entry is None else entry;self.execution_pc=None
        self.pending=None;self.first=None;self.tail=None;self.link=None
        self.read_ref=None;self.material=None;self.target=None
        self.kernel=kernel;self.done=False;self.grant=1
        self.charges=Counter();self.steps=0;self.accounted=Counter(objects=0,refs=0,bits=0)
        self.nodes=[]
        def node(op,i=0,z=None,o=None):
            n=InstallNode(op,i,z,z if o is None else o)
            if op=='bind-cursor':n.binding=C.pointer(kernel.wires[LAYOUTS[i][1]][0])
            self.nodes.append(n);return n
        def seq(actions,end):
            for op,i in reversed(actions):end=node(op,i,end)
            return end
        done=node('done');bad=node('fault')
        handoff=seq([('entry-read',0)],node('read-null',z=node('entry-write',z=done),o=bad))
        start1=node('material');start2=node('material')
        end2=seq([('next-read',0),('cursor-write',0)],node('cursor-null',z=start2,o=handoff))
        relocate={}
        for kind,(refs,bits) in enumerate(LAYOUTS):
            end=end2
            for i in reversed(range(refs)):
                write=node('target-ref-write',i,end)
                translated=node('forward-read',z=node('forward-null',z=write,o=bad))
                source_read=node('source-ref-read',i,node('read-null',z=translated,o=write))
                end=source_read
            relocate[kind]=seq([('target-read',0)],end)
        dispatch=bad
        for kind in reversed(range(len(LAYOUTS))):dispatch=node('shape',kind,dispatch,relocate[kind])
        start2.zero=start2.one=dispatch
        reset=seq([('head-read',0),('cursor-write',0)],start2)
        end1=seq([('next-read',0),('cursor-write',0)],node('cursor-null',z=start1,o=reset))
        allocate={}
        for kind,(refs,bits) in enumerate(LAYOUTS):
            end=end1
            if bits:
                end=seq([('clear-bit-done',0),('bind-source',0),('bind-target',0),('bind-cursor',kind)],node('bits',z=end1))
            for i in reversed(range(refs)):end=node('clear-ref',i,end)
            end=node('forward-write',z=end)
            end=seq([('tail-write',0)],end)
            linked=node('link-tail',z=end)
            head=node('link-head',z=end)
            end=node('tail-null',z=linked,o=head)
            allocate[kind]=seq([('allocate',kind),('grant',0),('clear-directory',0)],end)
        dispatch=bad
        for kind in reversed(range(len(LAYOUTS))):dispatch=node('shape',kind,dispatch,allocate[kind])
        start1.zero=start1.one=dispatch;self.pc=start1
        n=len(sources);seed_nodes=1+sum(len(x) for x in kernel.wires.values())+len(self.nodes)
        self.width=(3*n+seed_nodes+64).bit_length()
        # Source Block has its own next directory slot; each Source descriptor
        # has block/next/forward. The header includes every runtime/installer latch.
        self.base=dict(objects=2*n+1,refs=sum(len(x.block.refs)+4 for x in sources)+15,
                       bits=sum(x.block.width for x in sources)+2)
        self.peak=0
    @property
    def phase(self):return self.pc.action if self.pc is not None else 'done'
    def event(self,k):self.charges[k]+=1;self.steps+=1
    def tick(self,bit_budget=1):
        assert not self.done
        n=self.pc;op=n.action;i=n.literal;branch=0
        s=self.cursor
        if op=='material':self.material=s.block;self.event('source_material_reference_reads')
        elif op=='shape':branch=int((len(self.material.refs),self.material.width)==LAYOUTS[i]);self.event('arity_tests')
        elif op=='allocate':
            refs,bits=LAYOUTS[i]
            self.pending=Block([None]*refs,bytearray((bits+7)//8),bits)
            self.accounted.update(objects=1,refs=refs+1,bits=bits)
            self.event('allocations')
        elif op=='grant':self.pending.grant=1;self.event('grant_deliveries')
        elif op=='clear-directory':self.pending.next=None;self.event('directory_reference_writes')
        elif op=='tail-null':branch=int(self.tail is None);self.event('null_tests')
        elif op=='link-tail':self.tail.next=self.pending;self.event('directory_reference_writes')
        elif op=='link-head':self.first=self.pending;self.event('initializer_reference_writes')
        elif op=='tail-write':self.tail=self.pending;self.event('initializer_reference_writes')
        elif op=='forward-write':s.forward=self.pending;self.event('forward_reference_writes')
        elif op=='clear-ref':self.pending.refs[i]=None;self.event('reference_initializations')
        elif op=='clear-bit-done':self.kernel.state.done=0;self.event('initializer_bit_writes')
        elif op=='bind-source':
            self.kernel.source_array=(C.c_uint8*len(self.material.data)).from_buffer(self.material.data)
            self.kernel.state.source=self.kernel.source_array;self.event('initializer_reference_writes')
        elif op=='bind-target':
            self.kernel.target_array=(C.c_uint8*len(self.pending.data)).from_buffer(self.pending.data)
            self.kernel.state.target=self.kernel.target_array;self.event('initializer_reference_writes')
        elif op=='bind-cursor':
            self.kernel.state.cursor=n.binding;self.event('initializer_binding_reads')
        elif op=='bits':
            state=self.kernel.state;before=tuple(state.counts)
            executed=self.kernel.run(bit_budget);self.steps+=executed
            for j,k in enumerate(('bit_initializations','source_bit_reads','target_bit_writes','bit_scan_completions')):
                self.charges[k]+=(state.counts[j]-before[j])&((1<<64)-1)
            self.charges['initializer_bit_latch_writes']+=(state.counts[1]-before[1])&((1<<64)-1)
            for k in ('native_controller_selections','native_successor_reference_reads','native_control_reference_writes'):
                self.charges[k]+=executed
            if not state.done:return
            self.event('bit_scan_returns')
        elif op=='next-read':self.link=s.next;self.event('source_directory_reads')
        elif op=='head-read':self.link=self.head;self.event('initializer_reference_reads')
        elif op=='cursor-write':self.cursor=self.link;self.event('initializer_reference_writes')
        elif op=='cursor-null':branch=int(self.cursor is None);self.event('null_tests')
        elif op=='target-read':self.target=s.forward;self.event('forward_reference_reads')
        elif op=='source-ref-read':self.read_ref=self.material.refs[i];self.event('source_reference_reads')
        elif op=='read-null':branch=int(self.read_ref is None);self.event('null_tests')
        elif op=='forward-read':self.read_ref=self.read_ref.forward;self.event('forward_reference_reads')
        elif op=='forward-null':branch=int(self.read_ref is None);self.event('null_tests')
        elif op=='target-ref-write':self.target.refs[i]=self.read_ref;self.event('target_reference_writes')
        elif op=='entry-read':self.read_ref=self.entry_descriptor.forward;self.event('forward_reference_reads')
        elif op=='entry-write':self.execution_pc=self.read_ref;self.event('execution_control_reference_writes')
        elif op=='done':self.done=True;self.event('installation_completions')
        elif op=='fault':
            self.event('installation_refusals');self.charges['native_controller_selections']+=1
            raise RuntimeError('malformed public installation source')
        else:raise AssertionError(op)
        # A primitive read returns its value into the named retained latch.
        # These coordinates count effects of that same actual transition;
        # they do not introduce another primitive or increment the step count.
        if op in ('material','next-read','head-read','target-read','source-ref-read','forward-read','entry-read'):
            self.charges['initializer_reference_latch_writes']+=1
        if op=='bind-cursor':self.charges['initializer_reference_writes']+=1
        if op=='allocate':self.charges['allocation_reference_latch_writes']+=1
        for k in ('native_controller_selections','native_successor_reference_reads','native_control_reference_writes'):self.charges[k]+=1
        self.pc=n.one if branch else n.zero
    def run(self,checkpoint=None):
        while not self.done:
            self.tick(1<<30)
            if checkpoint:checkpoint(self)
    def seed_census(self):
        c=self.kernel.census(self.width);n=len(self.nodes)
        return dict(objects=c['objects']+n,refs=c['refs']+3*n,bits=c['bits']+21*n)
    def census(self):
        assert self.grant==1 and all(x.grant==1 for x in self.nodes)
        assert all(s.grant==s.block.grant==1 for s in self.sources)
        n=r=b=0;p=self.first;seen=set()
        while p is not None:
            assert id(p) not in seen;seen.add(id(p));n+=1;r+=len(p.refs)+1;b+=p.width;p=p.next
        if self.pending is not None and id(self.pending) not in seen:
            n+=1;r+=len(self.pending.refs)+1;b+=self.pending.width
        got=Counter(objects=n,refs=r,bits=b)
        assert got==self.accounted,('installation census mismatch',got,self.accounted)
        seed=self.seed_census()
        total={k:self.base[k]+got[k]+seed[k] for k in ('objects','refs','bits')}
        total['reference_width']=self.width
        total['total_bits']=total['bits']+self.width*total['refs']+total['objects']
        self.peak=max(self.peak,total['total_bits'])
        return total
    def verify(self):
        assert self.done and self.execution_pc is self.entry_descriptor.forward
        for s in self.sources:
            t=s.forward
            assert t.data==s.block.data and t.width==s.block.width and t.grant==1
            assert all(a is (b.forward if b is not None else None) for a,b in zip(t.refs,s.block.refs))
        return self.census()
