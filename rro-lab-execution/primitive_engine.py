"""Resumable code/operand/field phases for field_graph's selected finite graph.

The CPU consumes actual encoded instruction fields. The compiler's name maps
and reference interpreters are not execution inputs after installation.
Installation is still host preparation, not a represented paid installer.
"""
from collections import Counter
from dataclasses import dataclass
from field_graph import Cell, Node, OPS


@dataclass(eq=False,slots=True)
class Slot:
    value: object


class View:
    """External name binding. The CPU holds actual Slot references instead."""
    def __init__(self,values):self.slots={k:Slot(v) for k,v in values.items()}
    def __getitem__(self,k):return self.slots[k].value
    def __setitem__(self,k,v):self.slots[k].value=v
    def values(self):return (s.value for s in self.slots.values())
    def items(self):return ((k,s.value) for k,s in self.slots.items())
    def __len__(self):return len(self.slots)


@dataclass(eq=False,slots=True)
class Code:
    opcode: int
    immediate: int = 0
    a: object = None
    b: object = None
    zero: object = None
    one: object = None
    origin: object = None  # observer identity only, never consumed by CPU


@dataclass(slots=True)
class CPU:
    pc: object
    phase: str = 'opcode'
    index: int = 0
    opcode: int = 0
    immediate: int = 0
    a: object = None
    b: object = None
    zero: object = None
    one: object = None
    rx: object = None
    ry: object = None
    destination: object = None
    bx: int = 0
    by: int = 0
    predicate: int = 0
    occupied: int = 0
    waiting: int = 0
    ended: int = 0
    acknowledged: int = 0


class PrimitiveEngine:
    # phase is one of this finite set; it is not an unbounded label/counter.
    PHASES=('opcode','immediate','a','b','zero','one','dispatch','read-a-ref',
            'read-b-ref','read-a-bit','read-b-bit','literal','get','set','read',
            'write','copy','put','null','equal','test','page','emit','end','ack',
            'commit','indirect','wait','goto','fault','return')
    # CPU's 8 reference, 16+5+5+5+7 binary fields, plus its finite phase.
    CPU_REFS=8
    CPU_BITS=16+5+5+7+5

    def __init__(self,program,capacity,entry='start',bindings=None,store=None,capture=False):
        assert capacity>0
        self.program=program  # observer/linker map, not consumed by microtick
        self.store=store
        self.capture=capture
        self.registers=View({k:None for k in sorted(set(program.refs)|{'pool','nil'})})
        self.bits=View({k:0 for k in program.bits})
        self.cells=[Cell([None]*4) for _ in range(capacity)]
        for i,c in enumerate(self.cells):c.refs[0]=self.cells[i+1] if i+1<capacity else None
        self.registers['pool']=self.cells[0]
        if bindings:
            for k,v in bindings.items():self.registers[k]=v
        self.code=[Code(OPS.index(n.op),origin=n) for n in program.nodes]
        mapping=dict(zip(program.nodes,self.code))
        r,b=self.registers.slots,self.bits.slots
        for n,c in zip(program.nodes,self.code):
            op,args=n.op,n.args
            c.zero,c.one=mapping[n.zero],mapping[n.one]
            if op=='put':c.a=b[args[0]];c.immediate=args[1]
            elif op in ('test','emit'):c.a=b[args[0]]
            elif op in ('null','indirect'):c.a=r[args[0]]
            elif op in ('copy','equal'):c.a=r[args[0]];c.b=r[args[1]]
            elif op=='get':c.a=r[args[0]];c.b=r[args[1]];c.immediate=args[2]
            elif op=='set':c.a=r[args[0]];c.b=r[args[2]];c.immediate=args[1]
            elif op=='read':c.a=b[args[0]];c.b=r[args[1]];c.immediate=args[2]
            elif op=='write':c.a=r[args[0]];c.b=b[args[2]];c.immediate=args[1]
            elif op=='bind':c.a=r[args[0]];c.b=mapping[args[1]]
            elif op=='page':c.a=b[args[0]];c.b=r[args[1]];c.immediate=args[2]
            elif op=='edge':c.a=r[args[0]];c.b=r[args[1]];c.immediate=args[2]
        self.entries={k:mapping[v] for k,v in program.entries.items()}
        self.cpu=CPU(self.entries[entry])
        objects=capacity+len(self.code)+len(r)+len(b)+1+(store.census()['objects'] if store else 0)
        self.width=(objects+1).bit_length()
        self.charges=Counter()
        self.output=[];self.events=[];self.cuts=[]
        self.ticks=0;self.microsteps=0;self.last=None
        self.observer=None
        self.check_census()

    @property
    def waiting(self):return bool(self.cpu.waiting)
    @property
    def pc(self):return self.cpu.pc.origin
    @pc.setter
    def pc(self,n):
        # Author harness only; executing indirect branches consume CPU refs.
        self.cpu.pc=next(c for c in self.code if c.origin is n)
        self.cpu.phase='opcode';self.cpu.index=0

    def _event(self,kind):
        self.charges[kind]+=1
        if self.capture:self.events.append(kind)
        if self.observer:self.observer(kind)
    def _cread(self,name):self._event('code_reference_reads');return getattr(self.cpu.pc,name)
    def _rread(self,slot):self._event('reference_register_reads');return slot.value
    def _rwrite(self,slot,value):self._event('reference_register_writes');slot.value=value
    def _bread(self,slot):self._event('bit_register_reads');return slot.value
    def _bwrite(self,slot,value):
        assert value in (0,1)
        self._event('bit_register_writes');slot.value=value
    def _cpu_ref(self,name,value):self._event('cpu_reference_writes');setattr(self.cpu,name,value)
    def _cpu_bit(self,name,value):
        assert value in (0,1)
        self._event('cpu_bit_writes');setattr(self.cpu,name,value)

    def microtick(self):
        """One finite-controller subphase, with its read/return latch preserved.

        Index/phase transitions are the bounded native finite-controller graph.
        Their five-bit index and phase selection are explicit CPU storage. This
        method does not yet expand their multibit updates into single-bit cuts.
        """
        c=self.cpu
        assert not c.waiting or c.phase=='commit'
        self.microsteps+=1;self._event('controller_selections')
        phase=c.phase
        # Native finite control follows its current state. Its encoded state and
        # bounded scan index are counted; neither is a host unbounded scalar.
        if phase in ('opcode','immediate'):
            width=5 if phase=='opcode' else 16
            value=getattr(c.pc,phase)
            bit=(value>>c.index)&1
            self._event('code_bit_reads')
            old=getattr(c,phase)
            setattr(c,phase,(old&~(1<<c.index))|(bit<<c.index))
            self._event('cpu_bit_writes')
            if c.index+1==width:
                c.index=0;c.phase='immediate' if phase=='opcode' else 'a'
            else:c.index+=1
        elif phase in ('a','b','zero','one'):
            self._cpu_ref(phase,self._cread(phase))
            c.phase={'a':'b','b':'zero','zero':'one','one':'dispatch'}[phase]
        elif phase=='dispatch':
            assert c.opcode<len(OPS),'bad installed opcode'
            op=OPS[c.opcode]
            # The truth table is the CPU's fixed finite opcode dispatcher.
            assert (op not in ('get','set','edge') or c.immediate<(2 if op=='edge' else 4))
            assert (op not in ('read','write') or c.immediate<8)
            assert (op!='put' or c.immediate<2)
            assert (op!='page' or c.immediate<32768)
            self._cpu_bit('predicate',0);self._cpu_ref('destination',None)
            c.phase={'copy':'read-b-ref','null':'read-a-ref','equal':'read-a-ref',
                     'get':'read-b-ref','set':'read-a-ref','read':'read-b-ref',
                     'write':'read-a-ref','test':'read-a-bit','emit':'read-a-bit',
                     'indirect':'read-a-ref','bind':'literal','edge':'read-b-ref',
                     'page':'read-b-ref'}.get(op,op)
        elif phase=='read-a-ref':
            self._cpu_ref('rx',self._rread(c.a))
            c.phase={'null':'null','equal':'read-b-ref','set':'read-b-ref',
                     'write':'read-b-bit','indirect':'indirect'}[OPS[c.opcode]]
        elif phase=='read-b-ref':
            self._cpu_ref('ry',self._rread(c.b))
            c.phase={'copy':'copy','equal':'equal','get':'get','set':'set',
                     'read':'read','edge':'get','page':'page'}[OPS[c.opcode]]
        elif phase=='read-a-bit':
            self._cpu_bit('bx',self._bread(c.a))
            c.phase='test' if OPS[c.opcode]=='test' else 'emit'
        elif phase=='read-b-bit':self._cpu_bit('by',self._bread(c.b));c.phase='write'
        elif phase=='literal':self._cpu_ref('ry',c.b);self._event('code_reference_reads');c.phase='copy'
        elif phase=='copy':self._rwrite(c.a,c.ry);c.phase='commit'
        elif phase=='put':self._bwrite(c.a,c.immediate);c.phase='commit'
        elif phase=='null':
            self._event('null_tests');self._cpu_bit('predicate',int(c.rx is None));c.phase='commit'
        elif phase=='equal':
            self._event('reference_equality_tests');self._cpu_bit('predicate',int(c.rx is c.ry));c.phase='commit'
        elif phase=='test':self._event('bit_tests');self._cpu_bit('predicate',c.bx);c.phase='commit'
        elif phase=='get':
            if OPS[c.opcode]=='edge':
                self._event('directory_reference_reads');value=c.ry.left if c.immediate==0 else c.ry.right
            else:
                assert isinstance(c.ry,Cell)
                self._event('reference_field_reads');value=c.ry.refs[c.immediate]
            self._cpu_ref('ry',value);c.phase='copy'
        elif phase=='set':
            assert isinstance(c.rx,Cell)
            self._event('reference_field_writes');c.rx.refs[c.immediate]=c.ry;c.phase='commit'
        elif phase=='read':
            assert isinstance(c.ry,Cell)
            self._event('bit_field_reads');self._cpu_bit('bx',(c.ry.bits>>c.immediate)&1);c.phase='read.writeback'
        elif phase=='read.writeback':self._bwrite(c.a,c.bx);c.phase='commit'
        elif phase=='write':
            assert isinstance(c.rx,Cell)
            self._event('bit_field_writes');c.rx.bits=(c.rx.bits&~(1<<c.immediate))|(c.by<<c.immediate);c.phase='commit'
        elif phase=='page':
            self._event('page_bit_reads');self._cpu_bit('bx',(c.ry.data[c.immediate//8]>>(c.immediate%8))&1);c.phase='read.writeback'
        elif phase=='emit':
            if not c.occupied:self._cpu_bit('occupied',1)
            else:
                self.output.append(c.bx);self._event('output_bits');self._cpu_bit('occupied',0);c.phase='commit'
        elif phase=='end':
            self._event('output_ends');self._cpu_bit('ended',1);self._cpu_bit('acknowledged',0);c.phase='commit'
        elif phase=='ack':
            assert c.ended,'Ack before End'
            self._event('acknowledgements');self._cpu_bit('acknowledged',1);c.phase='commit'
        elif phase=='return':
            assert c.acknowledged,'Return before acknowledgement'
            self._event('returns');c.phase='commit'
        elif phase=='indirect':
            assert isinstance(c.rx,Code)
            self._cpu_ref('destination',c.rx);c.phase='commit'
        elif phase=='goto':c.phase='commit'
        elif phase=='wait':self._cpu_bit('waiting',1);c.phase='commit'
        elif phase=='fault':raise RuntimeError('field graph refused illegal state or exhausted supply')
        elif phase=='commit':
            self.last=c.pc.origin
            self._cpu_ref('pc',c.destination if c.destination is not None else c.one if c.predicate else c.zero)
            self.ticks+=1;c.index=0;c.phase='opcode'
        else:raise AssertionError(phase)
        # These writes are still compound finite-control transitions, rather
        # than an asserted one-bit primitive expansion. Both are accounted.
        self.charges['finite_control_index_write_bits']+=5
        self.charges['finite_control_phase_write_bits']+=5

    def tick(self):
        assert not self.waiting
        done=self.ticks+1
        while self.ticks<done:
            self.microtick()

    def run(self,bound=10000000):
        for _ in range(bound):
            if self.waiting:return
            self.tick()
        raise RuntimeError('external author-check budget reached')

    def resume(self):
        assert self.waiting
        self._cpu_bit('waiting',0);self._event('resume_signals')

    def bw(self,name,value):self._bwrite(self.bits.slots[name],value)
    def observe_word(self,name):return sum(self.bits[f'{name}.{i}']<<i for i in range(64))

    def census(self):
        objects=len(self.cells)+len(self.code)+len(self.registers)+len(self.bits)+1
        refs=4*len(self.cells)+4*len(self.code)+len(self.registers)+self.CPU_REFS
        bits=8*len(self.cells)+21*len(self.code)+len(self.bits)+self.CPU_BITS
        if self.store:
            rom=self.store.census();objects+=rom['objects'];refs+=rom['reference_fields'];bits+=rom['bit_fields']
        return dict(objects=objects,reference_slots=refs,bit_fields=bits,grants=objects,
                    reference_width=self.width,bits=bits+refs*self.width+objects,
                    reserved_pool_cells=len(self.cells),
                    occupied_pool_cells=sum((x.bits>>7)&1 for x in self.cells))

    def check_census(self):
        assert all(c.bits in range(256) and len(c.refs)==4 for c in self.cells)
        assert all(v in (0,1) for v in self.bits.values())
        assert self.cpu.phase in self.PHASES+('read.writeback',)
        c=self.census();assert c['objects']<(1<<self.width)
        return c
