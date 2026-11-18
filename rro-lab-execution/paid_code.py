"""Resumable paid decoder for the complete A block image.

The fixed controller uses binary registers and actual directory references.
It is a component, not an implementation of the complete ML-Q/LAB session.
No constructor_image reader is on its executing path. File decompression is
an offline preparation of the chosen uncompressed immutable installation.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import json

U64 = (1 << 64)-1
PAGE = 4096
OPS = ('const','copy','add','sub','shl','shr','eq','lt','load8','load32',
       'load64','branch','jump','finish','wait','fault')


@dataclass(eq=False, slots=True)
class Page:
    data: bytes


@dataclass(eq=False, slots=True)
class Fork:
    left: object
    right: object


class CodeStore:
    """Immutable initialized pages under an actual balanced pointer directory.

    Unused rightmost pages alias one initialized all-zero page. Each allocated
    page has 32768 binary fields. No null or third bit value occupies a field.
    Directory objects and initialization bits are included in the census.
    """
    def __init__(self, data: bytes):
        if not data:
            raise ValueError('empty installation')
        self.length = len(data)
        self.initialized_objects = 0
        self.initialized_bits = 0
        self.initialized_references = 0
        leaves = [Page(data[i:i+PAGE].ljust(PAGE,b'\0'))
                  for i in range(0,len(data),PAGE)]
        self.real_pages = len(leaves)
        power = 1 << (len(leaves)-1).bit_length()
        self.height = (power-1).bit_length()
        self.zero = None
        if power != len(leaves):
            self.zero = Page(bytes(PAGE))
            leaves += [self.zero]*(power-len(leaves))
        self.pages = tuple(dict.fromkeys(leaves))
        for page in self.pages:
            self.initialized_objects += 1
            self.initialized_bits += 8*len(page.data)
        self.forks = []
        while len(leaves)>1:
            parents = [Fork(leaves[i],leaves[i+1]) for i in range(0,len(leaves),2)]
            self.forks.extend(parents)
            self.initialized_objects += len(parents)
            self.initialized_references += 2*len(parents)
            leaves = parents
        self.root = leaves[0]
        self.raw_sha256 = hashlib.sha256(data).hexdigest()

    def census(self):
        # Inspection reads actual retained objects, not the incremental ledger.
        return {'logical_bytes':self.length,'allocated_pages':len(self.pages),
                'directory_objects':len(self.forks),'height':self.height,
                'objects':len(self.pages)+len(self.forks),
                'bit_fields':sum(8*len(p.data) for p in self.pages),
                'reference_fields':sum(2 for f in self.forks)}

    def bit_at_leaf(self, leaf, address, within_byte):
        # Bounded field selector: twelve address bits plus three bit-position
        # bits select one of the fixed 32768 physical fields of this page.
        return (leaf.data[address & (PAGE-1)] >> within_byte) & 1


class Graph:
    def __init__(self):
        self.rows=[]
        self.labels={}
        self.fixups=[]

    def label(self, name):
        if name in self.labels:
            raise ValueError(name)
        self.labels[name]=len(self.rows)

    def emit(self, op, dst='', a='', b='', imm=0, target=None):
        self.rows.append([op,dst,a,b,imm,None,target])

    def c(self,r,k):self.emit('const',r,imm=k)
    def cp(self,r,a):self.emit('copy',r,a)
    def add(self,r,a,b):self.emit('add',r,a,b)
    def sub(self,r,a,b):self.emit('sub',r,a,b)
    def shl(self,r,a,k):self.emit('shl',r,a,imm=k)
    def shr(self,r,a,k):self.emit('shr',r,a,imm=k)
    def eq(self,r,a,b):self.emit('eq',r,a,b)
    def lt(self,r,a,b):self.emit('lt',r,a,b)
    def load(self,r,a,width):self.emit('load'+str(width),r,a)
    def j(self,label):self.emit('jump',target=label)
    def br(self,r,label):self.emit('branch',a=r,target=label)

    def mul(self,r,a,k):
        # Distinct temporaries: source is never overwritten before its final use.
        self.c(r,0)
        for j in range(k.bit_length()):
            if (k>>j)&1:
                self.shl('mul_shift',a,j)
                self.add('mul_sum',r,'mul_shift')
                self.cp(r,'mul_sum')

    def plus(self,r,a,k):
        self.c('constant',k)
        self.add(r,a,'constant')

    def is_(self,a,k,target):
        self.c('constant',k)
        self.eq('test',a,'constant')
        self.br('test',target)

    def done(self,op,dest=None):
        self.c('out_op',op)
        if dest is None:self.c('out_dest',0)
        else:self.cp('out_dest',dest)
        self.emit('finish')

    def seal(self):
        for i,row in enumerate(self.rows):
            row[5] = i+1 if i+1<len(self.rows) else i
            if row[6] is not None:row[6]=self.labels[row[6]]
        self.rows=tuple(tuple(x) for x in self.rows)
        names={v for row in self.rows for v in row[1:4] if v}
        self.register_names=tuple(sorted(names | {'input_pc'}))
        return self

    def manifest(self):
        raw=json.dumps({'rows':self.rows,'registers':self.register_names,
                        'labels':self.labels},sort_keys=True,separators=(',',':')).encode()
        return {'states':len(self.rows),'registers':len(self.register_names),
                'serialized_bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),
                'opcode_set':list(OPS)}


def decoder_program():
    g=Graph()
    g.label('boot')
    g.c('address',8);g.load('n','address',32)
    g.c('address',12);g.load('nm','address',64)
    g.c('at',20);g.c('i',0)
    g.label('boot_loop')
    g.load('size','at',32)
    g.mul('twenty','size',20)
    g.plus('tmp','at',4);g.add('at','tmp','twenty')
    g.plus('i_next','i',1);g.cp('i','i_next')
    g.is_('i',16,'boot_end');g.j('boot_loop')
    g.label('boot_end')
    g.cp('base_start','at')
    g.plus('tmp','n',1);g.mul('eight','tmp',8)
    g.add('binding_start','base_start','eight')
    g.mul('halt','nm',13)
    g.plus('stopped','halt',1)
    g.emit('wait')
    g.label('decode')
    g.eq('test','input_pc','stopped');g.br('test','stopped')
    g.eq('test','input_pc','halt');g.br('test','last_mark')
    g.lt('test','input_pc','halt');g.br('test','divide')
    g.emit('fault')
    g.label('divide')
    g.c('q',0);g.c('j',0);g.c('rem',0)
    # Constant-size restoring division by 13, with every intermediate in a
    # named register. No host divmod supplies q or j.
    for bit in reversed(range(64)):
        g.shl('rem_shift','rem',1)
        g.shr('input_shift','input_pc',bit)
        g.shl('input_low','input_shift',63);g.shr('input_bit','input_low',63)
        g.add('rem','rem_shift','input_bit')
        g.c('thirteen',13);g.lt('test','rem','thirteen')
        keep='div_keep_'+str(bit);g.br('test',keep)
        g.sub('rem_new','rem','thirteen');g.cp('rem','rem_new')
        g.c('place',1<<bit);g.add('q_new','q','place');g.cp('q','q_new')
        g.label(keep)
    g.cp('j','rem')
    g.is_('j',0,'plain_right')
    g.is_('j',1,'inspect_jump')
    g.c('eight',8);g.lt('test','j','eight');g.br('test','read_zero')
    g.c('read',1);g.j('find')
    g.label('read_zero');g.c('read',0)
    g.label('find')
    g.c('lo',0);g.cp('hi','n')
    g.label('search')
    g.sub('gap','hi','lo')
    g.is_('gap',1,'found')
    g.add('sum','lo','hi');g.shr('mid','sum',1)
    g.mul('eight','mid',8);g.add('address','base_start','eight')
    g.load('base_mid','address',64)
    g.lt('test','q','base_mid');g.br('test','search_left')
    g.cp('lo','mid');g.j('search')
    g.label('search_left');g.cp('hi','mid');g.j('search')
    g.label('found')
    g.mul('eight','lo',8);g.add('address','base_start','eight')
    g.load('base','address',64)
    g.mul('seventeen','lo',17);g.add('binding','binding_start','seventeen')
    g.load('kind','binding',8)
    g.plus('address','binding',1);g.load('d0','address',64)
    g.plus('address','binding',9);g.load('d1','address',64)
    g.c('tpl',20);g.c('i',0)
    g.label('template')
    g.load('size','tpl',32)
    g.eq('test','i','kind');g.br('test','template_found')
    g.mul('twenty','size',20);g.plus('tmp','tpl',4);g.add('tpl','tmp','twenty')
    g.plus('i_next','i',1);g.cp('i','i_next')
    g.j('template')
    g.label('template_found')
    g.sub('offset','q','base')
    g.lt('test','offset','size');g.br('test','rule');g.emit('fault')
    g.label('rule')
    g.mul('twenty','offset',20);g.plus('tmp','tpl',4)
    g.add('address0','tmp','twenty');g.mul('ten','read',10)
    g.add('address','address0','ten');g.load('write','address',8)
    g.plus('address','address',1);g.load('direction','address',8)
    g.plus('address','address',1);g.load('destination','address',64)
    g.is_('destination',U64,'dest_zero')
    g.is_('destination',U64-1,'dest_one')
    g.add('dest_local','base','destination');g.cp('destination','dest_local');g.j('dest_ready')
    g.label('dest_zero');g.cp('destination','d0');g.j('dest_ready')
    g.label('dest_one');g.cp('destination','d1')
    g.label('dest_ready')
    g.mul('jump_target','destination',13)
    # Complete finite thirteen-entry block selector, including both move cases.
    for read in (0,1):
        for write in (0,1):
            for direction in (-1,1):
                label=f'block_{read}_{write}_{direction}'
                skip=label+'_skip'
                g.c('constant',read);g.eq('test','read','constant')
                g.c('zero',0);g.eq('test_inv','test','zero');g.br('test_inv',skip)
                g.c('constant',write);g.eq('test','write','constant')
                g.eq('test_inv','test','zero');g.br('test_inv',skip)
                g.c('constant',direction&255);g.eq('test','direction','constant')
                g.eq('test_inv','test','zero');g.br('test_inv',skip)
                if read==0:
                    block=(('R','M','M','M','M') if write==0 else ('M','R','M','M','M')) if direction==1 else (
                        ('L','L','L','M','M') if write==0 else ('M','L','L','L','M'))
                    start=2
                else:
                    block=('R','M','M','M') if direction==1 else ('L','L','L','M')
                    start=8
                for j,op in enumerate(block,start):
                    g.is_('j',j,label+'_'+str(j))
                g.done(4,'jump_target')
                for j,op in enumerate(block,start):
                    g.label(label+'_'+str(j));g.done({'L':1,'R':2,'M':3}[op])
                g.label(skip)
    g.emit('fault')
    g.label('inspect_jump')
    g.mul('target_base','q',13);g.plus('jump_target','target_base',8);g.done(4,'jump_target')
    g.label('plain_right');g.done(2)
    g.label('last_mark');g.done(3)
    g.label('stopped');g.done(0)
    return g.seal()


@dataclass(eq=False, slots=True)
class Register:
    value: int = 0


@dataclass(eq=False, slots=True)
class Instruction:
    op: str
    dst: object
    a: object
    b: object
    imm: int
    next: object = None
    alt: object = None


class Decoder:
    """Every tick is one resumable microstep of an installed finite graph.

    Audit lists and charges are removable. Runtime accesses only instruction
    fields, binary registers, directory links, actual initialized page fields,
    and the represented continuation fields listed in state_map().
    """
    CONTROL_WIDTHS = {
        'phase':4,'index':7,'carry':1,'decided':1,'predicate':1,
        'load_width':7,'load_bit':7,'page_depth':7,'within_byte':3,
        'increment_bit':7,'emit_bit':7,'wire_latch':1,'wire_occupied':1,
        'input_occupied':1,'ended':1,'acknowledged':1,
    }
    PHASES=('idle','execute','alu','load','increment','emit','end','ack','return','fault')

    def __init__(self,store, graph=None):
        self.store=store
        self.graph=graph or decoder_program()
        self.regs={n:Register() for n in self.graph.register_names}
        self.code=[]
        for op,d,a,b,imm,_,_ in self.graph.rows:
            self.code.append(Instruction(op,self.regs.get(d),self.regs.get(a),
                                         self.regs.get(b),imm))
        for ins,row in zip(self.code,self.graph.rows):
            ins.next=self.code[row[5]]
            ins.alt=self.code[row[6]] if row[6] is not None else None
        self.pc=self.code[self.graph.labels['boot']]
        self.entry=self.code[self.graph.labels['decode']]
        self.input_register=self.regs['input_pc']
        self.opcode_register=self.regs['out_op']
        self.destination_register=self.regs['out_dest']
        self.saved_pc=None
        self.node=None
        self.selected=None
        self.load_address=Register()
        self.phase='execute'
        for k in self.CONTROL_WIDTHS:
            if k!='phase':setattr(self,k,0)
        self.audit=[]
        self.transcript=[]
        self.charges=Counter()
        self.output_acknowledgements=0  # Audit counter only.
        self.ticks_last=0
        count=store.initialized_objects+len(self.code)+len(self.regs)+3
        self.reference_width=max(1,count.bit_length())
        self.accounted_objects=store.initialized_objects
        self.accounted_bit_fields=store.initialized_bits
        self.accounted_reference_slots=store.initialized_references
        for register in (*self.regs.values(),self.load_address):
            self.accounted_objects += 1
            self.accounted_bit_fields += 64
        for instruction in self.code:
            self.accounted_objects += 1
            self.accounted_bit_fields += 72
            self.accounted_reference_slots += 5
        # Store header: root, length and height. CPU header: nine bound
        # references, and all explicitly named continuation bit fields.
        self.accounted_objects += 2
        self.accounted_reference_slots += 10
        self.accounted_bit_fields += 71+sum(self.CONTROL_WIDTHS.values())
        self.accounted_bits=(self.accounted_bit_fields + self.accounted_objects +
                             self.accounted_reference_slots*self.reference_width)
        self.charges.update(self.installation_tariff())
        self.run()

    @property
    def ready(self):return self.phase=='idle'

    def installation_tariff(self):
        c=self.census()
        # Known fixed precompiled installation, copied in full. This is an
        # exact bulk initialization tariff, not a simulated runtime inflater.
        return {'installed_bit_initializations':c['bit_fields'],
                'installed_reference_initializations':c['reference_slots'],
                'installation_reference_bits':c['reference_slots']*self.reference_width,
                'grant_deliveries':c['objects'],
                'installed_code_payload_bits':self.store.census()['bit_fields']}

    def census(self):
        rom=self.store.census()
        objects=rom['objects']+sum(1 for _ in self.code)+sum(1 for _ in self.regs.values())+3
        refs=rom['reference_fields']+sum(5 for _ in self.code)+10
        bits=(rom['bit_fields']+sum(72 for _ in self.code)
              +sum(64 for _ in (*self.regs.values(),self.load_address))
              +71+sum(self.CONTROL_WIDTHS.values()))
        w=getattr(self,'reference_width',max(1,objects.bit_length()))
        return {'objects':objects,'reference_slots':refs,'bit_fields':bits,
                'grant_bits':objects,'reference_width':w,
                'bits':bits+refs*w+objects}

    def check_census(self):
        c=self.census()
        assert c['objects']==self.accounted_objects
        assert c['reference_slots']==self.accounted_reference_slots
        assert c['bit_fields']==self.accounted_bit_fields
        assert c['bits']==self.accounted_bits,('independent census',c,self.accounted_bits)
        assert all(0<=r.value<=U64 for r in (*self.regs.values(),self.load_address))
        assert self.phase in self.PHASES
        for n,w in self.CONTROL_WIDTHS.items():
            if n!='phase':assert 0<=getattr(self,n)<1<<w,(n,getattr(self,n))
        if self.phase not in ('idle','fault'):assert self.pc in self.code

    def state_map(self):
        return {'registers':list(self.graph.register_names),
                'register_width':64,'continuation_fields':self.CONTROL_WIDTHS,
                'reference_registers':['pc','saved_pc','node','selected','entry','input_register',
                                      'opcode_register','destination_register','load_address','store.root'],
                'instruction_fields':['opcode:8','immediate:64','dst','a','b','next','alt'],
                'immutable_pages':'32768 binary fields per initialized page',
                'directory':'two actual child references per Fork',
                'input':'64 paid bits into input_pc through occupied latch',
                'output':'67 actual latch bits, End, Ack, Return',
                'excluded_audit':['audit','transcript','charges','ticks_last',
                                  'output_acknowledgements','graph builder names/labels']}

    def _control(self,name,value):
        width=self.CONTROL_WIDTHS[name]
        if name=='phase':
            assert value in self.PHASES
        else:
            assert 0<=value<1<<width,(name,value,width)
        self.charges['control_bit_writes']+=width
        setattr(self,name,value)

    def _ref(self,name,value):
        self.charges['reference_register_write_bits']+=self.reference_width
        setattr(self,name,value)

    def _rb(self,r,i):
        self.charges['register_bit_reads']+=1
        self.charges['bounded_register_selection_tests']+=6
        return (r.value>>i)&1

    def _wb(self,r,i,v):
        assert v in (0,1)
        r.value=(r.value&~(1<<i))|(v<<i)
        self.charges['register_bit_writes']+=1
        self.charges['bounded_register_selection_tests']+=6

    def _next(self,ins):
        self.charges['code_reference_read_bits']+=self.reference_width
        self._ref('pc',ins.next)
        self._control('phase','execute')

    def _begin(self,ins):
        self._control('index',0)
        self._control('carry',0)
        self._control('decided',0)
        self._control('predicate',1 if ins.op=='eq' else 0)
        self._control('phase','alu')

    def start(self,pc):
        if not self.ready:raise RuntimeError('decoder busy')
        if not isinstance(pc,int) or not 0<=pc<=U64:raise ValueError('64-bit PC required')
        self._control('input_occupied',1)
        for i in range(64):
            self._control('wire_latch',(pc>>i)&1)
            self._control('wire_occupied',1)
            self.charges['input_bits']+=1
            self._wb(self.input_register,i,self.wire_latch)
            self._control('wire_occupied',0)
        self._control('input_occupied',0)
        self._control('ended',0);self._control('acknowledged',0)
        self._ref('pc',self.entry)
        self._control('phase','execute')
        self.ticks_last=0

    def tick(self):
        if self.phase=='idle':raise RuntimeError('no pending fetch')
        if self.phase=='fault':raise RuntimeError('invalid code address or image')
        self.ticks_last+=1
        self.charges['microsteps']+=1
        # Decode the represented finite phase before using its transition.
        self.charges['control_bit_reads']+=self.CONTROL_WIDTHS['phase']
        # Deliberate complete scans give a literal tariff even for a field
        # unused in this transition. No minimal read count is claimed.
        self.charges['control_bit_reads']+=sum(self.CONTROL_WIDTHS.values())-4
        self.charges['reference_register_read_bits']+=10*self.reference_width
        self.charges['store_header_bit_reads']+=71
        ins=self.pc
        if self.phase=='execute':
            self.charges['code_bit_reads']+=72
            self.charges['code_reference_read_bits']+=3*self.reference_width
            op=ins.op
            if op in ('const','copy','add','sub','shl','shr','eq','lt'):
                self._begin(ins)
            elif op.startswith('load'):
                self._control('load_width',int(op[4:]))
                self._control('load_bit',0);self._control('within_byte',0)
                # Copy is bitwise, with the caller PC retained.
                self._control('index',0)
                self._control('phase','load')
                self._ref('node',None)
            elif op=='branch':
                target=ins.alt if self._rb(ins.a,0) else ins.next
                self.charges['bit_tests']+=1
                self.charges['code_reference_read_bits']+=self.reference_width
                self._ref('pc',target)
            elif op=='jump':
                self.charges['code_reference_read_bits']+=self.reference_width
                self._ref('pc',ins.alt)
            elif op=='finish':
                self._control('emit_bit',0);self._control('phase','emit')
            elif op=='wait':
                self._control('phase','idle')
            elif op=='fault':
                self._control('phase','fault')
            else:raise AssertionError(op)
        elif self.phase=='alu':
            i=self.index;op=ins.op
            self.charges['control_bit_reads']+=7+1+1+1
            if op=='const':v=(ins.imm>>i)&1
            elif op=='copy':v=self._rb(ins.a,i)
            elif op in ('shl','shr'):
                at=i-ins.imm if op=='shl' else i+ins.imm
                v=self._rb(ins.a,at) if 0<=at<64 else 0
            elif op in ('add','sub'):
                x,y=self._rb(ins.a,i),self._rb(ins.b,i)
                if op=='add':
                    v=x^y^self.carry;c=(x&y)|(x&self.carry)|(y&self.carry)
                else:
                    v=x^y^self.carry;c=((1-x)&(y|self.carry))|(y&self.carry)
                self.charges['boolean_gate_evaluations']+=5
                self._control('carry',c)
            else:
                x,y=self._rb(ins.a,63-i),self._rb(ins.b,63-i)
                if op=='eq':
                    self._control('predicate',self.predicate & (1^(x^y)))
                elif not self.decided and x!=y:
                    self._control('predicate',int(x<y));self._control('decided',1)
                self.charges['boolean_gate_evaluations']+=3
                v=0
            if op not in ('eq','lt'):self._wb(ins.dst,i,v)
            if i==63:
                if op in ('add','sub') and self.carry:
                    self._control('phase','fault');return
                if op in ('eq','lt'):
                    for j in range(64):self._wb(ins.dst,j,self.predicate if j==0 else 0)
                self._next(ins)
            else:self._control('index',i+1)
        elif self.phase=='load':
            # Initial 64-tick address copy. The old address and partial copy
            # coexist and are both included in census at every cut.
            if self.index<64:
                self._wb(self.load_address,self.index,self._rb(ins.a,self.index))
                self._control('index',self.index+1)
                return
            if self.node is None:
                self._ref('node',self.store.root)
                self._control('page_depth',0)
                return
            if self.page_depth<self.store.height:
                bit=self._rb(self.load_address,12+self.store.height-1-self.page_depth)
                self.charges['bit_tests']+=1
                self.charges['directory_reference_read_bits']+=self.reference_width
                self.charges['navigation_links']+=1
                self._ref('node',self.node.right if bit else self.node.left)
                self._control('page_depth',self.page_depth+1)
                return
            # One actual bit selected from a fixed page. Selection has 15
            # finite tests; its address/bit inputs are read and charged.
            address=0
            for j in range(12):address |= self._rb(self.load_address,j)<<j
            self.charges['control_bit_reads']+=3
            self.charges['page_selection_tests']+=15
            value=self.store.bit_at_leaf(self.node,address,self.within_byte)
            self.charges['code_payload_bit_reads']+=1
            self._wb(ins.dst,self.load_bit,value)
            next_bit=self.load_bit+1
            self._control('load_bit',next_bit)
            if next_bit==self.load_width:
                # Unused high fields are binary zero, with actual writes.
                for j in range(self.load_width,64):self._wb(ins.dst,j,0)
                self._ref('node',None);self._next(ins)
            elif self.within_byte==7:
                self._control('within_byte',0);self._control('increment_bit',0)
                self._control('carry',1);self._control('phase','increment')
            else:
                self._control('within_byte',self.within_byte+1);self._ref('node',None)
        elif self.phase=='increment':
            i=self.increment_bit
            x=self._rb(self.load_address,i);c=self.carry
            self._wb(self.load_address,i,x^c)
            self._control('carry',x&c)
            if i==63:
                if self.carry:self._control('phase','fault');return
                self._ref('node',None);self._control('phase','load')
            else:self._control('increment_bit',i+1)
        elif self.phase=='emit':
            i=self.emit_bit
            bit=self._rb(self.opcode_register,i) if i<3 else self._rb(self.destination_register,i-3)
            self._control('wire_latch',bit);self._control('wire_occupied',1)
            self.charges['output_bits']+=1
            self.transcript.append(bit)  # Sink audit, never a runtime read.
            self._control('wire_occupied',0);self._control('emit_bit',i+1)
            if i==66:self._control('phase','end')
        elif self.phase=='end':
            self.charges['output_ends']+=1
            self._control('ended',1);self._control('phase','ack')
        elif self.phase=='ack':
            # The local ML-I acknowledgement is the input of this transition.
            self.charges['acknowledgements']+=1
            self._control('acknowledged',1);self._control('phase','return')
            self.output_acknowledgements+=1
        elif self.phase=='return':
            self.charges['returns']+=1
            self._ref('node',None);self._ref('selected',None)
            self._ref('saved_pc',None);self._control('phase','idle')

    def run(self):
        while not self.ready:self.tick()

    def result(self):
        # Independent observer view of retained result registers. This Python
        # tuple is not consumed by the executing component.
        if not self.ready or not self.acknowledged:
            raise RuntimeError('no acknowledged result')
        op=self.opcode_register.value
        return (None,'L','R','M')[op] if op<4 else ('J',self.destination_register.value)
