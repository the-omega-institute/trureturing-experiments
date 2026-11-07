"""Literal connection selection for the particular installed service.

Normal output and original-port bit output have distinct literal successors.
The source provider remains exactly session.OriginalPort.call. Connection state
is separate from service state; its finite received word and actual reply latch
are backed by the existing paid Buffer arena. It supplies no parser, cache,
classification, evaluator answer, instruction or service continuation.
"""
import native_engine as ne
from lab import Arena, Buffer
from field_service import service_graph


class ServiceCore(ne.Core):
    def __init__(self,r,b,endpoint):
        super().__init__(r,b);self.endpoint=endpoint
    def node(self,op,a=None,imm=0,zero=None,one=None):
        if op in ('send','end','ack'):
            ordinary=super().node(op,a,imm,zero,one)
            special=super().node({'send':'portsend','end':'sourcecall','ack':'sourceack'}[op],
                                 self.endpoint,0,zero,one)
            if op=='send':
                ordinary=self.seq([self.br('bx')],ordinary)
                special=self.seq([self.br('bx')],special)
            refused=super().node('fault')
            return self.choice(tuple(f'imm.{i}' for i in reversed(range(16))),
                               {0:ordinary,1:special},refused)
        return super().node(op,a,imm,zero,one)


class OriginalConnection:
    """The original finite-word call and its supplied physical latch operations.

    The installed graph assembles the transmitted native chain. invoke's word
    projection is its denotation at the supplied OriginalPort.call boundary,
    not a service parser/selector/continuation. Only that unchanged call obtains
    source information. Its own private provider operations stay unchanged.
    """
    def __init__(self,port,*,storage=None):
        self.port=port;self.a=Arena(registers=8,storage=storage);self.calls=0
        self.binding=self.a.alloc('service_original_connection')
        for key,value in (('reply',None),('source',port.root)):
            self.a.write(self.binding,key,value)
    def send(self,bit,target):
        assert bit in (0,1) and not target.bit(4)
        # A physical bit input transfers into the actual occupied peer latch.
        # Every append and occupancy release is a later installed field action.
        target.put(3,bit);target.put(4,1)
        self.a.charges['connection_operand_bits']+=1
    def invoke(self,target):
        assert not target.bit(2) and not target.bit(4) and self.port.grant
        # This projects the actual word that the native peer received. No host
        # value chooses a service continuation or replaces a transmitted bit.
        p=target.refs[1].refs[1];value=[]
        while p is not None:
            value.append('R' if p.bit(0) else 'L');p=p.refs[1]
        returned=self.port.call(''.join(value),self.a)
        # The supplied original operation returns its actual two-bit latch.
        first=self.a.read(returned.header,'head')
        a=self.a.read(first,'bit');second=self.a.read(first,'next')
        b=self.a.read(second,'bit')
        target.put(0,a);target.put(1,b);target.put(2,1)
        self.a.write(self.binding,'reply',returned.header);self.calls+=1
        self.a.charges['connection_reply_latch_bits']+=2
    def ack(self,target):
        assert target.bit(2)
        target.put(2,0)
        # One latch acknowledgement; no bulk clear/allocation is hidden here.
        # Provider-return and transmitted native backing remain retained/paid.
        self.a.charges['connection_acknowledgements']+=1


class ServiceEngine(ne.NativeEngine):
    native_operations=ne.NOPS+('portsend','sourcecall','sourceack')
    def __init__(self,*args,connection=None,**kwargs):
        self.connection=connection
        super().__init__(*args,**kwargs)
        for n,c in zip(self.program.nodes,self.code):
            if n.op=='emit' and len(n.args)==2:c.immediate=n.args[1]
            elif n.op in ('end','ack') and n.args:c.immediate=n.args[0]
    def make_core(self,r,b,registers):
        return ServiceCore(r,b,registers['port.cell']).build()
    def microtick(self):
        if not self.sealed:self.seal()
        n=self.native_pc;op=self.native_operations[n.bits&31]
        if op not in ('portsend','sourcecall','sourceack'):return super().microtick()
        assert self.connection is not None
        self.native_steps+=1;self.microsteps+=1;self.native_events[op]+=1
        try:
            if op=='portsend':self.connection.send(self.read_bit,n.refs[0].refs[0]);self._event('port_operand_bits')
            elif op=='sourcecall':
                self.connection.invoke(n.refs[0].refs[0]);self._event('original_calls')
            else:self.connection.ack(n.refs[0].refs[0]);self._event('original_acknowledgements')
        except (AssertionError,RuntimeError):
            self._event('native_controller_selections');self._event('native_refusals')
            raise
        for kind in ('native_controller_selections','native_successor_reference_reads','native_control_reference_writes'):
            self._event(kind)
        self.native_pc=n.refs[2]
    def deliver(self,bit=None,*,end=False):
        assert self.sealed and self.waiting and not self.bits['in.ready']
        assert not self.bits['closed'] and not self.bits['eval.busy']
        assert end or bit in (0,1)
        self.bw('in.bit',int(bit or 0));self.bw('in.end',int(end));self.bw('in.ready',1)
        self.resume();self.run()
    def feed(self,bits,*,seal=True):
        for b in bits:self.deliver(int(b))
        if seal:self.deliver(end=True)
    def service_census(self):
        # This is an actual common retained cut, not the sum of fixture peaks.
        out={'installation':self.census()}
        if self.connection:
            out['connection']=self.connection.a.census()
            out['same_original_source']=self.connection.port.source.a.census()
        return out
