"""Actual linked-word conversions and returned-head copier output scans.

No host constructor/conversion/decoder runs on the executing path. Accepted
arithmetic framing remains A's premise; the service parser is not supplied here.
"""
from field_graph import Graph, Engine, loader, driver, growth


def literal(g,header,word):
    for b in word:g.put('io.bit',int(b));g.append(header,'io.bit','io.cell')


def nwmg_bit(g,header,bit):
    g.put('io.one',1);g.append(header,'io.one','io.cell');g.append(header,bit,'io.cell')


def copier_literal(g,word):
    for b in word:g.put('io.bit',int(b));nwmg_bit(g,'converted','io.bit')


def call_sourcebit(g,bit):
    back=g.fresh('sourcebit.return')
    g.bitcopy('source.bit',bit);g.bind('source.return',back)
    g.j('convert.sourcebit');g.label(back)


def conversions(g,nregisters=361,done='converted.done'):
    g.label('convert')
    g.cp('received','input')
    g.allocate('converted');g.put('source.first',1)
    copier_literal(g,'10010000')
    g.put('source.literal',1)
    call_sourcebit(g,'source.literal');call_sourcebit(g,'source.literal')
    g.label('convert.input.loop')
    g.null('input','convert.separators')
    g.read('source.input','input');g.get('input','input',1)
    # Four-symbol bit0=01, bit1=10, with actual data-dependent control.
    one,after=g.fresh(),g.fresh()
    g.test('source.input',one)
    g.put('source.literal',0);call_sourcebit(g,'source.literal')
    g.put('source.literal',1);call_sourcebit(g,'source.literal');g.j(after)
    g.label(one)
    g.put('source.literal',1);call_sourcebit(g,'source.literal')
    g.put('source.literal',0);call_sourcebit(g,'source.literal')
    g.label(after);g.j('convert.input.loop')
    g.label('convert.separators')
    for _ in range(nregisters):
        g.put('source.literal',1)
        call_sourcebit(g,'source.literal');call_sourcebit(g,'source.literal')
    copier_literal(g,'00100000');g.j(done)
    g.label('convert.sourcebit')
    g.put('source.zero',0);g.put('source.one',1)
    for b in ('source.zero','source.one','source.zero','source.first'):
        nwmg_bit(g,'converted',b)
    g.boolean('source.inverse',('source.bit',),lambda b:1-b)
    for b in ('source.inverse','source.bit','source.zero','source.zero'):
        nwmg_bit(g,'converted',b)
    g.put('source.first',0);g.ret('source.return')


def conversion_graph(nregisters=361):
    g=Graph();g.label('start');g.j('convert')
    conversions(g,nregisters)
    g.label('converted.done');g.emit('wait');g.j('converted.done')
    g.label('exhausted');g.emit('fault','finite supplied pool exhausted')
    return g.seal()


def read_copier_block(g):
    # scan is an actual plus-arm B cell, at the first bit of an NWMG pair.
    # The two link steps preserve all intermediate references and bit latches.
    for i in range(8):
        first=g.fresh('pair.first')
        g.get('decode.bitref','scan',2);g.null('decode.bitref','fault')
        g.get('decode.mark','decode.bitref',1);g.null('decode.mark','fault')
        g.get('scan','scan',1)
        g.get('decode.bitref','scan',2);g.null('decode.bitref','fault')
        g.get('decode.mark','decode.bitref',1)
        one,done=g.fresh(),g.fresh();g.null('decode.mark',done)
        g.put(f'block.{i}',1);g.j(one)
        g.label(done);g.put(f'block.{i}',0)
        g.label(one);g.get('scan','scan',1)


def validate_bits(g,checks):
    for bit,value in checks:
        if value:g.test(bit,None,'fault')
        else:g.test(bit,'fault')


def data_block(g):
    read_copier_block(g)
    validate_bits(g,[(f'block.{i}',v) for i,v in ((0,0),(1,1),(2,0),(6,0),(7,0))])
    g.boolean('decode.valid',('block.4','block.5'),lambda a,b:a^b)
    g.test('decode.valid',None,'fault')


def output_decoder(g,done='decoded.done'):
    # It follows the actual returned active header, then the actual binary-head
    # flag. No retired-copy count or computed integer coordinate selects data.
    g.label('decode.output');g.allocate('decoded')
    g.cp('scan','head');read_copier_block(g)
    validate_bits(g,[(f'block.{i}',v) for i,v in ((0,1),(1,0),(2,0),(3,1),(4,0))])
    g.label('decode.find.head');data_block(g)
    g.test('block.3','decode.first');g.j('decode.find.head')
    g.label('decode.first');g.bitcopy('decode.firstbit','block.5')
    g.label('decode.second');data_block(g)
    # 01=bit0,10=bit1,11=the actual separator;00 is not output.
    g.test('decode.firstbit','decode.first.one')
    g.test('block.5','decode.zero','fault')
    g.label('decode.first.one');g.test('block.5','decode.separator')
    g.put('decode.value',1);g.j('decode.append')
    g.label('decode.zero');g.put('decode.value',0)
    g.label('decode.append');g.append('decoded','decode.value','decoded.cell')
    data_block(g);g.bitcopy('decode.firstbit','block.5');g.j('decode.second')
    g.label('decode.separator');g.j(done)


def print_buffer(g,header,done='printed'):
    g.get('print.cursor',header,1)
    loop=g.fresh('print.loop');end=g.fresh('print.end')
    g.label(loop);g.null('print.cursor',end)
    g.read('print.latch','print.cursor');g.emit('emit','print.latch')
    g.get('print.cursor','print.cursor',1);g.j(loop)
    g.label(end);g.emit('end');g.emit('ack');g.emit('return');g.j(done)


def output_graph():
    # Start at a separately installed comb fixture; this graph is a component
    # author check and not a full evaluator execution certificate.
    g=Graph();g.label('start');g.j('load')
    loader(g,ready='decode.output');output_decoder(g)
    g.label('decoded.done');print_buffer(g,'decoded')
    g.label('printed');g.emit('wait');g.j('printed')
    growth(g)
    g.label('fault');g.emit('fault','invalid actual returned record or pair')
    g.label('exhausted');g.emit('fault','finite supplied pool exhausted')
    return g.seal()


def observe_buffer(header):
    # Independent view for checks; neither converter nor print consumes it.
    p=header.refs[1];out=[]
    while p is not None:out.append(str(p.bits&1));p=p.refs[1]
    return ''.join(out)


def install_input(e,bits):
    # Input data installation is still host preparation, not Receive/End/Ack.
    p=e.registers['pool'];first=last=None
    for bit in bits:
        assert bit in '01' and p is not None
        c,p=p,p.refs[0];c.refs[0]=None;c.bits=128+int(bit)
        if last is None:first=c
        else:last.refs[1]=c
        last=c
    e.registers['pool']=p;e.registers['input']=first
