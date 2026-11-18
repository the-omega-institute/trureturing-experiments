"""Particular ML-F/ML-Q field graph, with literal primitive installation.

Graph construction is public preparation. All variable words, unary length
stacks, cache positions, receipt links and continuations are actual cells or
reference/bit registers. No service decision calls a host parser or law solver.
The original physical provider is a separately bound ML-I connection.
"""
from field_graph import Graph, loader, driver, growth
from field_io import conversions, output_decoder, print_buffer, literal
from field_decoder import lower, page_selector
from paid_code import decoder_program


class Service:
    def __init__(self,g):self.g=g
    def call(self,name):
        g=self.g;back=g.fresh(name+'.back')
        g.bind(name+'.return',back);g.j(name);g.label(back)
    def header(self,name):self.g.allocate(name)
    def copy(self,src,dst):
        g=self.g;self.header(dst);g.get('copy.cursor',src,1)
        loop,end=g.fresh(),g.fresh();g.label(loop);g.null('copy.cursor',end)
        g.read('copy.bit','copy.cursor');g.append(dst,'copy.bit')
        g.get('copy.cursor','copy.cursor',1);g.j(loop);g.label(end)
    def cat(self,src,dst):
        g=self.g;g.get('cat.cursor',src,1)
        loop,end=g.fresh(),g.fresh();g.label(loop);g.null('cat.cursor',end)
        g.read('cat.bit','cat.cursor');g.append(dst,'cat.bit')
        g.get('cat.cursor','cat.cursor',1);g.j(loop);g.label(end)
    def wrap(self,src,dst):
        g=self.g;g.get('wrap.cursor',src,1);g.put('wrap.one',1)
        loop,end=g.fresh(),g.fresh();g.label(loop);g.null('wrap.cursor',end)
        g.append(dst,'wrap.one');g.get('wrap.cursor','wrap.cursor',1);g.j(loop)
        g.label(end);literal(g,dst,'0');self.cat(src,dst)
    def end(self):
        g=self.g;ok=g.fresh();g.null('cursor',ok,'bad');g.label(ok)
    def take(self):self.call('take')
    def readv(self):self.call('readv')
    def readu(self):self.call('readu')
    def unary(self):self.call('unary')
    def shape(self):self.call('shape')
    def mask(self):self.call('mask')
    def cylinders(self):self.call('cylinders')
    def classify(self):self.call('classify')

    def routines(self):
        g=self.g
        g.label('take');g.null('cursor','take.missing')
        g.j('take.present')
        g.label('take.missing');g.test('strict','bad','need.more')
        g.label('take.present')
        g.read('taken','cursor');g.get('cursor','cursor',1);g.ret('take.return')
        # readv's marks stack and result chain survive every read/append cut.
        g.label('readv');self.header('vword');g.cp('vmarks','nil')
        g.label('v.unary');self.take();g.test('taken','v.mark','v.payload')
        g.label('v.mark');g.allocate('vtoken');g.set('vtoken',1,'vmarks')
        g.cp('vmarks','vtoken');g.j('v.unary')
        g.label('v.payload');g.null('vmarks','v.done')
        self.take();g.append('vword','taken');g.get('vmarks','vmarks',1);g.j('v.payload')
        g.label('v.done');g.ret('readv.return')
        # canonical U payload: one explicit zero, otherwise leading one.
        g.label('readu');self.readv();g.get('uhead','vword',1);g.null('uhead','bad')
        g.read('u.first','uhead');g.test('u.first','u.done')
        g.get('u.next','uhead',1);g.null('u.next','u.done','bad')
        g.label('u.done');g.ret('readu.return')
        # Unbounded unary value from an MSB-first canonical binary word. Doubling
        # scans a snapshot and appends one new token per old token.
        g.label('unary');self.header('unary.word');g.get('unary.cursor','vword',1)
        g.label('unary.bit');g.null('unary.cursor','unary.done')
        g.get('unary.snapshot','unary.word',2);g.get('unary.scan','unary.word',1)
        g.label('unary.double');g.null('unary.scan','unary.add')
        g.get('unary.next','unary.scan',1);g.put('unary.one',1)
        g.append('unary.word','unary.one')
        g.equal('unary.scan','unary.snapshot','unary.add')
        g.cp('unary.scan','unary.next');g.j('unary.double')
        g.label('unary.add');g.read('unary.value','unary.cursor')
        skip=g.fresh();g.test('unary.value',None,skip)
        g.put('unary.one',1);g.append('unary.word','unary.one');g.label(skip)
        g.get('unary.cursor','unary.cursor',1);g.j('unary.bit')
        g.label('unary.done');g.ret('unary.return')
        # shape constructs precisely 15N paired mask/value locations.
        g.label('shape');g.get('n.test','vword',1)
        for _ in range(3):g.null('n.test','bad');g.get('n.test','n.test',1)
        self.unary();self.header('shape.word');g.get('shape.n','unary.word',1)
        g.put('shape.zero',0)
        g.label('shape.loop');g.null('shape.n','shape.done')
        for _ in range(15):g.append('shape.word','shape.zero')
        g.get('shape.n','shape.n',1);g.j('shape.loop')
        g.label('shape.done');g.ret('shape.return')
        # prior stores the original numerator/denominator payloads and compares
        # them by length, then the first unequal MSB. No cancellation is done.
        g.label('prior');self.take();g.test('taken','bad')
        self.readu();g.cp('prior.a','vword');self.readu();g.cp('prior.b','vword')
        g.get('prior.bp','prior.b',1);g.read('prior.nonzero','prior.bp')
        g.test('prior.nonzero',None,'bad')
        g.get('prior.ap','prior.a',1);g.put('prior.decided',0);g.put('prior.greater',0)
        g.label('prior.compare');g.null('prior.ap','prior.aend');g.null('prior.bp','bad')
        g.read('prior.ab','prior.ap');g.read('prior.bb','prior.bp')
        cmpdone=g.fresh();g.test('prior.decided',cmpdone)
        g.boolean('prior.diff',('prior.ab','prior.bb'),lambda a,b:a^b)
        g.test('prior.diff',None,cmpdone)
        g.bitcopy('prior.greater','prior.ab');g.put('prior.decided',1)
        g.label(cmpdone);g.get('prior.ap','prior.ap',1);g.get('prior.bp','prior.bp',1)
        g.j('prior.compare')
        g.label('prior.aend');g.null('prior.bp','prior.same.length','prior.done')
        g.label('prior.same.length');g.test('prior.greater','bad')
        g.label('prior.done');g.ret('prior.return')
        # mask validates all 15N bits then the selected values. It uses a fresh
        # paired word so validation never overwrites the attained session cache.
        g.label('mask');self.header('mask.word');g.get('mask.shape','shape.word',1)
        g.label('mask.bits');g.null('mask.shape','mask.values.start')
        self.take();g.append('mask.word','taken')
        g.get('mask.shape','mask.shape',1);g.j('mask.bits')
        g.label('mask.values.start');g.get('mask.scan','mask.word',1)
        g.label('mask.values');g.null('mask.scan','mask.done')
        g.read('mask.selected','mask.scan');skip=g.fresh()
        g.test('mask.selected',None,skip);self.take();g.write('mask.scan',1,'taken')
        g.label(skip);g.get('mask.scan','mask.scan',1);g.j('mask.values')
        g.label('mask.done');g.ret('mask.return')
        g.label('cylinders');self.readu();self.unary();g.get('cyl.count','unary.word',1)
        g.label('cyl.loop');g.null('cyl.count','cyl.done')
        self.readv();g.cp('cyl.outer','cursor');g.get('cursor','vword',1)
        g.bitcopy('cyl.strict','strict');g.put('strict',1)
        self.take();g.test('taken','cyl.end')
        self.mask()
        g.label('cyl.end');g.null('cursor','cyl.next','bad')
        g.label('cyl.next');g.bitcopy('strict','cyl.strict');g.cp('cursor','cyl.outer');g.get('cyl.count','cyl.count',1);g.j('cyl.loop')
        g.label('cyl.done');g.ret('cylinders.return')
        g.label('arithmetic');self.readu();self.shape();self.call('prior');self.mask()
        self.cylinders();self.end();g.ret('arithmetic.return')
        self.classifier()

    def pair(self,a,b,slot=False):
        g=self.g
        for name,value in (('class.a',a),('class.b',b)):
            g.put(name+'0',(value>>1)&1);g.put(name+'1',value&1)
        if not slot:g.cp('class.slot','nil')
        g.ret('classify.return')
    def classifier(self):
        g=self.g;g.label('classify');g.get('class.cursor','operand',1)
        g.get('class.slot','cache',1)
        g.label('class.spine');g.null('class.slot','class.marker')
        g.null('class.cursor','class.branch')
        g.read('class.bit','class.cursor');g.get('class.cursor','class.cursor',1)
        g.test('class.bit','class.advance','class.suffix')
        g.label('class.advance');g.get('class.slot','class.slot',1);g.j('class.spine')
        g.label('class.branch');self.pair(2,2)
        g.label('class.suffix')
        table={'':(2,2),'0':(1,2),'1':(2,1),'00':(3,1),
               '01':(3,0),'10':(1,3),'11':(0,3)}
        marker={'':2,'0':2,'1':0,'00':2,'01':0,'000':1,'001':0}
        def trie(prefix,rows,depth):
            end=g.fresh();g.null('class.cursor',end)
            if depth:
                g.read('class.bit','class.cursor');g.get('class.cursor','class.cursor',1)
                one=g.fresh();g.test('class.bit',one)
                trie(prefix+'0',rows,depth-1);g.label(one);trie(prefix+'1',rows,depth-1)
            else:self.pair(3,3)
            g.label(end)
            value=rows.get(prefix,(3,3))
            if isinstance(value,int):value=(value,value)
            self.pair(*value,slot=value[0]!=value[1])
        trie('',table,2)
        g.label('class.marker');trie('',marker,3)

    def receive(self):
        g=self.g
        # Explicit bit/End channel. The external producer may fill only an empty
        # ready latch at a native waiting cut. Control is never replaced.
        g.label('receive');g.test('in.ready','receive.capture')
        g.emit('wait');g.j('receive')
        g.label('receive.capture');g.test('in.end','receive.seal')
        g.bitcopy('received.bit','in.bit');g.put('in.ready',0)
        g.append('request','received.bit');g.j('receive')
        g.label('receive.seal');g.put('strict',0);g.put('in.ready',0);g.get('cursor','request',1)
        g.test('initialized','command','initialize')
        g.label('need.more');g.j('receive')
        g.label('initialize');self.readu();self.shape();g.cp('cache','shape.word')
        self.call('prior');self.end()
        g.cp('public','request');g.put('initialized',1);g.j('entry')

    def commands(self):
        g=self.g;g.label('command');self.take();g.bitcopy('tag.first','taken')
        self.take();g.bitcopy('tag.second','taken')
        g.test('tag.first','command.high');g.test('tag.second','read','evaluate')
        g.label('command.high');g.test('tag.second','close','predict')
        g.label('evaluate');self.readv();g.cp('eval.input','vword')
        self.end();g.put('strict',1);g.get('cursor','eval.input',1)
        self.call('arithmetic');g.j('evaluate.begin')
        g.label('read');self.readv();g.cp('operand','vword');self.end()
        self.classify();self.header('port.transmitted')
        g.set('port.cell',1,'port.transmitted');g.put('port.zero',0)
        g.get('port.cursor','operand',1)
        g.label('port.send');g.null('port.cursor','port.call')
        g.read('port.bit','port.cursor');g.emit('emit','port.bit',1)
        g.read('port.occupied','port.cell',4);g.test('port.occupied',None,'fault')
        g.read('port.received','port.cell',3);g.append('port.transmitted','port.received')
        g.write('port.cell',4,'port.zero')
        g.get('port.cursor','port.cursor',1);g.j('port.send')
        g.label('port.call');g.emit('end',1)
        g.read('reply.0','port.cell',0);g.read('reply.1','port.cell',1)
        # Keep captured reply bits before acknowledging the provider's latch.
        g.emit('ack',1);g.j('receipt')
        g.label('receipt');self.header('receipt.record')
        self.wrap('operand','receipt.record')
        for b in ('reply.0','reply.1'):g.append('receipt.record',b)
        g.get('receipt.oldtail','receipts',2);g.null('receipt.oldtail','receipt.first')
        g.set('receipt.oldtail',3,'receipt.record');g.j('receipt.linked')
        g.label('receipt.first');g.set('receipts',1,'receipt.record')
        g.label('receipt.linked');g.set('receipts',2,'receipt.record')
        g.null('class.slot','cache.constant')
        # Actual returned code must be precisely one of the public table pair.
        g.boolean('reply.is.a',('reply.0','reply.1','class.a0','class.a1'),lambda a,b,c,d:(a==c)&(b==d))
        g.boolean('reply.is.b',('reply.0','reply.1','class.b0','class.b1'),lambda a,b,c,d:(a==c)&(b==d))
        g.test('reply.is.a','cache.zero');g.test('reply.is.b','cache.one','fault')
        g.label('cache.zero');g.put('cache.value',0);g.j('cache.update')
        g.label('cache.one');g.put('cache.value',1)
        g.label('cache.update');g.read('cache.known','class.slot',0);g.test('cache.known','cache.repeat')
        g.j('cache.write')
        g.label('cache.repeat');g.read('cache.old','class.slot',1)
        g.boolean('cache.agrees',('cache.old','cache.value'),lambda a,b:a==b)
        g.test('cache.agrees',None,'fault')
        g.label('cache.write');g.put('cache.one',1);g.write('class.slot',0,'cache.one')
        g.write('class.slot',1,'cache.value');g.j('read.reply')
        g.label('cache.constant')
        g.boolean('cache.agrees',('reply.0','reply.1','class.a0','class.a1'),lambda a,b,c,d:(a==c)&(b==d))
        g.test('cache.agrees',None,'fault')
        g.label('read.reply');self.header('answer');literal(g,'answer','10110')
        for b in ('reply.0','reply.1'):g.append('answer',b)
        g.j('print.answer')
        g.label('predict');self.take();g.bitcopy('predict.variant','taken')
        g.test('predict.variant','predict.cylinders')
        self.readv();g.cp('operand','vword');self.end();self.classify()
        self.header('pred.list');literal(g,'pred.list','1110100') # U(4)
        for code in range(4):
            self.header('pred.cylinder')
            g.boolean('pred.a',('class.a0','class.a1'),lambda a,b,c=code:2*a+b==c)
            g.boolean('pred.b',('class.b0','class.b1'),lambda a,b,c=code:2*a+b==c)
            present,empty,done=g.fresh(),g.fresh(),g.fresh()
            g.test('pred.a',present);g.test('pred.b',present,empty)
            g.label(empty);literal(g,'pred.cylinder','1');g.j(done)
            g.label(present);literal(g,'pred.cylinder','0')
            g.get('pred.scan','cache',1)
            loop,end=g.fresh(),g.fresh();g.label(loop);g.null('pred.scan',end)
            match,after=g.fresh(),g.fresh();g.equal('pred.scan','class.slot',match)
            literal(g,'pred.cylinder','0');g.j(after)
            g.label(match);literal(g,'pred.cylinder','1')
            g.label(after);g.get('pred.scan','pred.scan',1);g.j(loop)
            g.label(end);novalue=g.fresh();g.null('class.slot',novalue)
            # code b denotes one; code a denotes zero, in all informative pairs.
            g.append('pred.cylinder','pred.b');g.label(novalue)
            g.label(done);self.wrap('pred.cylinder','pred.list')
        g.j('predict.build')
        g.label('predict.cylinders')
        g.cp('pred.list.start','cursor');g.cp('shape.word','cache');self.cylinders()
        self.end();self.header('pred.list');g.cp('cat.cursor','pred.list.start')
        loop,end=g.fresh(),g.fresh();g.label(loop);g.null('cat.cursor',end)
        g.read('cat.bit','cat.cursor');g.append('pred.list','cat.bit')
        g.get('cat.cursor','cat.cursor',1);g.j(loop);g.label(end)
        g.label('predict.build');self.copy('public','eval.input')
        g.get('pred.scan','cache',1)
        loop,end=g.fresh(),g.fresh();g.label(loop);g.null('pred.scan',end)
        g.read('pred.bit','pred.scan');g.append('eval.input','pred.bit')
        g.get('pred.scan','pred.scan',1);g.j(loop);g.label(end)
        g.get('pred.scan','cache',1)
        loop,end=g.fresh(),g.fresh();g.label(loop);g.null('pred.scan',end)
        g.read('pred.bit','pred.scan');skip=g.fresh();g.test('pred.bit',None,skip)
        g.read('pred.value','pred.scan',1);g.append('eval.input','pred.value')
        g.label(skip);g.get('pred.scan','pred.scan',1);g.j(loop);g.label(end)
        self.cat('pred.list','eval.input');g.j('evaluate.begin')
        g.label('close');self.end();g.put('closed',1)
        self.header('answer');literal(g,'answer','110');g.j('print.answer')
        g.label('bad');self.header('answer');literal(g,'answer','000');g.j('print.answer')
        g.label('print.answer');print_buffer(g,'answer','service.return')
        g.label('service.return');g.test('closed','closed.wait','entry')
        g.label('closed.wait');g.emit('wait');g.j('closed.wait')
        g.label('entry');self.header('request');g.j('receive')


def service_graph(height=None,nregisters=361):
    """height=None is the frontend component graph; evaluation stops at its
    backed ready latch. Production height connects only the unchanged A path.
    A component stop supplies no numerical output or evaluator acceptance.
    """
    g=Graph();s=Service(g);g.label('start')
    g.put('initialized',0);g.put('closed',0);s.header('receipts');s.header('port.cell');g.j('entry')
    s.receive();s.commands();s.routines()
    g.label('evaluate.begin');g.put('eval.busy',1);g.get('input','eval.input',1)
    if height is None:
        g.put('eval.ready',1);g.emit('wait');g.j('evaluate.begin')
    else:
        graph=decoder_program();g.constant('bpc',0);g.j('decoder.'+str(graph.labels['boot']))
        lower(g,graph,height,'convert','image.fetched')
        g.label('image.fetched');g.ret('fetch.return')
        conversions(g,nregisters)
        g.label('converted.done');g.get('input','converted',1);g.j('load')
        loader(g);driver(g,image=True,halt_continuation='decode.output')
        # Autonomous boundary loop uses the same installed driver. An author
        # pause is not necessary for any B-instruction successor.
        ready=g.rows[g.labels['ready']];assert ready[0]=='wait'
        ready[:]=['goto',(), 'step',None]
        halt=g.rows[g.labels['halted']];assert halt[0]=='wait'
        halt[:]=['goto',(), 'decode.output',None]
        output_decoder(g)
        g.label('decoded.done');g.put('eval.busy',0);s.header('answer');literal(g,'answer','01')
        s.wrap('decoded','answer');g.j('print.answer')
        growth(g);page_selector(g)
    g.label('fault');g.emit('fault','invalid service/source phase')
    g.label('exhausted');g.emit('fault','finite supplied pool exhausted')
    return g.seal()
