"""Source-bound bounded consumers and mutation checks; no full C execution."""
from __future__ import annotations
import argparse
from collections import Counter
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

from packed_constructor import PackedStorage, PackedBlock, PackedSource, PackedObj
from native_install import Block, Source, Installer, BitKernel, LAYOUTS
from native_engine import NativeEngine
from service_native import ServiceEngine, OriginalConnection
from field_graph import Graph, b_graph, install_b, observe_b
from field_io import observe_buffer, output_graph
from field_service import service_graph
from lab import Arena, Buffer, reference_run
from session import prepare_source, OriginalPort
from service_checks import cache, receipts, refused
from framing import u, v, rational, parse_reply, REPLIES, words, cylinder_word

HERE=Path(__file__).resolve().parent
assert sys.version_info >= (3,12) and __debug__

def sha(data):return hashlib.sha256(data).hexdigest()
def atomic(path, obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');tmp.replace(path)

def bind(e):
    ids={};n=0
    for s in e.installation.sources:
        for x in (s,s.block,s.forward):
            if x is not None and x not in ids:ids[x]=n;n+=1
    return ids

def state(e, ids, deep=False):
    ref=lambda x:None if x is None else ids[x]
    out={'pc':ref(e.native_pc),'refs':{k:ref(x.refs[0]) for k,x in e.cr.items()},
         'bits':{k:x.bit(0) for k,x in e.cb.items()},
         'latches':[ref(e.read_ref),ref(e.saved_ref),ref(e.address),e.read_bit],
         'ticks':e.ticks,'microsteps':e.microsteps,'native_steps':e.native_steps,
         'charges':dict(e.charges),'events':dict(e.native_events),'output':e.output[:],
         'census':e.census(),'installation_charges':dict(e.installation.charges)}
    if deep:
        h=hashlib.sha256()
        for s in e.installation.sources:
            for b in (s.block,s.forward):
                row=(len(b.refs),b.width,b.grant,[ref(x) for x in b.refs],ref(b.next))
                h.update(json.dumps(row,separators=(',',':')).encode());h.update(b.data)
            h.update(json.dumps((ref(s.next),ref(s.forward),s.grant)).encode())
        out['all_block_fields_sha256']=h.hexdigest()
    return out

def compare(a,b,ia,ib,deep=False):
    x,y=state(a,ia,deep),state(b,ib,deep)
    assert x==y,[(k,x[k],y[k]) for k in x if x[k]!=y[k]][:2]
    return sha(json.dumps(x,sort_keys=True).encode())

def installer(k):
    sources=[]
    for i,(r,b) in enumerate(LAYOUTS):
        s=Source(Block([None]*r,bytearray((b+7)//8),b));s.block.grant=1
        for j in range(b):s.block.put(j,(i+j)&1)
        sources.append(s)
    for i,s in enumerate(sources):
        s.next=sources[i+1] if i+1<len(sources) else None
        for j in range(len(s.block.refs)):s.block.refs[j]=sources[(i+j)%len(sources)]
    storage=PackedStorage();ps,mapping=storage.import_sources(sources)
    a=Installer(sources,k.independent_header());b=Installer(ps,k.independent_header(),storage=storage)
    phases=Counter();cuts=0
    ai={x:i for i,s in enumerate(sources) for x in (s,)}
    bi={x:i for i,s in enumerate(ps) for x in (s,)}
    seeds={id(ins):ins.seed_census() for ins in (a,b)}
    def view(ins,ids):
        if ins.pending is not None:ids.setdefault(ins.pending,len(ids))
        seed=seeds[id(ins)]
        c={key:ins.base[key]+ins.accounted[key]+seed[key] for key in ('objects','refs','bits')}
        c['total_bits']=c['bits']+ins.width*c['refs']+c['objects']
        return {'phase':ins.phase,'done':ins.done,'steps':ins.steps,'charges':dict(ins.charges),
            'cursor':None if ins.cursor is None else ids[ins.cursor],
            'census':c,'forward_refs':[(None if x.forward is None else ids[x.forward], [] if x.forward is None else [None if z is None else ids[z] for z in x.forward.refs]) for x in ins.sources],'pending':None if ins.pending is None else (ins.pending.width,ins.pending.grant,bytes(ins.pending.data)),
            'bit_state':(ins.kernel.state.latch,ins.kernel.state.done,tuple(ins.kernel.state.counts))}
    while not a.done:
        phase=a.phase;a.tick(1);b.tick(1);phases[phase]+=1;cuts+=1
        assert view(a,ai)==view(b,bi),phase
    a.verify();b.verify()
    assert all(isinstance(s,PackedSource) and isinstance(s.block,PackedBlock) and isinstance(s.forward,PackedBlock) for s in ps)
    # An actual paid pending target has a zero grant until the next source step.
    assert phases['allocate']==phases['grant']==len(LAYOUTS)
    saved=ps[0].forward.refs[1];ps[0].forward.refs[1]=ps[0].forward
    assert ps[0].forward.refs[1] is not saved
    refused(b.verify);ps[0].forward.refs[1]=saved;b.verify()
    return dict(cuts=cuts,phases=dict(phases),steps=a.steps,modeled_census=a.census(),storage=storage.metrics())

def primitive(k):
    g=Graph();g.label('start');g.allocate('x');g.allocate('y');g.put('v',1)
    for i in range(4):
        g.set('x',i,'y');g.get('z','x',i);g.equal('z','y','ok'+str(i),'fault');g.label('ok'+str(i))
    for i in range(8):g.write('x',i,'v');g.read('b','x',i);g.test('b',None,'fault')
    g.bind('return','finish');g.ret('return');g.label('finish')
    g.emit('emit','v');g.emit('end');g.emit('ack');g.emit('return');g.emit('wait')
    g.label('fault');g.emit('fault');g.label('exhausted');g.emit('fault');p=g.seal()
    storage=PackedStorage();a=NativeEngine(p,3,kernel=k);b=NativeEngine(p,3,kernel=k,storage=storage)
    a.seal();b.seal();ia,ib=bind(a),bind(b);cuts=0;pending=Counter()
    while not a.waiting:
        a.microtick();b.microtick();cuts+=1
        # Every internal primitive cut compares all latches and retained census.
        # Full graph hashing is reserved for commit/terminal cuts below.
        assert a.native_pc.bits==b.native_pc.bits
        for key in ('read_ref','saved_ref','address'):
            x,y=getattr(a,key),getattr(b,key)
            assert (None if x is None else ia[x])==(None if y is None else ib[y])
        assert a.read_bit==b.read_bit and a.output==b.output and a.charges==b.charges
        if b.cb['occupied'].bit(0):pending['occupied']+=1
        if b.cb['ended'].bit(0) and not b.cb['acknowledged'].bit(0):pending['End-before-Ack']+=1
        if b.cb['acknowledged'].bit(0):pending['Ack-before-Return']+=1
    terminal=compare(a,b,ia,ib,True)
    assert a.output==[1] and all(pending.values())
    mutations=[]
    def detect(name,apply,undo):
        apply();detected=False
        try:compare(a,b,ia,ib,True)
        except (AssertionError,KeyError):detected=True
        finally:undo()
        assert detected,name
        mutations.append(name)
    cell=b.cells[0];v0=cell.refs[0];bits0=cell.bits;grant=cell.grant
    detect('reference',lambda:cell.refs.__setitem__(0,cell),lambda:cell.refs.__setitem__(0,v0))
    detect('value',lambda:setattr(cell,'bits',bits0^1),lambda:setattr(cell,'bits',bits0))
    detect('grant',lambda:setattr(cell,'grant',0),lambda:setattr(cell,'grant',grant))
    pc=b.native_pc
    detect('dispatch',lambda:setattr(b,'native_pc',pc.refs[2]),lambda:setattr(b,'native_pc',pc))
    old=b.read_bit
    detect('pending',lambda:setattr(b,'read_bit',old^1),lambda:setattr(b,'read_bit',old))
    detect('output',lambda:b.output.append(0),lambda:b.output.pop())
    for key in ('code_bit_reads','native_control_reference_writes'):
        detect(key,lambda key=key:b.charges.update({key:-1}),lambda key=key:b.charges.update({key:1}))
    assert storage.metrics()['reference_reads']>cuts and all(isinstance(c,PackedBlock) for c in b.cells)
    return dict(internal_cuts=cuts,pending=dict(pending),terminal_sha256=terminal,mutations=mutations,storage=storage.metrics())

def b_prefixes(k):
    rows=[];p=b_graph()
    cases=[(('L','M','R','R','M','L','L','L'),'',8),((('J',2),'M','R'),'',4),(('M',('J',1)),'1',5),((),'01',1)]
    for w,x,bound in cases:
        storage=PackedStorage()
        a,ac=install_b(w,x,512,p,lambda *a,**kw:NativeEngine(*a,kernel=k,**kw))
        b,bc=install_b(w,x,512,p,lambda *a,**kw:NativeEngine(*a,kernel=k,storage=storage,**kw))
        a.seal();b.seal();ac=tuple(a.forward_observer[c] for c in ac);bc=tuple(b.forward_observer[c] for c in bc)
        ia,ib=bind(a),bind(b);a.run();b.run();states=[observe_b(a,ac)]
        compare(a,b,ia,ib)
        for _ in range(bound):
            a.resume();b.resume();a.run();b.run();compare(a,b,ia,ib)
            if a.last is p.entries['halted']:break
            states.append(observe_b(a,ac));assert states[-1]==observe_b(b,bc)
        assert states==reference_run(w,x,bound)[0]
        deep=compare(a,b,ia,ib,True)
        rows.append(dict(program=w,input=x,prefix=states,native_steps=a.native_steps,terminal_sha256=deep,storage=storage.metrics()))
        del a,b;gc.collect()
    return rows

def arena_checks():
    storage=PackedStorage();a=Arena();b=Arena(storage=storage)
    x=Buffer(a);y=Buffer(b)
    for bit in '010101':
        x.append(int(bit));y.append(int(bit));assert a.census()==b.census() and a.charges==b.charges
    assert x.bits()==y.bits();x.clear();y.clear()
    assert a.census()==b.census() and a.charges==b.charges
    assert all(isinstance(o,PackedObj) for o in b.objects)
    assert any(o.fields.get('bit','missing') is None for o in b.objects)
    old=b.charges['reference_copy_bits'];b.charges['reference_copy_bits']-=1
    assert a.charges!=b.charges;b.charges['reference_copy_bits']=old
    return dict(census=b.census(),widening_charge_falsifier=True,cleanup_nullable_bit=True,storage=storage.metrics())

def session_case(k,n,bits,operands,prior):
    storage=PackedStorage();source=prepare_source(bits);psource=prepare_source(bits,storage=storage)
    assert source.logical()==psource.logical() and source.routing()==psource.routing()
    port=OriginalPort(source);pport=OriginalPort(psource)
    a=ServiceEngine(service_graph(),120000,kernel=k,connection=OriginalConnection(port))
    b=ServiceEngine(a.program,120000,kernel=k,connection=OriginalConnection(pport,storage=storage),storage=storage)
    a.seal();b.seal();ia,ib=bind(a),bind(b);a.run();b.run()
    a.feed(u(n)+prior);b.feed(u(n)+prior);compare(a,b,ia,ib)
    digests=[]
    for word in operands:
        command='01'+v(word.translate(str.maketrans('LR','01')))
        a.feed(command);b.feed(command)
        assert port.calls==pport.calls and receipts(a)==receipts(b) and cache(a)==cache(b)
        assert source.a.charges==psource.a.charges and a.connection.a.charges==b.connection.a.charges
        assert a.service_census()==b.service_census()
        digests.append(compare(a,b,ia,ib))
    assert set(x[1] for x in port.calls)>=set(REPLIES)
    # Original indices remain distinct even for equal observations, and repeats
    # retain separate receipts. A physical indexed-cache mutation is observable.
    cp=b.registers['cache'].refs[1];old=cp.bit(0);cp.put(0,old^1)
    assert cache(a)!=cache(b);cp.put(0,old)
    header=b.registers['cache'];head=header.refs[1];header.refs[1]=head.refs[1]
    assert cache(a)!=cache(b);header.refs[1]=head
    rec=b.registers['receipts'].refs[1].refs[2];old=rec.bit(0);rec.put(0,old^1)
    assert receipts(a)!=receipts(b)
    rec.put(0,old)
    before=len(port.calls)
    cylinders=[None,tuple(() for _ in range(n)),tuple(((0,1),) if j==0 else () for j in range(n))]
    command='101'+words([cylinder_word(c) for c in cylinders])
    a.feed(command);b.feed(command)
    assert a.bits['eval.ready'] and b.bits['eval.ready'] and a.bits['eval.busy']
    assert observe_buffer(a.registers['eval.input'])==observe_buffer(b.registers['eval.input'])
    assert len(port.calls)==len(pport.calls)==before
    refused(lambda:a.deliver(0));refused(lambda:b.deliver(0))
    terminal=compare(a,b,ia,ib,True)
    return dict(N=n,prior=prior,calls=port.calls,receipts=receipts(b),cache=cache(b),cut_hashes=digests,
        native_steps=b.native_steps,terminal_sha256=terminal,common_census=b.service_census(),storage=storage.metrics(),
        pending_evaluator=True,arithmetic_result=False,index_and_receipt_falsifiers=True)

def sessions(k):
    ops=['','LL','LRR','LLR','LL','LR','LLL','LRL','R'*15+'LL','R'*45+'LL','L'*130,'R'*145]
    n4=session_case(k,4,[0,1]*30,ops,'0'+u(2)+u(6))
    gc.collect()
    # Authentic outside-support source; no arithmetic evaluation is performed.
    bits=[int(c) for j in range(9) for c in format(j,'07b')+'0'*8]
    n9=session_case(k,9,bits,ops+['R'*j+'LL' for j in range(15*9)],'0'+u(0)+u(2))
    return [n4,n9]

def sessions_one(k,n):
    ops=['','LL','LRR','LLR','LL','LR','LLL','LRL','R'*15+'LL','R'*45+'LL','L'*130,'R'*145]
    if n==4:return session_case(k,4,[0,1]*30,ops,'0'+u(2)+u(6))
    bits=[int(c) for j in range(9) for c in format(j,'07b')+'0'*8]
    return session_case(k,9,bits,ops+['R'*j+'LL' for j in range(15*9)],'0'+u(0)+u(2))


def record_census():
    from packed_constructor import PackedState, PackedRefs
    state=PackedState(1);h=state.allocate(0,(4,8));state.grant(h,0)
    assert state.modeled_census()[0]['grant_bits']==1
    assert state.modeled_census()[0]['bits']==4+8+1
    values=[None,(0,0),(3,2**64),(2**81,2**129),(2**1024,2**1023)]
    v=PackedRefs(len(values))
    for i,x in enumerate(values):v[i]=x;assert list(v)[:i+1]==values[:i+1]
    return dict(pending_zero_grant_counted=True,arbitrary_integer_codec_cases=len(values),host_width_bytes=v.width)

def decoded_output(k):
    from field_graph import loader,growth
    from field_service import Service
    from field_io import output_decoder,print_buffer,literal
    from copier import initial_word
    from nwmg import encode
    g=Graph();service=Service(g);g.label('start');service.header('port.cell');g.j('load')
    loader(g,ready='decode.output');output_decoder(g)
    g.label('decoded.done');service.header('answer');literal(g,'answer','01');service.wrap('decoded','answer')
    print_buffer(g,'answer','finished');g.label('finished');g.emit('wait');g.j('finished')
    growth(g);g.label('fault');g.emit('fault');g.label('exhausted');g.emit('fault')
    p=g.seal();word='01';bits=encode(initial_word(''.join('01' if b=='0' else '10' for b in word)+'11'))
    storage=PackedStorage()
    a,_=install_b((),bits,24000,p,lambda *a,**kw:ServiceEngine(*a,kernel=k,**kw))
    b,_=install_b((),bits,24000,p,lambda *a,**kw:ServiceEngine(*a,kernel=k,storage=storage,**kw))
    a.seal();b.seal();ia,ib=bind(a),bind(b);a.run();b.run()
    terminal=compare(a,b,ia,ib,True)
    actual=''.join(map(str,b.output));assert parse_reply(actual)==(1,word)
    assert observe_buffer(b.registers['decoded'])==word
    assert b.charges['output_ends']==b.charges['acknowledgements']==b.charges['returns']==1
    b.charges['bit_field_reads']-=1
    refused(lambda:compare(a,b,ia,ib));b.charges['bit_field_reads']+=1
    return dict(actual_word=word,actual_frame=actual,terminal_sha256=terminal,native_steps=b.native_steps,
        decoder_charge_falsifier=True,storage=storage.metrics(),
        scope='actual returned-head/Arithmetic-framing component; supplied arbitrary record, no evaluator law result')

def image_decoder(k):
    import gzip,struct,tempfile
    from constructor_image import build_image,Image
    from paid_code import CodeStore
    from field_decoder import decoder_graph,request
    with tempfile.TemporaryDirectory(prefix='rro-packed-image-') as td:
        td=Path(td);raw=b'RROTABLE'+struct.pack('<II',2,1)+struct.pack('<BbI',0,-1,0)+struct.pack('<BbI',1,1,1)
        (td/'binary').write_bytes(raw);build_image(td/'binary',td)
        data=gzip.decompress((td/'W_U.block.gz').read_bytes());(td/'raw').write_bytes(data)
        oracle=Image(td/'raw');store=CodeStore(data);p=decoder_graph(store.height);storage=PackedStorage()
        a=NativeEngine(p,1,bindings={'image.root':store.root},store=store,kernel=k)
        b=NativeEngine(p,1,bindings={'image.root':store.root},store=store,kernel=k,storage=storage)
        a.seal();b.seal();ia,ib=bind(a),bind(b);a.run();b.run();rows=[]
        for pc in (0,7,13*oracle.nm,13*oracle.nm+1):
            x=request(a,pc);y=request(b,pc);assert x==y
            op,dest=y;value=(None,'L','R','M')[op] if op<4 else ('J',dest)
            assert value==oracle.fetch(pc)
            rows.append(dict(pc=pc,decoded=value,state_sha256=compare(a,b,ia,ib)))
        invalid=13*oracle.nm+2
        refused(lambda:request(a,invalid));refused(lambda:request(b,invalid))
        deep=compare(a,b,ia,ib,True)
        return dict(cases=rows,invalid=invalid,terminal_sha256=deep,native_steps=b.native_steps,
            image_sha256=sha(data),image_bytes=len(data),storage=storage.metrics(),
            scope='existing installed decoder consuming source-generated copier/B image; small table, no full C')


def guards(k):
    from framing import reply
    invalid=[]
    for bad in ('ack','return'):
        g=Graph();g.label('start')
        if bad=='return':g.emit('end')
        g.emit(bad);g.emit('wait');p=g.seal()
        a=NativeEngine(p,1,kernel=k)
        b=NativeEngine(p,1,kernel=k,storage=PackedStorage())
        a.seal();b.seal();ia,ib=bind(a),bind(b)
        refused(a.run);refused(b.run)
        invalid.append(dict(guard=bad,cut_sha256=compare(a,b,ia,ib,True)))
        del a,b;gc.collect()
    storage=PackedStorage();source=prepare_source([0,1]*30)
    psource=prepare_source([0,1]*30,storage=storage)
    port=OriginalPort(source);pport=OriginalPort(psource)
    a=ServiceEngine(service_graph(),20000,kernel=k,connection=OriginalConnection(port))
    b=ServiceEngine(a.program,20000,kernel=k,connection=OriginalConnection(pport,storage=storage),storage=storage)
    a.seal();b.seal();ia,ib=bind(a),bind(b);a.run();b.run();cuts=[]
    # Endpoint one, with an unreduced original rational representation.
    a.feed(u(4)+'0'+u(3)+u(3));b.feed(u(4)+'0'+u(3)+u(3))
    cuts.append(compare(a,b,ia,ib))
    command='01'+v('00')
    a.feed(command[:-1],seal=False);b.feed(command[:-1],seal=False)
    assert not port.calls and not pport.calls and not a.output and not b.output
    a.deliver(end=True);b.deliver(end=True);cuts.append(compare(a,b,ia,ib))
    assert not port.calls and not pport.calls and not a.output and not b.output
    a.feed(command[-1:]);b.feed(command[-1:]);cuts.append(compare(a,b,ia,ib))
    assert receipts(a)==receipts(b) and len(port.calls)==len(pport.calls)==1
    badframes=('110','01'+v('')+'1','00'+v(''),'00'+v(u(3)+rational(0)))
    for word in badframes:
        before=len(a.output);calls=port.calls[:]
        a.feed(word);b.feed(word);cuts.append(compare(a,b,ia,ib))
        assert ''.join(map(str,a.output[before:]))==reply(0)
        assert port.calls==pport.calls==calls
    before=len(a.output)
    a.feed('11');b.feed('11');cuts.append(compare(a,b,ia,ib,True))
    assert ''.join(map(str,b.output[before:]))==reply(3) and a.bits['closed'] and b.bits['closed']
    refused(lambda:a.feed('01'+v('00')));refused(lambda:b.feed('01'+v('00')))
    assert a.service_census()==b.service_census()
    assert a.connection.a.charges==b.connection.a.charges and source.a.charges==psource.a.charges
    assert a.native_steps==b.native_steps==sum(b.native_events.values())==b.charges['native_controller_selections']
    cuts.append(compare(a,b,ia,ib,True))
    return dict(invalid_primitive_order=invalid,cut_hashes=cuts,malformed_frames=len(badframes),
        endpoint_unreduced_prior='3/3',partial_received_read=True,original_calls=port.calls,
        actual_output=b.output,permanent_close_guard=True,native_steps=b.native_steps,
        common_census=b.service_census(),storage=storage.metrics(),arithmetic_result=False,
        scope='actual paired installed service, source and original port; framing/Close and primitive I/O guards, no evaluation')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--theory-root',type=Path,required=True)
    ap.add_argument('--groups',default='arena,installer,primitive,b,session')
    args=ap.parse_args();started=time.perf_counter();out={};k=BitKernel()
    from pins import source_manifest
    source_manifest(args.theory_root,HERE/'theory-candidate.json')
    try:
        for name in args.groups.split(','):
            t=time.perf_counter()
            out[name]={'arena':arena_checks,'installer':lambda:installer(k),'primitive':lambda:primitive(k),
                       'b':lambda:b_prefixes(k),'session':lambda:sessions(k),'decoder':lambda:decoded_output(k),'image':lambda:image_decoder(k),'guards':lambda:guards(k),'record-census':record_census,'session4':lambda:sessions_one(k,4),'session9':lambda:sessions_one(k,9)}[name]()
            atomic(args.output,{'completed_groups':out,'all_requested_groups_complete':False})
            print(json.dumps({'group':name,'status':'pass','seconds':time.perf_counter()-t}),flush=True)
    finally:k.close()
    out.update(schema=1,checks='pass',C_released=False,C_complete=False,D_released=False,independent_acceptance=False,
        host={'ru_maxrss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'unit':'bytes on Darwin; KiB on Linux',
              'seconds':time.perf_counter()-started},
        scope_descriptors_inherited=len([x for x in os.environ.get('WORKTREE_SCOPE_FDS','').split(',') if x]),
        source_hashes={p.name:sha(p.read_bytes()) for p in sorted(HERE.glob('*.py'))},
        native_bits_sha256=sha((HERE/'native_bits.c').read_bytes()),
        plan_sha256=sha((HERE/'packed-consumer-plan.json').read_bytes()),
        candidate_manifest_sha256=sha((HERE/'theory-candidate.json').read_bytes()),
        theory_sha256=sha((args.theory_root/'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md').read_bytes()))
    atomic(args.output,out)
    print(json.dumps({'result':str(args.output),'sha256':sha(args.output.read_bytes())}),flush=True)

if __name__=='__main__':main()
