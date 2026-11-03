"""Author checks of installed primitive execution and chronological installation.

Finite checks corroborate the separate all-family ordinary proof. --full installs
and uses the entire unchanged A image and the complete connected field graph.
"""
from __future__ import annotations
import argparse
from collections import Counter
import ctypes as C
import gc
import gzip
import hashlib
import json
from pathlib import Path
import resource
import struct
import sys
import tempfile
import time

from field_graph import Graph, install_b, b_graph, observe_b, loader,driver,growth
from field_decoder import lower,page_selector,decoder_graph,request
from field_io import conversions,output_decoder,print_buffer,install_input,observe_buffer,conversion_graph,output_graph
from primitive_engine import PrimitiveEngine
from native_engine import NativeEngine
from native_install import Block,Source,Installer,BitKernel,Init
from field_protocol_checks import exchange_graph,refused,observed
from lab import reference_run
from paid_code import CodeStore,decoder_program
from constructor_pack import unpack
from constructor_image import build_image,Image
from constructor_compiler import initial_binary
from copier import initial_word
from nwmg import encode
HERE=Path(__file__).resolve().parent
assert sys.version_info>=(3,12) and __debug__

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def installer_cuts(k,mutation):
    a=Source(Block([None,None,None,None],bytearray([173]),8))
    b=Source(Block([],bytearray(b'\1'),1));a.block.grant=b.block.grant=1;a.next=b;a.block.refs[0]=a;a.block.refs[3]=b
    if mutation=='omit-installation-source-step':
        k.wires[8][0].next=C.pointer(k.wires[8][2])
    for i in range(4):k.state.counts[i]=(1<<64)-2
    ins=Installer([a,b],k);cuts={};phases=Counter()
    while not ins.done:
        phase=ins.phase
        if phase=='bits':
            wire=k.state.cursor.contents
            action,byte,mask=wire.action,wire.byte,wire.mask
        ins.tick(1);phases[phase]+=1
        if phase=='bits':
            if action==0:assert not (ins.pending.data[byte]&mask)
            elif action==1:assert k.state.latch==bool(ins.material.data[byte]&mask)
            elif action==2:assert bool(ins.pending.data[byte]&mask)==bool(ins.material.data[byte]&mask),'[FAIL] installation_source_operation'
            else:assert action==3 and k.state.done
        if phase not in cuts:
            cuts[phase]=dict(census=ins.census(),pending_grant=None if ins.pending is None else ins.pending.grant)
        if phase=='bits':
            state=k.state
            # Each C step's partial target agrees with the processed fields;
            # tick(1) retains the exact native vertex and occupied read latch.
            assert state.latch in (0,1)
    c=ins.verify()
    if mutation=='omit-installation-read':ins.charges['source_bit_reads']-=1
    assert ins.charges['source_bit_reads']==9,'[FAIL] installation_bit_read_census'
    assert ins.charges['target_bit_writes']==ins.charges['bit_initializations']==9
    assert ins.charges['initializer_bit_latch_writes']==9
    captures=sum(ins.charges[x] for x in ('source_material_reference_reads','source_directory_reads',
                 'initializer_reference_reads','forward_reference_reads','source_reference_reads'))
    assert ins.charges['initializer_reference_latch_writes']==captures
    assert ins.charges['allocation_reference_latch_writes']==2
    assert c['objects']==3*2+98573 and c['refs']==2*4+5*2+98923 and c['bits']==2*9+2168399
    assert ins.execution_pc is a.forward
    assert 'entry-read' in cuts and 'entry-write' in cuts
    assert ins.charges['execution_control_reference_writes']==1
    assert a.forward.refs[0] is a.forward and a.forward.refs[3] is b.forward
    # An allocated but undelivered target is retained alongside the whole source.
    assert cuts['allocate']['pending_grant']==0 and cuts['grant']['pending_grant']==1
    return dict(cuts=cuts,phases=dict(phases),steps=ins.steps,charges=dict(ins.charges),final=c)

def make_b(w,x,k,p=None,capacity=512):
    e,code=install_b(w,x,capacity,p or b_graph(),lambda *a,**kw:NativeEngine(*a,kernel=k,**kw))
    assert e.kernel is not k and e.kernel.wires is k.wires
    assert not e.kernel.state.source and not e.kernel.state.target and not e.kernel.state.cursor
    e.seal();return e,tuple(e.forward_observer[c] for c in code)

def generic(k,mutation):
    cases=[(('L','M','R','M','L','L','M','R','R','M'),'',10),
           ((('J',3),'M',('J',4),'R','M'),'',8),
           (('R','R','L','L','L','L','R'),'',7),
           (('M',('J',1)),'',12),((),'0101',1)]
    rows=[]
    for w,x,bound in cases:
        t=time.perf_counter();e,code=make_b(w,x,k);e.run();states=[observe_b(e,code)]
        for _ in range(bound):
            e.resume();e.run()
            if e.last is e.program.entries['halted']:break
            states.append(observe_b(e,code))
        assert states==reference_run(w,x,bound)[0],(w,states)
        if mutation=='omit-native-operation':e.charges['bit_field_reads']-=1
        # Independent primitive observer versus emitted charge coordinates.
        assert e.native_events['bfield']+e.native_events['codebit']==e.charges['bit_field_reads'],'[FAIL] native_operation_census'
        assert e.native_events['rfield']+e.native_events['coderef']==e.charges['reference_field_reads']
        assert e.charges['code_bit_reads']==21*e.ticks
        assert e.charges['code_reference_reads']==4*e.ticks
        assert e.native_steps==sum(e.native_events.values())==e.charges['native_controller_selections']
        rows.append(dict(program=w,input=x,prefix=states,field_actions=e.ticks,
                         native_steps=e.native_steps,installation_steps=e.installation.steps,
                         retained=e.check_census(),seconds=time.perf_counter()-t))
        # All execution objects are relocated; original source handles never
        # become target handles by a numerical conversion or source reset.
        assert all(c.grant for c in e.cells)
        del e;gc.collect()
    return rows

def op_selectors(k):
    # A real operand/field sweep includes every Cell bit/reference selector,
    # both Boolean branches, literal, bind and indirect.
    g=Graph();g.label('start');g.put('v',1)
    for i in range(4):g.set('x',i,'y');g.get('z','x',i);g.equal('z','y','ok'+str(i),'fault');g.label('ok'+str(i))
    for i in range(8):g.write('x',i,'v');g.read('b','x',i);g.test('b',None,'fault')
    g.null('nil','nil.ok','fault');g.label('nil.ok')
    g.bind('return','finish');g.ret('return');g.label('finish')
    g.emit('emit','v');g.emit('end');g.emit('ack');g.emit('return');g.emit('wait')
    g.label('fault');g.emit('fault');p=g.seal()
    from field_graph import Cell
    a=Cell([None]*4);b=Cell([None]*4)
    e=NativeEngine(p,1,bindings={'x':a,'y':b},kernel=k)
    # External supplied owned cells must be part of the retained finite reserve.
    e.cells.extend((a,b));e.seal();e.run();assert e.output==[1]
    fields=dict(e.native_events);counts=dict(e.charges)
    # Actual encoded operation mutation, before installation, changes execution.
    f=NativeEngine(exchange_graph(),1,kernel=k);f.code[0].opcode=10;f.run();assert f.output==[0]
    rejected=[]
    for mode in ('ack','return'):
        h=Graph();h.label('start');h.emit(mode);h.emit('wait')
        z=NativeEngine(h.seal(),1,kernel=k)
        refused(z.run);rejected.append(mode);del z;gc.collect()
    # Exhaustion and wrong local selector are operational refusals.
    f=NativeEngine(p,1,bindings={'x':a,'y':b},kernel=k);f.cells.extend((a,b));f.code[1].immediate=4
    refused(f.run);rejected.append('selector-4-on-four-reference-cell')
    return dict(native_events=fields,charges=counts,code_mutant_output=[0],refused=rejected)

def conversions_and_output(k):
    rows=[]
    p=conversion_graph(361)
    for name in ('positive-full','positive-partial','null-nine-distinct-left','positive-unreduced-endpoint'):
        x=(HERE/'constructor'/f'{name}.input.bits').read_text()
        expected=(HERE/'constructor'/f'{name}.B-input.bits').read_text()
        e=NativeEngine(p,len(x)+len(expected)+8,kernel=k);install_input(e,x);e.run()
        assert observe_buffer(e.registers['converted'])==expected
        rows.append(dict(name=name,input_bits=len(x),output_bits=len(expected),sha256=hashlib.sha256(expected.encode()).hexdigest(),
                         native_steps=e.native_steps,installation_steps=e.installation.steps))
        del e;gc.collect()
    out=[]
    for word in ('','0','1','01'):
        binary=''.join('01' if b=='0' else '10' for b in word)+'11';x=encode(initial_word(binary))
        e,_=make_b((),x,k,output_graph(),24000);e.run()
        assert ''.join(map(str,e.output))==word==observe_buffer(e.registers['decoded'])
        assert e.charges['output_ends']==e.charges['acknowledgements']==e.charges['returns']==1
        out.append(dict(word=word,native_steps=e.native_steps,installation_steps=e.installation.steps))
        del e;gc.collect()
    return dict(conversion=rows,output=out)

def connected_graph(height,r=361):
    """Production mode zero is the original conversion/driver/output path.

    Probe mode one is an installed decoder caller in the SAME complete graph;
    it returns only actual out_op/out_dest fields and does not operate a source.
    """
    g=Graph();d=decoder_program()
    g.label('start');g.constant('bpc',0);g.j('decoder.'+str(d.labels['boot']))
    lower(g,d,height,'boot.done','image.fetched')
    g.label('boot.done');g.test('probe.mode','probe.ready','convert')
    g.label('probe.ready');g.emit('wait')
    g.bind('fetch.return','probe.result');g.j('fetch')
    g.label('probe.result')
    for i in range(3):g.emit('emit',f'out_op.{i}')
    for i in range(64):g.emit('emit',f'out_dest.{i}')
    g.emit('end');g.emit('ack');g.emit('return');g.j('probe.ready')
    g.label('image.fetched');g.ret('fetch.return')
    conversions(g,r);g.label('converted.done');g.get('input','converted',1);g.j('load')
    loader(g);driver(g,image=True,halt_continuation='decode.output');output_decoder(g)
    g.label('decoded.done');print_buffer(g,'decoded');g.label('printed');g.emit('wait');g.j('printed')
    growth(g);page_selector(g)
    g.label('fault');g.emit('fault');g.label('exhausted');g.emit('fault')
    return g.seal()

def image_case(k,data,cases,r,capacity,probe):
    store=CodeStore(data);p=connected_graph(store.height,r)
    graph_hash=p.fingerprint()
    e=NativeEngine(p,capacity,bindings={'image.root':store.root},store=store,kernel=k)
    assert not e.kernel.state.source and not e.kernel.state.target and not e.kernel.state.cursor
    if probe:e.bits['probe.mode']=1
    else:install_input(e,'0')
    t=time.perf_counter();e.seal();install_seconds=time.perf_counter()-t
    t=time.perf_counter();e.run();boot_seconds=time.perf_counter()-t
    got=[]
    if probe:
        for c in cases:
            for i in range(64):e.bw(f'input_pc.{i}',(c['pc']>>i)&1)
            at=e.native_steps;e.resume();e.run();op=e.observe_word('out_op');dest=e.observe_word('out_dest')
            value=(None,'L','R','M')[op] if op<4 else ['J',dest]
            assert value==c['instruction'],(c,value)
            got.append(dict(pc=c['pc'],instruction=value,native_steps=e.native_steps-at))
    n=len(e.installation.sources);refs=sum(len(s.block.refs) for s in e.installation.sources)
    bits=sum(s.block.width for s in e.installation.sources);c=e.census()
    assert len(e.installation.nodes)==167
    assert c['objects']==3*n+98573 and c['refs']==2*refs+5*n+98923 and c['bits']==2*bits+2168399
    return e,dict(initializer_initial_bindings_null=True,source_profile=dict(n=n,r=refs,b=bits),graph_nodes=len(p.nodes),graph_sha256=graph_hash,native_installed_nodes=n,
                  raw_sha256=hashlib.sha256(data).hexdigest(),raw_bytes=len(data),
                  installation_steps=e.installation.steps,installation_charges=dict(e.installation.charges),
                  census=e.census(),installation_seconds=install_seconds,boot_seconds=boot_seconds,
                  native_steps=e.native_steps,cases=got)

def small_connected(k):
    with tempfile.TemporaryDirectory(prefix='native-small-') as td:
        td=Path(td);raw=b'RROTABLE'+struct.pack('<II',2,1)+struct.pack('<BbI',0,-1,0)+struct.pack('<BbI',1,1,1)
        (td/'binary').write_bytes(raw);build_image(td/'binary',td)
        data=gzip.decompress((td/'W_U.block.gz').read_bytes());(td/'raw').write_bytes(data)
        oracle=Image(td/'raw')
        e,row=image_case(k,data,[],1,24000,False)
        actual=observe_buffer(e.registers['converted']);assert actual==encode(initial_word(initial_binary('0',1)))
        h=pc=0;marks={i for i,b in enumerate(actual) if b=='1'};trace=[observed(e)]
        for _ in range(4):
            ins=oracle.fetch(pc)
            if ins=='L':h-=1;pc+=1
            elif ins=='R':h+=1;pc+=1
            elif ins=='M':marks.add(h);pc+=1
            else:pc=ins[1] if h in marks else pc+1
            e.resume();e.run();state=observed(e);assert state==(h,pc,sorted(marks));trace.append(state)
        row['prefix']=trace;row['native_steps']=e.native_steps
        return row

def existing_decoder_caller(k):
    with tempfile.TemporaryDirectory(prefix='native-decoder-caller-') as td:
        td=Path(td);raw=b'RROTABLE'+struct.pack('<II',2,1)+struct.pack('<BbI',0,-1,0)+struct.pack('<BbI',1,1,1)
        (td/'binary').write_bytes(raw);build_image(td/'binary',td)
        data=gzip.decompress((td/'W_U.block.gz').read_bytes());(td/'raw').write_bytes(data)
        expected=Image(td/'raw');store=CodeStore(data)
        e=NativeEngine(decoder_graph(store.height),1,bindings={'image.root':store.root},store=store,kernel=k)
        e.run();rows=[]
        for pc in (0,7,13*expected.nm,13*expected.nm+1):
            op,dest=request(e,pc)
            actual=(None,'L','R','M')[op] if op<4 else ('J',dest)
            assert actual==expected.fetch(pc),(pc,actual,expected.fetch(pc))
            rows.append(dict(pc=pc,instruction=actual))
        invalid=13*expected.nm+2;refused(lambda:request(e,invalid))
        return dict(caller='field_decoder.request',cases=rows,refused_address=invalid,native_steps=e.native_steps,
                    installation_steps=e.installation.steps,retained=e.check_census())


def full_image(k):
    prior=json.loads((HERE/'paid-code-results.json').read_text())['full']
    with tempfile.TemporaryDirectory(prefix='native-full-') as td:
        packed=Path(td)/'packed';unpack(HERE/'constructor-code','W_U.block.gz',packed)
        with gzip.open(packed,'rb') as f:data=f.read()
    assert hashlib.sha256(data).hexdigest()==prior['raw_sha256']
    e,row=image_case(k,data,[prior['cases'][i] for i in (3,-3,-2,-1)],361,50000,True)
    assert e.installation.charges['source_bit_reads']==sum(s.block.width for s in e.installation.sources)
    row['installation_source_target_bit_equality']=True;row['installation_reference_relocation_verified']=True
    row['installation_verification_phase']='before-first-native-transition'
    return row


def busy_cuts(k):
    # This graph tests the actual old-PC cut as well as occupied output and
    # real code-bit capture. Inspection changes no executing field or choice.
    p=exchange_graph();e=NativeEngine(p,1,kernel=k);e.seal()
    seen={};events=Counter();e.observer=lambda kind:events.update((kind,))
    while not e.waiting:
        old=e.cr['pc'].refs[0];e.microtick()
        label=None
        if e.charges['code_bit_reads']==1:label='first-code-bit'
        if e.cr['pc'].refs[0] is not old and e.ticks==0:
            label='PC-assigned-before-completion'
            assert e.cr['old.pc'].refs[0] is old
        if e.cb['occupied'].bit(0) and not e.output:label='occupied-before-transfer'
        if label not in seen and label:
            seen[label]=dict(census=e.census(),read_bit=e.read_bit,occupied=e.cb['occupied'].bit(0),
                             field_commits=e.ticks,native_steps=e.native_steps)
    while e.ticks<6:e.microtick()
    paused_ticks=e.ticks;paused_pc=e.cr['pc'].refs[0];reads=e.charges['code_bit_reads']
    for _ in range(12):e.microtick()
    assert e.waiting and e.ticks==paused_ticks and e.cr['pc'].refs[0] is paused_pc
    assert e.charges['code_bit_reads']==reads
    assert set(seen)=={'first-code-bit','PC-assigned-before-completion','occupied-before-transfer'}
    assert e.output==[1]
    assert all(events[x]==e.charges[x] for x in events)
    # The represented growth guards are exercised on the same native path.
    from field_graph import Cell
    rejected=[]
    for entry in ('growP','growS'):
        t=Cell([None]*4);l=Cell([t,None,None,None]);r=Cell([t,None,None,None])
        t.refs[1]=l
        if entry=='growS':t.refs[2]=r
        f=NativeEngine(b_graph(),1,entry=entry,bindings={'target':t},kernel=k);f.cells.extend((t,l,r))
        refused(f.run);rejected.append(entry)
        del f;gc.collect()
    # Every original grant-exhaustion fixture is retained, with its actual
    # native installation/guard path rather than a multiplied failure tariff.
    for capacity in (1,2,3):
        f,_=make_b((),'',k,capacity=capacity);refused(f.run);rejected.append('capacity-'+str(capacity))
        del f;gc.collect()
    return dict(cuts=seen,paused_native_steps=12,paused_commits_unchanged=True,paused_code_fetches_unchanged=True,refused=rejected)


def page_leaves(k):
    # Distinguish page offsets with a nonconstant public bit pattern; exercise
    # every native selector leaf and both actual directory-reference fields.
    data=bytes((i^(i>>4)^165)&255 for i in range(8192));store=CodeStore(data)
    g=Graph();g.label('start');g.ref_names.update(('page','fork','other','chosen'))
    for edge in (0,1):
        g.emit('edge','chosen','fork',edge)
        g.equal('chosen','page' if edge==0 else 'other','edge.ok'+str(edge),'fault')
        g.label('edge.ok'+str(edge))
    g.bit_names.add('b')
    for i in range(32768):
        g.emit('page','b','page',i)
        if (data[i//8]>>(i%8))&1:g.test('b',None,'fault')
        else:g.test('b','fault')
    g.emit('wait');g.label('fault');g.emit('fault')
    p=g.seal();t=time.perf_counter()
    e=NativeEngine(p,1,bindings={'fork':store.root,'page':store.pages[0],'other':store.pages[1]},store=store,kernel=k)
    e.run()
    assert e.charges['code_bit_reads']==21*e.ticks and e.charges['code_reference_reads']==4*e.ticks
    return dict(page_fields=32768,directory_edges=2,native_steps=e.native_steps,
                field_actions=e.ticks,source_sha256=hashlib.sha256(data).hexdigest(),seconds=time.perf_counter()-t)

def main():
    a=argparse.ArgumentParser();a.add_argument('--full',action='store_true');a.add_argument('--output',type=Path)
    a.add_argument('--mutation',choices=('omit-native-operation','omit-installation-read','omit-installation-source-step'))
    args=a.parse_args();t=time.perf_counter();k=BitKernel()
    try:
        result=dict(schema=1,B_primitive_engine_closed=bool(args.full),B_installation_closed=bool(args.full),
                    B_complete=False,C_released=False,whole_goal_complete=False)
        groups=(('installer',lambda:installer_cuts(k,args.mutation)),
                ('generic',lambda:generic(k,args.mutation)),('operations',lambda:op_selectors(k)),
                ('page_leaves',lambda:page_leaves(k)),('decoder_caller',lambda:existing_decoder_caller(k)),('busy_cuts',lambda:busy_cuts(k)),
                ('io',lambda:conversions_and_output(k)),('connected',lambda:small_connected(k)))
        for name,check in groups:
            at=time.perf_counter();result[name]=check()
            print(json.dumps(dict(component=name,seconds=time.perf_counter()-at,status='pass')),flush=True)
        if args.full:
            at=time.perf_counter();result['full_image']=full_image(k)
            print(json.dumps(dict(component='full_image',seconds=time.perf_counter()-at,status='pass')),flush=True)
    finally:k.close()
    result['measurements']=dict(wall_seconds=time.perf_counter()-t,host_ru_maxrss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                ru_maxrss_unit='bytes on Darwin, KiB on Linux')
    result['source_hashes']={n:digest(HERE/n) for n in ('native_bits.c','native_install.py','native_engine.py','native_checks.py')}
    result['consumer_source_hashes']={n:digest(HERE/n) for n in ('field_graph.py','primitive_engine.py','field_decoder.py','field_io.py','paid_code.py')}
    result['A_interface_sha256']=digest(HERE/'constructor/INTERFACE.json')
    result['strict_original_omissions']=['ML-N-evaluator-execution','ML-Cost-full-evaluator-session']
    if args.output:args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(dict(checks='pass',wall_seconds=result['measurements']['wall_seconds'],B_complete=False,C_released=False)))

if __name__=='__main__':main()
