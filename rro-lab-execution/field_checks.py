"""Author checks for field-graph execution; no C execution or B acceptance.

--full includes selected full-A-image primitive-engine fetches. Existing A law
runs and original strict failures remain immutable/reused, not rerun or cleared.
"""
from __future__ import annotations
import argparse
from collections import Counter
import gc
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import tempfile

from field_graph import Graph, Engine, OPS, b_graph, install_b, observe_b
from primitive_engine import PrimitiveEngine
from field_io import (conversion_graph, output_graph, observe_buffer, install_input)
from field_decoder import decoder_graph, image_b_graph, request
from lab import reference_run
from constructor_compiler import initial_binary
from copier import initial_word
from nwmg import encode
from paid_code import CodeStore
from constructor_pack import unpack

HERE=Path(__file__).resolve().parent
assert sys.version_info>=(3,12) and __debug__


def b_checks(mutation=None):
    graph=b_graph();counts=Counter();cuts=Counter()
    inputs=[''.join(x) for n in range(3) for x in product('01',repeat=n)]
    for q in range(4):
        instructions=('L','R','M')+tuple(('J',j) for j in range(q))
        for w in product(instructions,repeat=q):
            for bits in inputs:
                e,code=install_b(w,bits,256,graph,PrimitiveEngine)
                e.run();actual=[observe_b(e,code)]
                for _ in range(8):
                    e.resume();e.run()
                    if e.last is graph.entries['halted']:break
                    actual.append(observe_b(e,code))
                expected=reference_run(w,bits,8)[0]
                assert actual==expected,(w,bits,actual,expected)
                counts['instances']+=1;counts['commits']+=len(actual)-1
                counts['graph_actions']+=e.ticks;counts['cpu_subphases']+=e.microsteps
                counts['halted' if e.last is graph.entries['halted'] else 'bounded_prefixes']+=1
    assert counts['instances']==1722
    # One run observes every CPU subphase, not just instruction boundaries.
    w=('L','M','R','M','L','L','M','R','R','M')
    e,code=install_b(w,'',256,graph,PrimitiveEngine)
    observed=Counter();old=None
    def observer(kind):observed[kind]+=1
    e.observer=observer
    while not e.waiting:
        e.microtick();cuts[e.cpu.phase]+=1
    # microtick's wait is committed before the external run scheduler returns.
    e.microtick()
    trace=[observe_b(e,code)]
    while len(trace)<=10:
        e.resume()
        while not e.waiting:
            e.microtick();cuts[e.cpu.phase]+=1
        e.microtick()
        if e.last is graph.entries['halted']:break
        trace.append(observe_b(e,code))
    assert trace==reference_run(w,'',10)[0]
    if mutation=='omit-code-reads':e.charges['code_bit_reads']-=1
    assert dict(observed)=={k:e.charges[k] for k in observed},'[FAIL] primitive_event_census'
    assert e.charges['code_bit_reads']==21*e.ticks,'[FAIL] actual_code_scan'
    assert set(e.bits.values())<={0,1}
    # Independent count walks actual object/field collections, no tariff terms.
    c=e.census()
    assert c['reference_slots']==sum(len(x.refs) for x in e.cells)+sum(4 for x in e.code)+len(e.registers.slots)+8
    assert c['bit_fields']==sum(8 for x in e.cells)+sum(21 for x in e.code)+len(e.bits.slots)+38
    assert c['objects']==len(e.cells)+len(e.code)+len(e.registers.slots)+len(e.bits.slots)+1
    # Removing external observations changes neither control nor material.
    e.output.clear();e.events.clear();e.cuts.clear();before=observe_b(e,code)
    e.resume();e.run();assert before==observe_b(e,code)
    # Bad guards and inadequate finite supply refuse, without inventing a step.
    failed=0
    for capacity in (1,2,3):
        f,_=install_b((),'',capacity,graph,PrimitiveEngine)
        try:f.run()
        except RuntimeError:failed+=1
    assert failed==3
    return {'finite_cases':dict(counts),'subphase_cuts':dict(sorted(cuts.items())),
            'graph_nodes':len(graph.nodes),'graph_sha256':graph.fingerprint(),
            'cross_census':c,'cross_primitive_events':dict(sorted(observed.items())),
            'falsifiers':['exact-B-prefix-comparison','primitive-code-event-census',
                          'inadequate-grant-refuses','binary-bit-domain','audit-erasure']}


def arithmetic_checks():
    cases=[('add',13,29,42),('sub',29,13,16),('eq',0,0,1),('eq',1,0,0),
           ('lt',1,2,1),('lt',2,1,0),('lt',2,2,0),('shl',1,5,32),
           ('shr',128,3,16),('copy',(1<<64)-1,0,(1<<64)-1)]
    result=[]
    for op,a,b,want in cases:
        g=Graph();g.label('start');g.constant('a',a);g.constant('b',b)
        g.alu(op,'a','a','b' if op in ('add','sub','eq','lt') else None,b)
        g.emit('wait');g.label('fault');g.emit('fault')
        p=g.seal();e=PrimitiveEngine(p,1);e.run()
        assert e.observe_word('a')==want,(op,a,b,e.observe_word('a'))
        result.append(dict(op=op,a=a,b=b,result=want,graph_actions=e.ticks,cpu_subphases=e.microsteps))
    return result


def conversion_checks():
    p=conversion_graph(361);result=[]
    for name in ('positive-full','positive-partial','null-nine-distinct-left','positive-unreduced-endpoint'):
        bits=(HERE/'constructor'/f'{name}.input.bits').read_text()
        expected=(HERE/'constructor'/f'{name}.B-input.bits').read_text()
        e=PrimitiveEngine(p,len(bits)+len(expected)+8);install_input(e,bits);e.run()
        actual=observe_buffer(e.registers['converted'])
        assert actual==expected==encode(initial_word(initial_binary(bits,361)))
        result.append(dict(name=name,input_bits=len(bits),B_input_bits=len(actual),
                           sha256=hashlib.sha256(actual.encode()).hexdigest(),
                           graph_actions=e.ticks,cpu_subphases=e.microsteps,
                           occupied=e.census()['occupied_pool_cells']))
    return {'nodes':len(p.nodes),'graph_sha256':p.fingerprint(),'cases':result,
            'scope':'Actual conversion of published raw inputs; no law evaluation, paid Receive or LAB load in these cases.'}


def output_checks():
    p=output_graph();result=[]
    for word in ('','0','1','01'):
        binary=''.join('01' if b=='0' else '10' for b in word)+'11'
        bits=encode(initial_word(binary))
        e,_=install_b((),bits,24000,p,PrimitiveEngine);e.run()
        actual=''.join(map(str,e.output))
        assert actual==word==observe_buffer(e.registers['decoded'])
        assert e.charges['output_bits']==len(word)
        assert e.charges['output_ends']==e.charges['acknowledgements']==e.charges['returns']==1
        result.append(dict(word=word,graph_actions=e.ticks,cpu_subphases=e.microsteps,
                           retained=e.census()))
    return {'graph_nodes':len(p.nodes),'graph_sha256':p.fingerprint(),'cases':result,
            'scope':'Actual returned-head scan of installed comb fixtures and emitted separator-terminated words; not an executed full evaluator.'}


def full_image_checks():
    # Reuse published exact image and expected fetches; never call A reader or
    # rerun the unchanged stack law executions. Inflation is static preparation.
    saved=json.loads((HERE/'paid-code-results.json').read_text())
    with tempfile.TemporaryDirectory(prefix='field-image-') as td:
        packed=Path(td)/'packed';unpack(HERE/'constructor-code','W_U.block.gz',packed)
        with gzip.open(packed,'rb') as stream:data=stream.read()
    assert hashlib.sha256(data).hexdigest()==saved['full']['raw_sha256']
    store=CodeStore(data);del data;gc.collect()
    graph=decoder_graph(store.height)
    graph_hash=graph.fingerprint()
    e=PrimitiveEngine(graph,1,bindings={'image.root':store.root},store=store)
    e.run();boot=e.microsteps;out=[]
    for index in (3,-3,-2,-1):
        case=saved['full']['cases'][index];start=e.microsteps
        op,dest=request(e,case['pc'])
        got=(None,'L','R','M')[op] if op<4 else ['J',dest]
        assert got==case['instruction'],(case,got)
        out.append(dict(pc=case['pc'],instruction=got,cpu_subphases=e.microsteps-start))
    return dict(nodes=len(graph.nodes),graph_sha256=graph_hash,boot_cpu_subphases=boot,
                cases=out,retained=e.census(),primitive_events=dict(sorted(e.charges.items())),
                raw_image_sha256=saved['full']['raw_sha256'],
                scope='Full image installed by host preparation; selected actual primitive-engine code fetches, not chronological installation or law execution.')


def build(full=False,mutation=None):
    result={'schema':1,'B_complete':False,'C_released':False,
            'whole_runtime_state_map_complete':False,
            'B':b_checks(mutation),'ALU':arithmetic_checks(),
            'input_conversion':conversion_checks(),'output_scan':output_checks()}
    if full:result['full_image']=full_image_checks()
    result['source_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ('field_graph.py','primitive_engine.py','field_decoder.py','field_io.py','field_checks.py')}
    result['immutable_A_interface']=hashlib.sha256((HERE/'constructor/INTERFACE.json').read_bytes()).hexdigest()
    result['limits']=['Host installation has no paid initializer continuation.',
                      'Finite-control multibit updates and field-selector decoding still need literal primitive refinement.',
                      'No particular original-port service/controller or joint source/observer timeline.',
                      'Fixed finite supplied pool is fully retained; no dynamic widening, cleanup or shared-job controller.',
                      'No full compiled evaluator run; both original strict failures remain.']
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true');parser.add_argument('--output',type=Path)
    parser.add_argument('--check',type=Path);parser.add_argument('--mutation',choices=['omit-code-reads'])
    args=parser.parse_args();result=build(args.full,args.mutation)
    data=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if args.output:args.output.write_text(data)
    if args.check:assert args.check.read_text()==data,'saved field results differ'
    print(json.dumps(dict(checks='pass',B_complete=False,C_released=False,
                         sha256=hashlib.sha256(data.encode()).hexdigest())))

if __name__=='__main__':main()
