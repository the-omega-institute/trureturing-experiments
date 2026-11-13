"""Discriminating CPU/I-O checks and connected finite-image LAB prefixes.

These are author component checks; they neither supply the original-port
service nor run the full acquired law evaluator.
"""
import gzip
import hashlib
import json
from pathlib import Path
import struct
import sys
import tempfile
from field_graph import Graph, OPS
from primitive_engine import PrimitiveEngine
from field_decoder import execution_graph
from field_io import install_input,observe_buffer
from paid_code import CodeStore
from constructor_image import build_image,Image
from constructor_compiler import initial_binary
from copier import initial_word
from nwmg import encode

HERE=Path(__file__).resolve().parent


def exchange_graph():
    g=Graph();g.label('start');g.put('x',1);g.emit('emit','x')
    g.emit('end');g.emit('ack');g.emit('return');g.emit('wait')
    return g.seal()


def refused(call):
    try:call()
    except (AssertionError,RuntimeError):return True
    raise AssertionError('forbidden protocol completed')


def protocol():
    p=exchange_graph();e=PrimitiveEngine(p,1)
    refused(e.resume)
    while not e.cpu.occupied:e.microtick()
    assert e.output==[] and e.cpu.bx==1
    occupied=(e.cpu.phase,e.cpu.bx,e.cpu.occupied,e.microsteps)
    e.microtick();assert e.output==[1] and not e.cpu.occupied
    e.run();assert e.waiting and e.charges['returns']==1
    f=PrimitiveEngine(p,1);f.code[0].opcode=OPS.index('goto');f.run()
    assert f.output==[0] and e.output!=f.output
    for bad in ('ack','return'):
        g=Graph();g.label('start')
        if bad=='return':g.emit('end')
        g.emit(bad);g.emit('wait')
        refused(PrimitiveEngine(g.seal(),1).run)
    return dict(occupied_cut=occupied,original_output=e.output,code_bit_mutant_output=f.output,
                refused=['resume-before-wait','Ack-before-End','Return-before-Ack'],
                code_opcode_mutation_changes_execution=True,
                primitive_events=dict(sorted(e.charges.items())))


def observed(e):
    root=e.registers['root'];marked=[]
    for arm,key in ((0,2),(1,1)):
        c=root.refs[key];j=0
        while c.refs[2] is not None:
            if c.refs[2].refs[1] is not None:marked.append(-j-1 if arm else j)
            c=c.refs[1];j+=1
    h=e.registers['head'];depth=0
    while h is not root:depth+=1;h=h.refs[0]
    h=-depth if e.bits['arm'] else depth-1
    return h,e.observe_word('bpc'),sorted(marked)


def connected_prefix():
    with tempfile.TemporaryDirectory(prefix='field-prefix-') as td:
        td=Path(td)
        raw=b'RROTABLE'+struct.pack('<II',2,1)+struct.pack('<BbI',0,-1,0)+struct.pack('<BbI',1,1,1)
        (td/'binary').write_bytes(raw);build_image(td/'binary',td)
        data=gzip.decompress((td/'W_U.block.gz').read_bytes());(td/'raw').write_bytes(data)
        oracle=Image(td/'raw');store=CodeStore(data)
        # Same composed graph builder as the full-A path. This check uses a
        # one-register serialization only to exercise its connected mechanics.
        p=execution_graph(store.height,1)
        e=PrimitiveEngine(p,24000,bindings={'image.root':store.root},store=store)
        bits='0';install_input(e,bits);e.run()
        actual=observe_buffer(e.registers['converted'])
        assert actual==encode(initial_word(initial_binary(bits,1)))
        expected_head=expected_pc=0;marked={j for j,b in enumerate(actual) if b=='1'}
        trace=[observed(e)]
        for _ in range(4):
            # Independent reference used only for comparison after code setup;
            # its fetched value is never written into e or its installed code.
            ins=oracle.fetch(expected_pc)
            if ins=='L':expected_head-=1;expected_pc+=1
            elif ins=='R':expected_head+=1;expected_pc+=1
            elif ins=='M':marked.add(expected_head);expected_pc+=1
            else:expected_pc=ins[1] if expected_head in marked else expected_pc+1
            e.resume();e.run();got=observed(e)
            assert got==(expected_head,expected_pc,sorted(marked)),(got,expected_head,expected_pc)
            trace.append(got)
        return dict(graph_nodes=len(p.nodes),graph_sha256=p.fingerprint(),
                    converted_bits=len(actual),prefix=trace,cpu_subphases=e.microsteps,
                    graph_actions=e.ticks,retained=e.census(),
                    scope='Connected conversion/fetch/guarded LAB prefix on a finite component image; no full law evaluation or session.')


def main():
    assert sys.version_info>=(3,12) and __debug__
    result=dict(schema=1,protocol=protocol(),connected=connected_prefix(),B_complete=False,C_released=False)
    names=('field_graph.py','primitive_engine.py','field_decoder.py','field_io.py','field_protocol_checks.py')
    result['source_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    data=json.dumps(result,sort_keys=True,indent=2)+'\n'
    path=Path(sys.argv[1]);path.write_text(data)
    print(json.dumps(dict(checks='pass',sha256=hashlib.sha256(data.encode()).hexdigest())))

if __name__=='__main__':main()
