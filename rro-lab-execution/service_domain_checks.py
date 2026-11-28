"""Actual installed source-family continuations and refusal falsifiers.

Source preparations below are explicit test inputs. Full injections are applied
once in the preparer and never passed to a service, decoder or evaluator.
"""
from fractions import Fraction as F
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from field_service import service_graph
from service_native import ServiceEngine,OriginalConnection
from service_checks import cache,receipts,refused,source_unchanged
from session import OriginalPort,prepare_source
from framing import u,rational,v,parse_reply,REPLIES,parse_arithmetic
from field_io import observe_buffer
from native_install import BitKernel
from pins import suppliers
HERE=Path(__file__).resolve().parent
assert sys.version_info>=(3,12) and __debug__


def exercise(k,n,bits,recipe,grant=True):
    source=prepare_source(bits);port=OriginalPort(source);port.grant=grant
    initial=(source.logical(),tuple(source.events),source.routing())
    e=ServiceEngine(service_graph(),22000,kernel=k,connection=OriginalConnection(port))
    cuts={}
    def observe(kind):
        if kind=='original_calls' and 'actual_return_latch' not in cuts:
            cell=e.registers['port.cell']
            assert cell.bit(2) and not receipts(e)
            assert cache(e)==[(0,0)]*(15*n)
            cuts['actual_return_latch']=dict(occupied=cell.bit(2),bits=[cell.bit(0),cell.bit(1)],
                captured=[e.bits['reply.0'],e.bits['reply.1']],retained=e.service_census())
        if kind=='bit_register_writes' and e.cb['occupied'].bit(0):
            label='port_bit_occupied' if len(e.pc.args)==2 else 'output_bit_occupied'
            if label not in cuts:cuts[label]=dict(occupied=1,bx=e.cb['bx'].bit(0),native_steps=e.native_steps)
    e.observer=observe
    e.run();e.feed(u(n)+rational(F(1,3)))
    if not grant:
        refused(lambda:e.feed('01'+v('00')))
        assert not port.calls and not receipts(e) and cache(e)==[(0,0)]*(15*n)
        assert e.native_steps==e.charges['native_controller_selections'] and e.charges['native_refusals']==1
        return dict(recipe=recipe,grant=False,source_calls=0,receipts=0,refusal=True,
                    retained=e.service_census())
    addresses=[0,15,0,15*n-1,2,1,15*n-1]
    expected=[]
    for b in addresses:
        w='R'*b+'LL';before=len(e.output);e.feed('01'+v(w.translate(str.maketrans('LR','01'))))
        result=REPLIES[1 if bits[b]==0 else 2]
        assert parse_reply(''.join(map(str,e.output[before:])))==(2,format(REPLIES.index(result),'02b'))
        expected.append((w,result));assert receipts(e)==expected
        assert [(a,b) for a,b,_,_ in port.calls]==expected
    for b in set(addresses):assert cache(e)[b]==(1,bits[b])
    first=e.registers['cache'].refs[1];other=first
    for _ in range(15):other=other.refs[1]
    assert cache(e)[0]==cache(e)[15] and first is not other
    before_calls=len(port.calls);e.feed('100'+v('00'))
    assert len(port.calls)==before_calls and e.bits['eval.ready']
    _,_,d,_=parse_arithmetic(observe_buffer(e.registers['eval.input']))
    assert all(dict(d[j]).get(h)==bits[15*j+h] for j in range(n) for h in range(15) if 15*j+h in addresses)
    source_unchanged(port,initial)
    return dict(recipe=recipe,N=n,source_bits_sha256=hashlib.sha256(bytes(bits)).hexdigest(),
                invocations=len(expected),unique_addresses=len(set(addresses)),
                actual_ordered_receipts=expected,attained_cache=cache(e),
                predict_new_calls=0,predict_actual_input_sha256=hashlib.sha256(observe_buffer(e.registers['eval.input']).encode()).hexdigest(),
                retained=e.service_census(),pending_cuts=cuts,native_steps=e.native_steps,field_actions=e.ticks,
                installation_steps=e.installation.steps,source_unchanged=True)


def returned_output(k):
    from field_graph import Graph,install_b,loader,growth
    from field_service import Service
    from field_io import output_decoder,print_buffer,literal
    from copier import initial_word
    from nwmg import encode
    rows=[]
    for word in ('','01'):
        g=Graph();s=Service(g);g.label('start');s.header('port.cell');g.j('load')
        loader(g,ready='decode.output');output_decoder(g)
        g.label('decoded.done');s.header('answer');literal(g,'answer','01');s.wrap('decoded','answer')
        print_buffer(g,'answer','finished');g.label('finished');g.emit('wait');g.j('finished')
        growth(g);g.label('fault');g.emit('fault');g.label('exhausted');g.emit('fault')
        p=g.seal();bits=encode(initial_word(''.join('01' if b=='0' else '10' for b in word)+'11'))
        e,_=install_b((),bits,24000,p,lambda *a,**kw:ServiceEngine(*a,kernel=k,**kw));e.run()
        actual=''.join(map(str,e.output))
        assert parse_reply(actual)==(1,word) and observe_buffer(e.registers['decoded'])==word
        assert e.charges['output_ends']==e.charges['acknowledgements']==e.charges['returns']==1
        rows.append(dict(actual_word=word,actual_frame=actual,native_steps=e.native_steps,
                         field_actions=e.ticks,installation_steps=e.installation.steps,retained=e.census()))
        del e;gc.collect()
    return dict(cases=rows,scope='installed actual comb/header decoder plus service Arithmetic framing; arbitrary output-word component fixture, not evaluator execution')


def compatibility(k):
    import ast,copy,subprocess
    from native_engine import NativeEngine
    from field_protocol_checks import exchange_graph
    path=HERE/'native_engine.py'
    prior=subprocess.check_output(['git','-C',str(HERE.parent),'show','HEAD:rro-lab-execution/native_engine.py'],
        pass_fds=tuple(int(x) for x in os.environ.get('WORKTREE_SCOPE_FDS','').split(',') if x))
    current=path.read_bytes();tree=ast.parse(current)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='NativeEngine')
    factory=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='make_core')
    assert ast.dump(factory.body[0].value)==ast.dump(ast.parse('Core(r,b).build()',mode='eval').body)
    cls.body=[n for n in cls.body if n is not factory and not (isinstance(n,ast.Assign) and
               isinstance(n.targets[0],ast.Name) and n.targets[0].id=='native_operations')]
    class Normalize(ast.NodeTransformer):
        def visit_Attribute(self,node):
            node=self.generic_visit(node)
            if isinstance(node.value,ast.Name) and node.value.id=='self' and node.attr=='native_operations':
                return ast.Name(id='NOPS',ctx=ast.Load())
            return node
        def visit_Call(self,node):
            if isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id=='self' and node.func.attr=='make_core':
                assert [x.id for x in node.args]==['cpu_r','cpu_b','source_r']
                return ast.parse('Core(cpu_r,cpu_b).build()',mode='eval').body
            return self.generic_visit(node)
    normalized=Normalize().visit(tree)
    assert ast.dump(normalized)==ast.dump(ast.parse(prior))
    e=NativeEngine(exchange_graph(),1,kernel=k);e.run();assert e.output==[1]
    assert e.charges['output_ends']==e.charges['acknowledgements']==e.charges['returns']==1
    assert e.charges['code_bit_reads']==21*e.ticks and e.charges['code_reference_reads']==4*e.ticks
    saved={}
    for name in ('native-results.json','field-results.json','field-protocol-results.json','paid-code-results.json'):
        raw=(HERE/name).read_bytes();saved[name]=dict(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
    return dict(default_native_engine_fully_normalized_AST_equal=True,
                old_source_sha256=hashlib.sha256(prior).hexdigest(),current_source_sha256=hashlib.sha256(current).hexdigest(),
                fresh_actual_output=e.output,fresh_native_steps=e.native_steps,
                unchanged_component_results=saved,
                scope='default engine factory refactor only; specialized service has separate current author checks')


def main():
    old,_,_=suppliers();k=BitKernel();start=time.perf_counter();rows=[]
    try:
        compatible=compatibility(k);print('default native compatibility pass',flush=True)
        for branch in (0,1):
            graph=old.literal(branch)
            # The full maps are identity injections on all literal endpoints.
            left_labels=list(range(72));right_labels=list(range(243))
            edges=[(0,min(graph[0])),(0,min(graph[0])),(1,min(graph[1])),(2,min(graph[2]))]
            bits=[int(b) for a,c in edges for b in format(left_labels[a],'07b')+format(right_labels[c],'08b')]
            recipe=dict(branch=branch,edges=edges,left_full_injection=left_labels,right_full_injection=right_labels,
                        injection_samples=1,observer_receives_maps=False)
            rows.append(exercise(k,4,bits,recipe));print('shared full injection '+str(branch)+' pass',flush=True);gc.collect()
        # Nine distinct left words, all with the same right word. A right literal
        # endpoint has only eight left neighbors in either supplier; injectivity
        # forces nine distinct left endpoints at one right endpoint, impossible.
        bits=[int(b) for j in range(9) for b in format(j,'07b')+'0'*8]
        rows.append(exercise(k,9,bits,dict(all_promise_array_outside_both_prior_supports=True,
                    witness='nine distinct complete left labels at one common complete right label; right degree eight')))
        print('outside-prior supplied source pass',flush=True);gc.collect()
        output=returned_output(k);print('actual returned-head service framing pass',flush=True)
        rows.append(exercise(k,4,[0]*60,dict(initial_missing_Read_grant=True),grant=False))
        print('missing grant refusal pass',flush=True)
    finally:k.close()
    result=dict(schema=1,checks='pass',default_engine_compatibility=compatible,returned_output=output,instances=rows,seconds=time.perf_counter()-start,
                B_complete=False,C_released=False,D_complete=False,independent_acceptance=False,
                strict_original_omissions=['ML-N-evaluator-execution','ML-Cost-full-evaluator-session'],
                scope='actual installed Read and attained Predict input on specified same-source configurations; no full evaluator arithmetic return')
    result['source_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
        ('service_domain_checks.py','service_checks.py','field_service.py','service_native.py','native_engine.py','session.py','framing.py')}
    path=Path(sys.argv[1]);tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');tmp.replace(path)
    print(json.dumps(dict(checks='pass',result_sha256=hashlib.sha256(path.read_bytes()).hexdigest())))

if __name__=='__main__':main()
