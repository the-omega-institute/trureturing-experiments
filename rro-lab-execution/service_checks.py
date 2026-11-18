"""Exact author checks of the installed received ML-F/original-port path.

The frontend-only stop is identified explicitly; it never returns an arithmetic
answer. --full-image enters the actual unchanged A service path through its
converted/initial-loader prefix. Neither is a full evaluator/session C run.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import gc
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import resource
import sys
import tempfile
import time
from field_graph import Engine,OPS
from field_service import service_graph
from service_native import ServiceEngine,OriginalConnection
from native_install import BitKernel
from session import OriginalPort,prepare_source,classify
from framing import u,v,rational,arithmetic_input,parse_arithmetic,parse_command,reply,REPLIES,parse_reply,cylinder_word
from field_io import observe_buffer
HERE=Path(__file__).resolve().parent
assert sys.version_info>=(3,12) and __debug__


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cache(e):
    p=e.registers['cache'].refs[1];out=[]
    while p is not None:
        out.append((p.bit(0),p.bit(1)) if hasattr(p,'bit') else (p.bits&1,(p.bits>>1)&1))
        p=p.refs[1]
    return out

def receipts(e):
    p=e.registers['receipts'].refs[1];out=[];ids=[]
    from framing import Reader
    while p is not None:
        ids.append(id(p));r=Reader(observe_buffer(p));word=r.v().translate(str.maketrans('01','LR'))
        code=r.take(2);r.end();out.append((word,REPLIES[int(code,2)]));p=p.refs[3]
    assert len(set(ids))==len(ids)
    return out

def refused(fn):
    try:fn()
    except (RuntimeError,AssertionError):return True
    raise AssertionError('expected lawful refusal')


def source_unchanged(port,initial):
    assert port.root is port.source.root and port.grant
    assert initial==(port.source.logical(),tuple(port.source.events),port.source.routing())


def installed_reads(k,bit,mutation=None):
    p=service_graph();source=prepare_source([bit]*60);port=OriginalPort(source)
    initial=(source.logical(),tuple(source.events),source.routing())
    prep=Counter(source.a.charges);connection=OriginalConnection(port)
    e=ServiceEngine(p,20000,kernel=k,connection=connection)
    if mutation=='omit-port-bit':
        n=next(i for i,n in enumerate(p.nodes) if n.op=='emit' and len(n.args)==2)
        e.code[n].opcode=OPS.index('goto')
    if mutation=='omit-cache-write':
        n=next(i for i,n in enumerate(p.nodes) if n.op=='write' and n.args==('class.slot',0,'cache.one'))
        e.code[n].opcode=OPS.index('goto')
    cuts={}
    def observe(kind):
        if kind=='field_commits':
            for name in ('port.call','receipt','receipt.linked','read.reply','print.answer'):
                if name not in cuts and e.pc is p.entries[name]:
                    cuts[name]=dict(field_actions=e.ticks,native_steps=e.native_steps,
                        calls=connection.calls,receipts=len(receipts(e)),
                        captured_reply=[e.bits['reply.0'],e.bits['reply.1']],
                        current_cache=cache(e),retained=e.service_census())
    e.observer=observe
    e.run();e.feed(u(4)+'0'+u(0)+u(2))
    assert e.bits['initialized']==1 and not e.output
    assert cache(e)==[(0,0)]*60
    # Proper received prefix and occupied input latch are actual installed cuts.
    word='01'+v('00');e.feed(word[:-1],seal=False)
    assert not port.calls and not e.output
    e.deliver(end=True);assert not port.calls and not e.output
    if mutation=='omit-port-bit':
        try:e.feed(word[-1:])
        except RuntimeError:
            assert not port.calls and not receipts(e) and not e.registers['port.cell'].bit(4)
            raise AssertionError('actual_port_transfer_occupied')
    else:e.feed(word[-1:])
    expected=[('LL',REPLIES[1 if bit==0 else 2])]
    assert receipts(e)==expected==[(w,x) for w,x,_,_ in port.calls],'authentic_operand_receipt'
    assert cache(e)[0]==(1,bit),'actual_cache_write'
    assert ''.join(map(str,e.output))==reply(2,format(REPLIES.index(expected[0][1]),'02b'))
    operands=['','LL','LRR','LLR','LL','LR','LLL','LRL']
    operands+=['R'*j+'LL' for j in range(1,15)]+['R'*15+'LL','R'*45+'LL','LL']
    operands+=[''.join(w) for n in range(5) for w in product('LR',repeat=n)]
    operands+=['R'*60+s for s in ('','L','R','LL','LR','LLL','LLR','RRR','LLLL')]
    operands+=['L'*90,'R'*103,'R'*59+'LRRRR']
    outputs=[]
    for w in operands:
        before=len(e.output);e.feed('01'+v(w.translate(str.maketrans('LR','01'))))
        got=''.join(map(str,e.output[before:]));kind,value=parse_reply(got)
        b,pair=classify(4,w);val=pair[bit] if b is not None else pair[0]
        assert kind==2 and value==format(REPLIES.index(val),'02b'),(w,got,val)
        expected.append((w,val));assert receipts(e)==expected
        assert expected==[(a,b) for a,b,_,_ in port.calls]
        outputs.append(got)
    c=cache(e);assert all(c[j]==(1,bit) for j in range(15))
    assert c[15]==c[45]==(1,bit) and c[30]==(0,0)
    assert len(expected)==connection.calls==e.charges['original_calls']
    assert connection.a.charges['original_calls']==len(expected) # every repeat/constant/absence
    assert set(x for _,x in expected)==set(REPLIES)
    assert e.charges['output_ends']==e.charges['acknowledgements']==e.charges['returns']==len(expected)
    # Bad sealed syntax returns an actual BadFrame and cannot call the source.
    for bad in ('110','01'+v('')+'1','00'+v(''),'00'+v(u(3)+rational(0))):
        before=len(e.output);calls=connection.calls;e.feed(bad)
        assert ''.join(map(str,e.output[before:]))==reply(0) and connection.calls==calls
    # Close goes through actual output, End/Ack/Return, then a permanent gate.
    before=len(e.output);e.feed('11')
    assert ''.join(map(str,e.output[before:]))==reply(3) and e.bits['closed']
    refused(lambda:e.feed('01'+v('00')))
    source_unchanged(port,initial)
    assert e.native_steps==sum(e.native_events.values())==e.charges['native_controller_selections']
    assert e.charges['code_bit_reads']==21*e.ticks and e.charges['code_reference_reads']==4*e.ticks
    assert e.native_events['portsend']==sum(len(w) for w,_ in expected)
    return dict(source_bit=bit,graph_sha256=p.fingerprint(),graph_nodes=len(p.nodes),
                invocations=len(expected),distinct_addresses=len(set(w for w,_ in expected)),
                replies=sorted(set(x for _,x in expected)),receipt_payload_bits=sum(2*len(w)+3 for w,_ in expected),
                partial_and_equal_index_cache=c,pending_cuts=cuts,receipts=expected,output_sha256=hashlib.sha256(bytes(e.output)).hexdigest(),
                field_actions=e.ticks,native_steps=e.native_steps,installation_steps=e.installation.steps,
                native_events=dict(e.native_events),operations=dict(e.charges),
                provider_operations=dict(source.a.charges-prep),connection_operations=dict(connection.a.charges),
                actual_common_retained_cut=e.service_census(),
                prefix_wait=True,malformed_sealed_rejected=True,close_gate=True,source_unchanged=True)


class FieldInput:
    # Independent field-transition corroboration; never supplies a production
    # native service operation or numerical result.
    def __init__(self):self.e=Engine(service_graph(),50000)
    def run(self):self.e.run()
    def feed(self,word,seal=True):
        e=self.e
        for b in word:
            assert e.waiting and not e.bits['in.ready']
            e.bits['in.bit']=int(b);e.bits['in.end']=0;e.bits['in.ready']=1;e.resume();e.run()
        if seal:
            e.bits['in.end']=1;e.bits['in.ready']=1;e.resume();e.run()


def parser_reference():
    full=tuple(tuple((h,(h+j)%2) for h in range(15)) for j in range(4))
    args=[(4,F(0),tuple(() for _ in range(4)),[]),(4,F(1),full,[None,full]),
          (5,F(1,3),tuple(((14,j%2),) for j in range(5)),[None]),
          (73,F(1,7),tuple(((0,j%2),) for j in range(73)),[])]
    words=[arithmetic_input(*a) for a in args]
    # Preserve actual >64-bit and unreduced rational representations.
    words += [u(4)+'0'+u(2**80)+u(2**81)+'0'*60+u(0),
              u(4)+'0'+u(6)+u(6)+'0'*60+u(0)]
    valid=[]
    for word in words:
        q=FieldInput();q.run();q.feed(u(4)+rational(F(1,3)))
        q.feed('00'+v(word));e=q.e
        assert e.bits['eval.ready'] and observe_buffer(e.registers['eval.input'])==word
        assert parse_arithmetic(word)
        valid.append(dict(bits=len(word),sha256=hashlib.sha256(word.encode()).hexdigest()))
    # Rejected complete nested words, canonical naturals, denominator and range.
    invalid=[u(3)+rational(0),v('00')+rational(0),u(4)+'0'+u(0)+u(0),
             u(4)+'0'+u(2)+u(1),u(4)+'1'+u(1)+u(2),
             u(4)+rational(0)+'0'*60+u(1)+v(''),
             u(4)+rational(0)+'0'*60+u(1)+v('1'+'0'),
             u(4)+rational(0)+'0'*60+u(1)+v('0')]
    for word in invalid:
        q=FieldInput();q.run();q.feed(u(4)+rational(0));q.feed('00'+v(word))
        assert ''.join(map(str,q.e.output))==reply(0) and not q.e.bits.get('eval.ready',0)
    return dict(scope='independent field-transition grammar checks, not native/evaluator execution',
                accepted=valid,rejected_complete_nested=len(invalid))


def installed_predict(k,variant):
    connection=OriginalConnection(OriginalPort(prepare_source([0]*60))) if not variant else None
    e=ServiceEngine(service_graph(),15000,kernel=k,connection=connection);e.run();e.feed(u(4)+rational(F(1,3)))
    if not variant:
        for w in ('LL','R'*15+'LL','LL'):e.feed('01'+v(w.translate(str.maketrans('LR','01'))))
    before_calls=e.charges['original_calls']
    cylinders=[None,tuple(() for _ in range(4)),(((0,1),),(),(),())]
    if variant:
        from framing import words
        command='101'+words([cylinder_word(c) for c in cylinders])
    else:command='100'+v('0'*2)
    e.feed(command);assert e.bits['eval.ready'] and e.bits['eval.busy']
    framed=observe_buffer(e.registers['eval.input']);n,theta,d,cs=parse_arithmetic(framed)
    assert (n,theta)==(4,F(1,3))
    assert d==(tuple(() for _ in range(4)) if variant else (((0,0),),((0,0),),(),()))
    if variant:assert cs==cylinders
    else:
        assert cs==[None,(((0,0),),(),(),()),(((0,1),),(),(),()),None]
    assert e.charges['original_calls']==before_calls
    refused(lambda:e.deliver(0))
    return dict(variant=variant,input_bits=len(framed),input_sha256=hashlib.sha256(framed.encode()).hexdigest(),
                actual_evaluator_frame=framed,prior_actual_calls=before_calls,predict_new_calls=0,busy_refusal=True,
                field_actions=e.ticks,native_steps=e.native_steps,installation_steps=e.installation.steps,
                retained=e.service_census(),scope='installed Predict frontend reaches backed evaluator entry; no numerical answer returned')


def full_image(k):
    from paid_code import CodeStore
    from constructor_pack import unpack
    from constructor_compiler import initial_binary
    from copier import initial_word
    from nwmg import encode
    with tempfile.TemporaryDirectory(prefix='rro-service-image-') as td:
        td=Path(td);unpack(HERE/'constructor-code','W_U.block.gz',td/'W_U.block.gz')
        raw=gzip.decompress((td/'W_U.block.gz').read_bytes())
        expected=json.loads((HERE/'constructor/INTERFACE.json').read_text())['W_U']
        assert len(raw)==expected['uncompressed_image_bytes']
        assert hashlib.sha256(raw).hexdigest()==expected['uncompressed_image_sha256']
        store=CodeStore(raw);p=service_graph(store.height)
        e=ServiceEngine(p,50000,bindings={'image.root':store.root},store=store,kernel=k)
        e.run();e.feed(u(4)+rational(F(1,3)))
        data=tuple(() for _ in range(4));word=arithmetic_input(4,F(1,3),data,[])
        e.feed('00'+v(word),seal=False)
        e.bw('in.bit',0);e.bw('in.end',1);e.bw('in.ready',1);e.resume()
        while e.pc is not p.entries['converted.done']:e.tick()
        converted=observe_buffer(e.registers['converted'])
        assert converted==encode(initial_word(initial_binary(word,361)))
        assert observe_buffer(e.registers['eval.input'])==word
        # Advance the real loader's root setup, then stop before its first cell.
        while e.pc is not p.entries['load.loop']:e.tick()
        root=e.registers['root'];assert root.refs[1] is not None and root.refs[2] is not None
        assert not e.output and e.bits['eval.busy'] and e.charges['original_calls']==0
        return dict(graph_nodes=len(p.nodes),graph_sha256=p.fingerprint(),
                    unchanged_image_bytes=len(raw),unchanged_image_sha256=hashlib.sha256(raw).hexdigest(),
                    actual_received_input=word,converted_bits=len(converted),
                    converted_sha256=hashlib.sha256(converted.encode()).hexdigest(),
                    field_actions=e.ticks,native_steps=e.native_steps,installation_steps=e.installation.steps,
                    retained=e.service_census(),operations=dict(e.charges),
                    scope='actual complete A image, received validation, connected conversions and initial real loader prefix; not a full evaluator run')


def main():
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path);a.add_argument('--full-image',action='store_true')
    a.add_argument('--mutation',choices=('omit-port-bit','omit-cache-write'));args=a.parse_args()
    start=time.perf_counter();k=BitKernel()
    try:
        if args.mutation:
            installed_reads(k,0,args.mutation)
            raise AssertionError('mutation unexpectedly survived')
        result={'schema':1,'stage':'author implementation','B_complete':False,'C_released':False,
                'D_complete':False,'independent_acceptance':False,'whole_goal_complete':False,
                'strict_original_omissions':['ML-N-evaluator-execution','ML-Cost-full-evaluator-session']}
        result['field_grammar']=parser_reference();print('field grammar pass',flush=True)
        result['native_reads']=[]
        for bit in (0,1):
            t=time.perf_counter();result['native_reads'].append(installed_reads(k,bit))
            print(json.dumps(dict(native_reads=bit,seconds=time.perf_counter()-t,status='pass')),flush=True);gc.collect()
        result['native_predict']=[]
        for variant in (0,1):
            result['native_predict'].append(installed_predict(k,variant));print('native predict '+str(variant)+' pass',flush=True);gc.collect()
        if args.full_image:result['production_full_image_prefix']=full_image(k);print('full A service prefix pass',flush=True)
    finally:k.close()
    result['seconds']=time.perf_counter()-start
    result['host_peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['source_hashes']={n:digest(HERE/n) for n in ('field_service.py','service_native.py','service_checks.py','native_engine.py','native_install.py','native_bits.c','field_graph.py','field_decoder.py','field_io.py','session.py','framing.py','constructor/INTERFACE.json')}
    if args.output:
        temp=args.output.with_suffix('.json.tmp');temp.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');temp.replace(args.output)
    print(json.dumps(dict(checks='pass',seconds=result['seconds'],result_sha256=digest(args.output) if args.output else None)),flush=True)

if __name__=='__main__':main()
