"""One bounded IS2 service-interface and prospective-certificate bundle.

The installed numerical witness reads N from an accepted original frame and
computes N+1. It is NOT the fixed UJB evaluator or a full acquisition session.
The complete fixed evaluator is generated and linked only; never installed.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent)]
import immutable_spine as im
import constructor_ir as ci
import evaluator_program as ep
from field_graph import Node,Program
from field_service import service_graph
from field_io import observe_buffer
from native_install import BitKernel
from service_native import OriginalConnection
from session import OriginalPort,prepare_source
from framing import Reader,u,v,parse_arithmetic,parse_reply,REPLIES

CAPS=dict(wall_seconds=600,native_steps=4_000_000,installation_events=20_000_000,
          reserve_cells=16_384,installed_source_blocks=100_000,
          static_field_nodes=100_000,accepted_input_bits=1024,
          source_bits=60,source_calls=8,retained_output_bits=4096,
          modeled_objects=500_000,monitored_host_peak_bytes=2*1024**3)
NATIVE_COST={'copy':2,'put':2,'null':2,'equal':4,'test':2,'get':4,'set':4,
             'read':4,'write':4,'bind':2,'indirect':2,'emit':2,'goto':0,
             'end':1,'ack':1,'return':1,'wait':1,'fault':1}


def host_cap():
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    assert rss<=CAPS['monitored_host_peak_bytes'],'monitored resident memory cap'
    return rss


def compose(evaluator):
    """Literal native edges to the unchanged service handlers, shared pool only.

    Namespace renaming changes no value rule. End/Ack/Return use the same
    ordinary output connection as the frontend. Return reads output.return;
    a local acknowledged bit alone cannot choose the service continuation.
    """
    service=service_graph()
    def rename(x):return x if x in ('pool','nil') else 'is2.'+x
    def clone(p,names):
        by={n:Node(n.op,()) for n in p.nodes}
        for n,c in by.items():
            c.args=tuple(by[x] if isinstance(x,Node) else names(x) if isinstance(x,str) else x for x in n.args)
            c.zero=by[n.zero];c.one=by[n.one]
        return list(by.values()),{k:by[n] for k,n in p.entries.items()}
    sn,se=clone(service,lambda x:x);en,ee=clone(evaluator,rename)
    returned=Node('put',('eval.busy',0),se['service.return'],se['service.return'])
    reset=ee['start']
    # Every new invocation starts with empty logical stacks. Published old
    # frames are retained and remain immutable; no clearing of backing cells.
    stack_names=sorted(x for x in evaluator.refs if x.startswith('stack.'))
    added=[returned]
    for name in reversed(stack_names):
        reset=Node('copy',(rename(name),'nil'),reset,reset);added.append(reset)
    bound=Node('bind',('output.return',returned),reset,reset);added.append(bound)
    begin=se['evaluate.begin']
    assert begin.op=='put' and begin.args==('eval.busy',1)
    cursor=begin.zero;assert cursor.op=='get' and cursor.args==('input','eval.input',1)
    ready=cursor.zero;assert ready.op=='put' and ready.args==('eval.ready',1)
    ready.op='get';ready.args=('is2.received','eval.input',1);ready.zero=ready.one=bound
    old_print=se['print.answer']
    print_bind=Node('bind',('output.return',se['service.return']),old_print,old_print)
    for n in sn:
        if n.zero is old_print:n.zero=print_bind
        if n.one is old_print:n.one=print_bind
    se['print.answer']=print_bind;added.append(print_bind)
    nodes=tuple(sn+en+added)
    entries={**se,**{'is2.'+k:n for k,n in ee.items()},'is2.returned':returned}
    refs=tuple(sorted(set(service.refs)|{rename(x) for x in evaluator.refs}|{'output.return'}))
    bits=tuple(sorted(set(service.bits)|{rename(x) for x in evaluator.bits}))
    p=Program(nodes,entries,refs,bits)
    ids={n:i for i,n in enumerate(nodes)}
    im.GRAPH_LABELS[id(p)]={k:ids[n] for k,n in entries.items()}
    assert all(n.op in NATIVE_COST for n in nodes)
    return p


def static_demand(raw,states,registers,routines):
    n,theta,cache,cyl=parse_arithmetic(raw);t=len(cyl);l=len(raw)
    # The represented denominator is read without rational cancellation.
    rd=Reader(raw);assert rd.u()==n;assert rd.take(1)=='0';a=rd.u();b=rd.u()
    # D <= 2^(26 N); every summand denominator divides D. Mixtures have
    # denominator dividing b D and conditional ratios have width <= 2m.
    m=64*(n+l+t+2)
    y=64*(t+2)*(m+1)
    arrays=(16*n+3)*(1<<min(n,243))
    data=l+y+(n+1)**2+arrays+m+1024
    word_bound=16**(routines+1)*(data+1)
    exponent=registers*(word_bound+1)+(states-1).bit_length()+1
    source=[dict(row).get(h,0) for row in cache for h in range(15)]
    schedule=['','LL','LRR','LLR','LL','LR','LLL','LRL']
    schedule+=['R'*i+'LL' for i in range(15*n) if (i%15,source[i]) in cache[i//15]]
    schedule+=['LL','R'*15+'LL','R'*45+'LL','L'*130,'R'*(15*n+10),'LL']
    # Add real post-evaluation old-slot Read and Close, preserving acquisition.
    post='LL';q=len(schedule)+1;s=sum(map(len,schedule))+len(post)
    # Exact received bit count for initialization, Read commands, Arithmetic,
    # then appended old-slot Read and Close; prior prefix is raw[:rd.pos].
    input_bits=rd.pos+sum(2+len(v(w.translate(str.maketrans('LR','01')))) for w in schedule+[post])+2+len(v(raw))+2
    # Frontend intermediate words: unary N/shape, framing/validation words,
    # pending request, receipt and answer lengths. The prefix routine type
    # invariant in the theory supplies this common length envelope.
    frontend_word=64*(l+input_bits+s+15*n+t+1024)**2
    frontend=service_graph();frontend_states=len(frontend.nodes)
    frontend_registers=len(frontend.refs)+len(frontend.bits)
    # Encode the active typed chains, named cursors, cache and receipt links;
    # dead retained material is resource demand, not hidden control state.
    frontend_exponent=(frontend_registers+4)*(32*frontend_word**2+1)+(frontend_states-1).bit_length()+input_bits.bit_length()+8
    return dict(N=n,L=l,t=t,represented_prior=[a,b],input_sha256=hashlib.sha256(raw.encode()).hexdigest(),
        acquisition_calls=len(schedule),acquisition_operand_bits=sum(map(len,schedule)),
        appended_continuation=['Read(LL)','Close'],complete_calls=q,complete_operand_bits=s,
        total_input_bits=input_bits,m=m,U_out=y,logical_word_bound=word_bound,
        U_step_expression='2^E',E=exponent,U_push_expression='2^E',
        frontend_word_bound=frontend_word,frontend_bound_expression='2^E_F',E_F=frontend_exponent,
        K_expression='5+L+2^E+2^E_F',
        whole_native_work_expression='I_install + F_native*(2^E+2^E_F+2*L+2*U_out+64) + provider_work',
        owner_backing=dict(H_reserved=None,W_backed=None,
            absent_owner='existing host/resource owner of complete chronology, including inherited overlapping jobs',
            host_layout_upper_bounds=None,existing_live_original_session_manifest=None))


class CappedInstaller(im.Installer):
    def run(self,checkpoint=None):
        assert len(self.sources)<=CAPS['installed_source_blocks']
        while not self.done:
            assert self.steps<CAPS['installation_events'],'installation cap'
            if self.steps%1000==0:host_cap()
            assert self.accounted['objects']<=CAPS['installed_source_blocks']
            self.tick(min(100_000,CAPS['installation_events']-self.steps))
        return None


class ServiceDirect(im.Direct):
    def __init__(self,p,kernel,connection):
        super().__init__(p,'',CAPS['reserve_cells'],kernel,connection=connection,service_return=True)
        self.handoffs=[];self.return_cuts=[];self.boundary_observer=self.observe
    def ref(self,n):return self.regs[n].refs[0]
    def value(self,n):return self.bitregs[n].bit(0)
    def observe(self,e,label):
        if label=='is2.start':
            actual=self.ref('eval.input').refs[1]
            assert self.ref('is2.received') is actual and self.value('eval.busy')==1
            cursor=actual;cells=[]
            while cursor is not None:
                assert cursor.grant;cells.append(cursor);cursor=cursor.refs[1]
            self.handoffs.append(dict(same_chain=True,current_read_grants=True,bits=''.join(str(x.bit(0)) for x in cells),
                native_steps=self.steps))
            self.frozen.update({x:tuple(x.refs) for x in cells})
        if label=='is2.returned':
            assert self.ended and self.acked and self.value('eval.busy')==1
            assert self.ref('output.return') is self.pc
            self.return_cuts.append(dict(actual_indirect_return=True,End=True,Ack=True,
                busy_before_cleanup=1,output=''.join(self.output),native_steps=self.steps))
        if label.startswith('is2.q.') or label=='is2.print':
            for name in self.regs:
                if name.startswith('is2.stack.'):
                    cell=self.ref(name)
                    while cell is not None and cell not in self.frozen:
                        self.frozen[cell]=tuple(cell.refs);cell=cell.refs[0]
            for name in ('work.root','token.one','token.zero','token.marker','route.header'):
                cell=self.ref('is2.'+name)
                if cell is not None:self.frozen[cell]=tuple(cell.refs)
    def tick(self):
        assert self.steps<CAPS['native_steps'],'native step cap'
        assert len(self.output)<=CAPS['retained_output_bits'],'output cap'
        super().tick()
    def drive(self):
        self.run(CAPS['native_steps']);assert self.done and not self.refused
    def deliver(self,bit=None,end=False):
        assert self.done and not self.value('in.ready') and not self.value('closed') and not self.value('eval.busy')
        for name,value in [('in.bit',int(bit or 0)),('in.end',int(end)),('in.ready',1)]:
            assert self.bitregs[name].grant;self.bitregs[name].put(0,value)
            self.events['paid_input_latch_writes']+=1
        self.done=False;self.events['paid_resume_signals']+=1;self.drive()
    def feed(self,bits):
        for b in bits:self.deliver(int(b))
        self.deliver(end=True)
    def cache(self):
        p=self.ref('cache').refs[1];out=[]
        while p is not None:out.append((p.bit(0),p.bit(1)));p=p.refs[1]
        return out
    def receipts(self):
        p=self.ref('receipts').refs[1];out=[];ids=[]
        while p is not None:
            ids.append(id(p));r=Reader(observe_buffer(p));w=r.v().translate(str.maketrans('01','LR'))
            code=r.take(2);r.end();out.append((w,REPLIES[int(code,2)]));p=p.refs[3]
        assert len(set(ids))==len(ids)
        return out,ids


def tiny_program():
    with im.shared_assembler():
        a=ci.Assembler();ci.primitives(a)
        def main(a):
            a.call('readu',['actual.n','input'],'input')
            a.op('add','computed.n','actual.n',1)
            a.g.copy('computed.n','input')
        return im.prune(a.compile(main))


def native_interface(kernel,raw):
    n,theta,cache,cyl=parse_arithmetic(raw)
    source_bits=[dict(row).get(h,0) for row in cache for h in range(15)]
    assert len(source_bits)<=CAPS['source_bits'] and len(raw)<=CAPS['accepted_input_bits']
    source=prepare_source(source_bits);port=OriginalPort(source);connection=OriginalConnection(port)
    initial=(source.logical(),tuple(source.events),source.routing())
    p=compose(im.field_graph(tiny_program()));e=ServiceDirect(p,kernel,connection);e.drive()
    r=Reader(raw);r.u();r.rational();e.feed(raw[:r.pos])
    # Actual pre-evaluation old slots, four reply classes and a partial cache.
    before_words=['','LL','LRR','LLR','LRL','R'*15+'LL']
    for w in before_words:e.feed('01'+v(w.translate(str.maketrans('LR','01'))))
    before=e.receipts();cache_before=e.cache();out_before=len(e.output)
    assert set(x[1] for x in before[0])==set(REPLIES)
    e.feed('00'+v(raw))
    assert len(e.handoffs)==len(e.return_cuts)==1
    assert e.handoffs[0]['bits']==raw
    computed=''.join(e.output[out_before:]);assert parse_reply(computed)==(1,bin(n+1)[:1:-1])
    assert e.receipts()==before and e.cache()==cache_before and not e.value('eval.busy')
    assert e.connection is connection and connection.port is port and port.root is source.root
    e.feed('01'+v('00'));after=e.receipts()
    assert after[0]==before[0]+[before[0][1]] and after[1][:-1]==before[1]
    assert after[1][-1] not in before[1] and e.cache()==cache_before
    before_close=len(e.output);e.feed('11')
    assert ''.join(e.output[before_close:])=='110' and e.value('closed')
    assert initial==(source.logical(),tuple(source.events),source.routing())
    assert [(w,x) for w,x,_,_ in port.calls]==after[0]
    assert len(port.calls)<=CAPS['source_calls']
    assert e.events['end']==e.events['ack']==e.events['return']==len(port.calls)+2
    c=e.census();assert c['objects']<=CAPS['modeled_objects']
    # Actual current-right refusal at the same input field instruction. Restore
    # no state and do not resume: this object is terminal after Close.
    return dict(scope='native interface witness, computes N+1; no fixed evaluator run or full acquisition schedule',
        graph_sha256=p.fingerprint(),field_nodes=len(p.nodes),native_nodes=len(e.native),
        input_sha256=hashlib.sha256(raw.encode()).hexdigest(),computed_frame=computed,
        native_steps=e.steps,installation_events=e.inst.steps,census=c,
        handoffs=e.handoffs,returns=e.return_cuts,occurrences=after[0],fresh_old_slot_receipt=True,
        cache_before_and_after=cache_before,original_root_preserved=True,source_unchanged=True,
        close=True,events=dict(e.events),source_census=source.a.census(),connection_census=connection.a.census())


def main():
    start=time.monotonic()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('prospective wall cap')))
    signal.alarm(CAPS['wall_seconds'])
    result=dict(caps=CAPS,candidate=im.IDENTITY,status='incomplete',full_fixed_installations=0,
                full_evaluator_executions=0,full_original_or_replacement_sessions=0,checks=[])
    kernel=None
    try:
        with im.shared_assembler():fixed=ep.build()
        p=im.prune(fixed);g=im.field_graph(p);linked=compose(g)
        assert len(linked.nodes)<=CAPS['static_field_nodes']
        a=ci.Assembler();ci.primitives(a);ep.define_evaluator(a);routines=len(a.functions)+1
        result['static']=dict(stack_sha256=hashlib.sha256(im.canonical(p)).hexdigest(),stack_states=len(p['code']),
            stack_registers=len(p['registers']),routines_including_main=routines,
            composed_field_nodes=len(linked.nodes),composed_refs=len(linked.refs),composed_bits=len(linked.bits),
            composed_native_nodes=len(linked.nodes)+sum(NATIVE_COST[n.op] for n in linked.nodes),
            composed_field_sha256=linked.fingerprint())
        result['prospective']=[static_demand((HERE.parent/'constructor'/(name+'.input.bits')).read_text(),
            len(p['code']),len(p['registers']),routines) for name in
            ('positive-full','positive-partial','null-nine-distinct-left')]
        result['checks'].append('static full-source linking and three prospective expressions')
        im.Installer=CappedInstaller
        kernel=BitKernel()
        result['interface']=native_interface(kernel,(HERE.parent/'constructor/positive-full.input.bits').read_text())
        result['checks'].append('actual tiny native service input/computation/End/Ack/indirect Return/old Read/receipt/cache/Close')
        result['bundle_exit']=0;result['status']='bounded-interface-complete; full-session-resource-offer-absent'
    except BaseException as exc:
        result['bundle_exit']=1;result['error']=dict(type=type(exc).__name__,message=str(exc),
            locations=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(exc.__traceback__)])
    finally:
        if kernel is not None:kernel.close()
        result['owned_jobs_joined']=True;result['elapsed_seconds']=time.monotonic()-start
        result['monitored_peak_host_bytes']=host_cap()
        result['source_hashes']={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__),HERE/'immutable_spine.py',HERE.parent/'field_service.py',HERE.parent/'service_native.py',
                      HERE.parent/'evaluator_program.py',HERE.parent/'constructor_ir.py',HERE.parent/'native_install.py']}
        dest=HERE/'is2-service-certificate-results.json';tmp=dest.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');tmp.replace(dest)
    print(json.dumps({k:result[k] for k in ('status','bundle_exit','checks','elapsed_seconds','owned_jobs_joined')}))
    raise SystemExit(result['bundle_exit'])


if __name__=='__main__':main()
