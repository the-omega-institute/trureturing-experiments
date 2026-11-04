"""Native guarded-star component and changed-image source certificate.

The component driver supplies public test arguments, never numerical answers.
Its first two calls derive partitions from the cache actually acquired through
OriginalPort and the native service. The final general-entry cut is a stopped
prefix, not a completed fallback. Production retains the full original inject.
"""
from __future__ import annotations
from collections import Counter
import ast
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path[:0] = [str(HERE), str(BASE)]
import constructor_ir as ci
import evaluator_program as ep
import immutable_spine as im
import is2_service_certificate as sc
import is2_whole_session_demand as whole
from framing import Reader, mask_word, parse_arithmetic, v
from is2_chronology import acquisition
from native_install import BitKernel
from service_native import OriginalConnection
from session import OriginalPort, prepare_source

IMAGE = '326af82c315be2e1fb83cf6df5bcffdbb75888543473a3d3261f7d3d4a8d8d30'
CAPS = dict(seconds=600, cases=12, native_including_installation=50_000_000,
            monitored_resident_bytes=4*1024**3, reserve_cells=140_000)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def monitor(enforce=True):
    scale = 1 if sys.platform == 'darwin' else 1024
    own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*scale
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*scale
    if enforce:assert own+child <= CAPS['monitored_resident_bytes'], 'resident ceiling'
    return dict(self_peak_bytes=own, child_peak_bytes=child,
                conservative_sum_of_peaks_bytes=own+child)


class CountAssembler(ci.Assembler):
    def __init__(self):
        super().__init__(); self.own={}; self.edges={}
    def call(self, name, outs, *args):
        self.edges.setdefault(self.context, Counter())[name] += 1
        return super().call(name, outs, *args)
    def counted(self, name, body):
        def wrapped(a):
            start=len(self.g.code); body(a); rows=self.g.code[start:]
            self.own[name]=dict(instructions=len(rows)+(1 if name=='main' else 12),
                                pushes=sum(x[0]=='push' for x in rows))
        return wrapped
    def define(self, name, args, outs, body):
        return super().define(name,args,outs,self.counted(name,body))


def full_source():
    with im.shared_assembler():
        a=CountAssembler(); ci.primitives(a); ep.define_evaluator(a)
        raw=a.compile(a.counted('main',ep.main))
    fixed=im.prune(raw); linked=sc.compose(im.field_graph(fixed))
    assert linked.fingerprint()==IMAGE
    f=len(linked.nodes)+sum(sc.NATIVE_COST[x.op] for x in linked.nodes)
    r,b=len(linked.refs),len(linked.bits)
    static=dict(stack_sha256=hashlib.sha256(im.canonical(fixed)).hexdigest(),
                stack_states=len(fixed['code']), stack_registers=len(fixed['registers']),
                field_nodes=len(linked.nodes), reference_registers=r, bit_registers=b,
                native_code_vertices=f, composed_sha256=linked.fingerprint(),
                entries={k:linked.nodes.index(linked.entries[k]) for k in
                         ('is2.start','is2.print','is2.returned','service.return','evaluate.begin')},
                inject_result_register=fixed['registers']['inject.J'],
                own_sites=a.own, syntactic_calls={k:dict(v) for k,v in a.edges.items()})
    rows=[]
    for name in whole.PINS:
        old=whole.row(name,a.own,a.edges,len(sc.service_graph().nodes))
        ev=old['evaluator']; lam=ev['multipliers']; edges=ev['calls_after_guard']
        lam['inject']=2*old['N']+2
        # These children occur only inside the unreached general-map loop.
        edges['inject']['incident']=0; edges['inject']['odometer']=0
        memo={}
        def cost(fn):
            if fn in memo:return memo[fn]
            children=[(m,cost(g)) for g,m in edges.get(fn,{}).items() if m]
            memo[fn]={u:lam[fn]*(a.own[fn][u]+sum(m*c[u] for m,c in children))
                      for u in ('instructions','pushes')}
            return memo[fn]
        z=cost('main'); ev['body_costs']=memo
        ev['evaluator_instructions_upper']=z['instructions']
        ev['evaluator_pushes_upper']=z['pushes']
        ev['maps_per_inject_upper']=0
        k=5+old['L']+old['frontend']['exact_extractions']+z['pushes']
        n0=f+r+b+1; rr0=4*(f+1)+r; bb0=8+21*f+b
        nn=k+n0; rr=4*k+rr0; bb=8*k+bb0
        # ML-I's retained source/target/descriptors and unchanged initializer.
        census=dict(O=3*nn+1+98572, A=2*rr+5*nn+15+98908,
                    C=2*bb+2+2168397)
        common={u:census[u]+old['original_source'][u]+old['connection'][u]
                +dict(O=1281587,A=5224821,C=18724355)[u]
                -dict(O=98404,A=98404,C=22*98404)[u] for u in census}
        execution=128*(z['instructions']+old['L']+2*ev['Y']
                       +len(fixed['registers'])+64)+5*old['frontend']['field_actions_upper']
        old.update(K=k,K_bit_length=k.bit_length(),source_blocks=nn,
                   native_census=census,common_with_B=common,
                   installation_transitions_upper=128*nn+16*rr+4*bb+64,
                   execution_transitions_upper=execution,
                   chosen_reserve_payload_bytes=500*k,
                   frozen_entries_upper=2*old['L']+z['pushes']+5,
                   observer_name_visits_upper=r*(z['instructions']+2),
                   observer_stack_head_visits_upper=len(fixed['registers'])*(z['instructions']+2),
                   chosen_K_fits_Py_ssize_t=k<=sys.maxsize,
                   chosen_K_pointer_array_bytes_fit_Py_ssize_t=8*k<=sys.maxsize,
                   chosen_common_objects_fit_64bit_names=common['O']<2**64)
        rows.append(old)
    static['layout_constants']=dict(n=n0,refs=rr0,bits=bb0,
        O=3*n0+1+98572,A=2*rr0+5*n0+15+98908,C=2*bb0+2+2168397)
    return static,rows


def test_cases():
    # Each pair is (p,q); repeated incidences are intentionally retained.
    out=[dict(name='acquired-overlap-G'+str(i),source=i,acquired=True) for i in (0,1)]
    for name,p,q in [('left-two',(0,0,0,0),(0,1,0,1)),
                     ('right-two',(0,1,0,1),(0,0,0,0))]:
        out += [dict(name=name+'-G'+str(i),p=p,q=q,source=i) for i in (0,1)]
    out += [dict(name='left-degree-27',p=(0,)*27,q=tuple(range(27)),source=0),
            dict(name='left-overdegree-28',p=(0,)*28,q=tuple(range(28)),source=1),
            dict(name='right-degree-8',p=tuple(range(8)),q=(0,)*8,source=1),
            dict(name='right-overdegree-9',p=tuple(range(9)),q=(0,)*9,source=0),
            dict(name='conjunction-conflict',conflict=True,source=0),
            dict(name='conjunction-nonstar-prefix',fallback=True,source=1)]
    assert len(out)==CAPS['cases']
    return out


def component_program(cases, call_edges):
    with im.shared_assembler():
        a=ci.Assembler(); ci.primitives(a); ep.define_evaluator(a)
        # Compile only the component's actual routine closure. The original
        # assembler otherwise emits return alternatives from unused bodies.
        needed=set(); todo=['readu','readcache','side','blocks','conjoin','inject']
        while todo:
            name=todo.pop()
            if name in needed:continue
            needed.add(name)
            todo.extend(g for g,m in call_edges.get(name,{}).items()
                        if m and not (name=='inject' and g in ('incident','odometer')))
        a.functions={name:value for name,value in a.functions.items() if name in needed}
        ins,outs,body=a.functions['inject']
        class PrefixCut(Exception):pass
        def prefix(a):
            original_mark=a.g.mark
            def mark(name):
                original_mark(name)
                if name=='inject.general':raise PrefixCut
            a.g.mark=mark
            try:body(a)
            except PrefixCut:a.g.emit('halt',False)
            finally:a.g.mark=original_mark
        a.functions['inject']=(ins,outs,prefix)
        def main(a):
            v0=a.var
            a.call('readu',[v0('actual.n'),'input'],'input')
            ep.bit(a,'input',v0('prior.sign'))
            a.call('readu',[v0('prior.a'),'input'],'input')
            a.call('readu',[v0('prior.b'),'input'],'input')
            a.call('readcache',[v0('actual.cache'),'input'],'input',v0('actual.n'))
            for i,c in enumerate(cases):
                a.g.mark('component.case.'+str(i))
                if c.get('acquired'):
                    a.set(v0('n'),v0('actual.n')); a.set(v0('data'),v0('actual.cache'))
                    a.call('side',[v0('p'),v0('fullp')],v0('data'),v0('n'),7,0)
                    a.call('side',[v0('q'),v0('fullq')],v0('data'),v0('n'),8,7)
                elif c.get('conflict'):
                    a.set(v0('n'),v0('actual.n'))
                    mask=tuple(((0,1),) if j==0 else () for j in range(4))
                    a.raw(v0('fixture'),map(int,mask_word(mask)))
                    a.call('readcache',[v0('extra'),v0('fixture')],v0('fixture'),v0('n'))
                    a.call('conjoin',[v0('data'),v0('valid')],v0('actual.cache'),v0('extra'),v0('n'))
                elif c.get('fallback'):
                    a.set(v0('n'),4)
                    left=(0,0,1,1);right=(0,1,0,1)
                    base=tuple(tuple((h,0 if h!=14 else right[j]) for h in range(15)
                                     if h!=6 or j<2) for j in range(4))
                    extra=tuple(((6,left[j]),) if j>=2 else () for j in range(4))
                    for name,data in [('base',base),('extra',extra)]:
                        a.raw(v0('fixture'),map(int,mask_word(data)))
                        a.call('readcache',[v0(name),v0('fixture')],v0('fixture'),v0('n'))
                    a.call('conjoin',[v0('data'),v0('valid')],v0('base'),v0('extra'),v0('n'))
                    a.call('side',[v0('p'),v0('fullp')],v0('data'),v0('n'),7,0)
                    a.call('side',[v0('q'),v0('fullq')],v0('data'),v0('n'),8,7)
                else:
                    a.set(v0('n'),len(c['p']))
                    for name in ('p','q'):a.raw(v0(name),ep.encode_numbers(c[name]))
                if c.get('conflict'):
                    a.set(v0('J'),0)
                else:
                    a.call('blocks',[v0('r')],v0('p'),v0('n'))
                    a.call('blocks',[v0('s')],v0('q'),v0('n'))
                    a.call('inject',[v0('J')],v0('p'),v0('q'),v0('n'),v0('r'),v0('s'),c['source'])
                a.g.mark('component.return.'+str(i))
            a.g.copy(v0('J'),'input')
        raw=a.compile(main)
    # Literal source prefix: no production successor beyond this cut is claimed.
    cut=raw['labels']['inject.general']; assert raw['code'][cut]==('halt',False)
    seen,todo=set(),[raw['start']]
    while todo:
        q=todo.pop()
        if q not in seen:seen.add(q);todo.extend(im.successors(raw['code'][q]))
    ids={q:i for i,q in enumerate(sorted(seen))}
    labels={name:ids[q] for name,q in raw['labels'].items() if q in ids}
    fixed=im.prune(raw);graph=sc.compose(im.field_graph(fixed))
    return fixed,graph,labels


def natural(bits):
    return sum(int(x)<<i for i,x in enumerate(bits))


def run_batch(static):
    start=time.monotonic(); cases=test_cases(); result=dict(cases=[],unexecuted=[],
        native_installations=0,full_evaluator_executions=0,full_chronologies=0,
        caps=CAPS,batch_count=1,acquisition=None)
    state=dict(installation=0,native=0,installer=None,current=None,branches=[],map_entries=0,
               source_instruction_entries=0,push_entries=0)
    def check():
        assert state['installation']+state['native'] < CAPS['native_including_installation'], 'aggregate native cap'
        monitor()
    class Installer(im.Installer):
        def run(self,checkpoint=None):
            state['installer']=self
            while not self.done:
                check(); remaining=CAPS['native_including_installation']-state['native']-self.steps
                assert remaining>2,'installation native cap'
                self.tick(min(10000,remaining-2));state['installation']=self.steps
            return None
    class Component(sc.ServiceDirect):
        def __init__(self,program,kernel,connection,labels,fixed):
            self.test_labels={}
            for name,idx in labels.items():self.test_labels.setdefault('is2.q.'+str(idx),[]).append(name)
            self.fixed=fixed;self.install_bound=None
            im.Direct.__init__(self,program,'',CAPS['reserve_cells'],kernel,
                               connection=connection,service_return=True)
            self.handoffs=[];self.return_cuts=[];self.boundary_observer=self.observe
        def tick(self):
            state['native']=self.steps
            if self.steps%4096==0:check()
            assert self.steps+state['installation'] < CAPS['native_including_installation']
            im.Direct.tick(self);state['native']=self.steps
        def drive(self):
            self.run(CAPS['native_including_installation']-state['installation'])
            assert self.done and not self.refused
        def observe(self,e,label):
            sc.ServiceDirect.observe(self,e,label)
            if label.startswith('is2.q.'):
                idx=int(label.rsplit('.',1)[1]);state['source_instruction_entries']+=1
                state['push_entries']+=self.fixed['code'][idx][0]=='push'
            for name in self.test_labels.get(label,()):
                if name.startswith('component.case.'):
                    i=int(name.rsplit('.',1)[1]);state['current']=i
                    state['case_start']=(self.steps,state['source_instruction_entries'],state['push_entries'])
                    state['branches']=[]
                elif name in ('inject.star.left','inject.star.right'):
                    state['branches'].append(name)
                elif name=='inject.general':
                    state['map_entries']+=1
                    result['cases'].append(dict(index=state['current'],name=cases[state['current']]['name'],
                        status='stopped-before-general-map-initialization',fallback_completed=False,
                        native_steps=self.steps-state['case_start'][0],branches=list(state['branches']),
                        r=self.value_word('main.r'),s=self.value_word('main.s')))
                elif name.startswith('component.return.'):
                    i=int(name.rsplit('.',1)[1]);c=cases[i];got=self.value_word('main.J')
                    r,s=(1,1) if c.get('acquired') else (0,0) if c.get('conflict') else (max(c['p'])+1,max(c['q'])+1)
                    expected=0 if c.get('conflict') else 72*fall(27,s) if r==1 else 243*fall(8,r)
                    assert got==expected,(c['name'],got,expected)
                    if c.get('conflict'):assert self.value_word('main.valid')==0
                    else:assert (self.value_word('main.r'),self.value_word('main.s'))==(r,s)
                    result['cases'].append(dict(index=i,name=c['name'],status='computed-native-return' if not c.get('conflict') else 'native-conjunction-rejected',
                        J=got,r=r,s=s,source=c['source'],native_steps=self.steps-state['case_start'][0],
                        source_instruction_entries=state['source_instruction_entries']-state['case_start'][1],
                        push_entries=state['push_entries']-state['case_start'][2],branches=list(state['branches']),
                        general_map_initializations=0,expected_used_only_after_return=True))
        def value_word(self,name):
            return natural(self.word('is2.stack.'+str(self.fixed['registers'][name])))
    kernel=None;old_inst=im.Installer;e=None
    def timeout(*_):raise TimeoutError('one component batch wall ceiling')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(CAPS['seconds'])
    try:
        fixed,graph,labels=component_program(cases,static['syntactic_calls'])
        result['component_image']=dict(stack_sha256=hashlib.sha256(im.canonical(fixed)).hexdigest(),
            composed_sha256=graph.fingerprint(),stack_states=len(fixed['code']),
            stack_registers=len(fixed['registers']),field_nodes=len(graph.nodes),
            reference_registers=len(graph.refs),bit_registers=len(graph.bits),
            native_code_vertices=len(graph.nodes)+sum(sc.NATIVE_COST[n.op] for n in graph.nodes),
            cut='inject.general',production_fallback_retained=True)
        raw=(BASE/'constructor/positive-full.input.bits').read_text()
        n,theta,cache,cyl=parse_arithmetic(raw)
        source=prepare_source([dict(row).get(h,0) for row in cache for h in range(15)])
        port=OriginalPort(source);conn=OriginalConnection(port)
        im.Installer=Installer;kernel=BitKernel();e=Component(graph,kernel,conn,labels,fixed)
        result['native_installations']=1;result['material']=dict(reserve_cells=CAPS['reserve_cells'],
            source_blocks=len(e.sources),native_census=e.census(),initializer_seed=e.inst.seed_census(),
            all_supplied_grants=True,all_installed_relocations_verified=True,
            scope='one retained component owner, actual source/connection and kernel; no whole-session offer')
        e.drive();rd=Reader(raw);rd.u();rd.rational();e.feed(raw[:rd.pos])
        for word in acquisition(n,cache):e.feed('01'+v(word.translate(str.maketrans('LR','01'))))
        paired=e.cache();attained=tuple(tuple((h,paired[15*j+h][1]) for h in range(15)
                    if paired[15*j+h][0]) for j in range(n))
        assert attained==cache
        receipts,ids=e.receipts();assert len(receipts)==74 and len(set(ids))==74
        result['acquisition']=dict(input_sha256=hashlib.sha256(raw.encode()).hexdigest(),
            original_calls=len(port.calls),receipt_occurrences=len(ids),cache_matches_original_frame=True,
            source_census=source.a.census(),connection_census=conn.a.census(),
            native_steps_before_component=e.steps,receipts_preserved_through_component=False)
        before=e.receipts();e.feed('00'+v(raw))
        assert e.receipts()==before;result['acquisition']['receipts_preserved_through_component']=True
        result['status']='bounded-component-complete-with-stopped-general-prefix';result['exit']=0
    except BaseException as exc:
        result['status']='unresolved-bounded-component';result['exit']=1
        result['error']=dict(type=type(exc).__name__,message=str(exc),
            locations=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name)
                       for x in traceback.extract_tb(exc.__traceback__)])
    finally:
        im.Installer=old_inst
        if state['installer'] is not None:state['installation']=state['installer'].steps
        if e is not None:state['native']=e.steps
        result['executed_native_transitions']=state['native']
        result['executed_installation_transitions']=state['installation']
        result['aggregate_native_transitions']=state['native']+state['installation']
        result['general_entry_prefixes']=state['map_entries']
        done={x['index'] for x in result['cases']}
        result['unexecuted']=[dict(index=i,name=c['name']) for i,c in enumerate(cases) if i not in done]
        if e is not None:
            result['native_events']=dict(e.events);result['scoped_material_used']=dict(
                pool_cells_consumed=CAPS['reserve_cells']-remaining_pool(e),
                frozen_entries=len(e.frozen),output_bits=len(e.output))
        if kernel is not None:kernel.close()
        signal.alarm(0);result['elapsed_seconds']=time.monotonic()-start
        result['resident_monitor']=monitor(False);result['owned_jobs_joined']=True
    return result


def remaining_pool(e):
    cursor=e.regs['pool'].refs[0];n=0
    while cursor is not None:n+=1;cursor=cursor.refs[0]
    return n


def fall(c,n):
    z=1
    for i in range(n):z*=max(0,c-i)
    return z


def main():
    assert sys.version_info>=(3,12) and __debug__
    static,rows=full_source()
    result=dict(schema=1,identity='native-guarded-star-component/v1',static=static,
        complete_session_rows=rows,whole_session_prices_complete=False,
        H_reserved=None,W_backed=None,full_run_released=False,all_six_outcomes_complete=False)
    if '--static' in sys.argv:
        cases=test_cases();p,g,labels=component_program(cases,static['syntactic_calls'])
        print(json.dumps(dict(static=static,rows=[{k:r[k] for k in ('case','K','execution_transitions_upper','chosen_reserve_payload_bytes')} for r in rows],
            component=dict(states=len(p['code']),field=len(g.nodes),refs=len(g.refs),bits=len(g.bits),
                           native=len(g.nodes)+sum(sc.NATIVE_COST[x.op] for x in g.nodes))),sort_keys=True))
        return
    result['component']=run_batch(static)
    result['source_hashes']={str(p.relative_to(BASE)):sha(p) for p in
        (BASE/'evaluator_program.py',HERE/'is2_chronology.py',Path(__file__),
         BASE/'constructor_ir.py',HERE/'immutable_spine.py',HERE/'is2_service_certificate.py',
         BASE/'native_install.py',BASE/'native_bits.c',BASE/'field_service.py',BASE/'session.py')}
    target=HERE/'guarded-star-component-results.json';temporary=target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');temporary.replace(target)
    print(json.dumps(dict(status=result['component']['status'],exit=result['component']['exit'],
        cases=result['component']['cases'],unexecuted=result['component']['unexecuted'],
        aggregate_native_transitions=result['component']['aggregate_native_transitions'],
        elapsed_seconds=result['component']['elapsed_seconds']),sort_keys=True),flush=True)
    raise SystemExit(result['component']['exit'])


if __name__=='__main__':main()
