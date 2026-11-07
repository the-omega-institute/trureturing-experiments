"""Native guarded-star component and changed-image source certificate.

The component driver supplies public test arguments, never numerical answers.
Its first two calls derive partitions from the cache actually acquired through
OriginalPort and the native service. The final general-entry cut is a stopped
prefix, not a completed fallback. Production retains the full original inject.
"""
from __future__ import annotations
from collections import Counter
import argparse
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
from framing import Reader, REPLIES, mask_word, parse_arithmetic, v
from is2_chronology import acquisition
from native_install import BitKernel
from service_native import OriginalConnection
from session import OriginalPort, prepare_source

IMAGE = '326af82c315be2e1fb83cf6df5bcffdbb75888543473a3d3261f7d3d4a8d8d30'
CAPS = dict(seconds=600, cases=12, native_including_installation=50_000_000,
            monitored_resident_bytes=4*1024**3, reserve_cells=140_000)
SOURCE_PATHS = tuple(BASE/p for p in (
    'constructor_ir.py','evaluator_program.py','field_graph.py','field_service.py',
    'field_io.py','framing.py','lab.py','session.py','service_native.py',
    'native_engine.py','primitive_engine.py','paid_code.py','native_install.py','native_bits.c',
    'native-prerequisite/immutable_spine.py','native-prerequisite/is2_chronology.py',
    'native-prerequisite/is2_service_certificate.py',
    'native-prerequisite/is2_whole_session_demand.py',
    'native-prerequisite/guarded_star_component.py'))


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
               source_instruction_entries=0,push_entries=0,phase='source-construction')
    observer=Counter();before=None;source=port=conn=None;initial_source=None
    def observation_check():
        # Diagnostics share the original wall/resident ceilings, including finalization.
        assert time.monotonic()-start < CAPS['seconds'], 'diagnostic wall ceiling'
        monitor()
    def source_snapshot():
        objects=source.a.objects;ids={o:i for i,o in enumerate(objects)}
        digest=hashlib.sha256()
        for i,o in enumerate(objects):
            if i%4096==0:observation_check()
            fields=[(k,value if k=='bit' else ids[value] if value is not None else None)
                    for k,value in sorted(o.fields.items())]
            digest.update(im.canonical((o.kind,fields,o.writes)))
            observer['original_source_objects_scanned']+=1
            observer['original_source_fields_scanned']+=len(fields)
        return dict(material_sha256=digest.hexdigest(),objects=len(objects),
                    root_index=ids[source.root],events_sha256=hashlib.sha256(im.canonical(source.events)).hexdigest())
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
            state['component']=self
            self.handoffs=[];self.return_cuts=[];self.boundary_observer=self.observe
            # Descriptor ordinals are stable installation identities, never dispatch keys.
            self.diag_ids={s.forward:i for i,s in enumerate(self.sources)}
            self.diag_nodes={n:i for i,n in enumerate(program.nodes)}
            self.diag_aliases={}
            for name,node in program.entries.items():
                self.diag_aliases.setdefault(self.diag_nodes[node],[]).append(name)
            self.diag_fields={}
            for entry,idx in self.boundary.items():
                cursor=entry
                while True:
                    self.diag_fields[cursor]=idx;observer['native_site_mapping_visits']+=1
                    cursor=cursor.refs[2]
                    if cursor in self.boundary:break
                if idx%4096==0:observation_check()
            self.diag_instructions={}
            stops=sorted((self.diag_nodes[node],int(name.rsplit('.',1)[1]))
                         for name,node in program.entries.items() if name.startswith('is2.q.'))
            for j,(lo,q) in enumerate(stops):
                hi=stops[j+1][0] if j+1<len(stops) else self.diag_nodes[program.entries['is2.print']]
                for idx in range(lo,hi):self.diag_instructions[idx]=q
            self.diag_last_boundary=None;self.diag_last_source=None;self.diag_last_decision=None
            self.diag_pool_test=None;self.diag_previous=None;self.diag_attempt=None
            self.diag_offer=hashlib.sha256();self.diag_delivery=None
            self.diag_initial_native=self.native_identity()
        def diagnostic_ref(self,cell):
            if cell is None:return None
            if cell in self.diag_ids:return self.diag_ids[cell]
            return dict(unmapped_process_id=id(cell))
        def field_identity(self,idx):
            row=self.p.nodes[idx];q=self.diag_instructions.get(idx)
            return dict(index=idx,op=row.op,
                args=[self.diag_nodes[x] if x in self.diag_nodes else x for x in row.args],
                aliases=sorted(self.diag_aliases.get(idx,())),
                zero=self.diag_nodes[row.zero],one=self.diag_nodes[row.one],
                instruction=None if q is None else dict(index=q,row=self.fixed['code'][q],
                    aliases=sorted(self.test_labels.get('is2.q.'+str(q),()))))
        def pc_identity(self,pc):
            idx=self.diag_fields.get(pc)
            return dict(installed_source_ordinal=self.diagnostic_ref(pc),
                native_bits=None if pc is None else pc.bits,
                native_op=None if pc is None else im.OPS[pc.bits&31],
                native_literal=None if pc is None else pc.bits>>5,
                field=None if idx is None else self.field_identity(idx))
        def native_identity(self):
            source_digest=hashlib.sha256();installed_digest=hashlib.sha256()
            targets={s.block:s.forward for s in self.sources}
            for i,block in enumerate(self.native):
                if i%4096==0:observation_check()
                # Source references were relocated to retained Source descriptors by Direct.
                source_row=(block.bits,[None if x is None else self.diag_ids[x.forward] for x in block.refs])
                # Locate the target via the source directory, not by decoding an opcode.
                # The directory is in the same order as source_blocks.
                installed=targets[block]
                installed_row=(installed.bits,[self.diagnostic_ref(x) for x in installed.refs])
                source_digest.update(im.canonical(source_row));installed_digest.update(im.canonical(installed_row))
                observer['native_image_rows_scanned']+=1
            return dict(source_sha256=source_digest.hexdigest(),installed_sha256=installed_digest.hexdigest(),
                        matches=source_digest.digest()==installed_digest.digest(),rows=len(self.native))
        def pool_observation(self):
            cursor=self.regs['pool'].refs[0]
            # A suffix count is an observation of the original linked layout, not a traversal.
            head=self.diagnostic_ref(cursor);base=self.diag_ids[self.pool[0]]
            offset=head-base if isinstance(head,int) else None
            return dict(head=head,is_null=cursor is None,
                initial_suffix_cells=0 if cursor is None else
                    len(self.pool)-offset if offset is not None and 0<=offset<len(self.pool) else None,
                suffix_count_requires_unchanged_pool_links=True)
        def deliver(self,bit=None,end=False):
            observer['input_offers']+=1
            self.diag_offer.update(b'E' if end else str(int(bit or 0)).encode())
            self.diag_delivery=dict(offer=observer['input_offers'],end=end,bit=None if end else int(bit or 0))
            return sc.ServiceDirect.deliver(self,bit,end)
        def tick(self):
            state['native']=self.steps
            if self.steps%4096==0:check()
            assert self.steps+state['installation'] < CAPS['native_including_installation']
            n=self.pc;op=im.OPS[n.bits&31];steps=self.steps
            controls=self.events['control_write'];decision=None
            observer['native_attempts']+=1
            self.diag_attempt=(n,steps)
            if op in ('null','equal','test'):
                idx=self.diag_fields[n];row=self.p.nodes[idx]
                branch=int(self.rr is None) if op=='null' else int(self.rr is self.saved) if op=='equal' else self.bit
                decision=dict(pc=self.pc_identity(n),native_steps_before=steps,
                    rr=self.diagnostic_ref(self.rr),saved=self.diagnostic_ref(self.saved),bit=self.bit,
                    branch=branch,selected_successor=self.diagnostic_ref(n.refs[3 if branch else 2]),
                    pool=self.pool_observation() if row.op=='null' and row.args==('pool',) else None)
                if decision['pool'] is not None:
                    decision['pool']['tested_rr_is_pool_head']=self.rr is self.regs['pool'].refs[0]
                observer['conditional_decisions']+=1
            try:im.Direct.tick(self)
            finally:
                state['native']=self.steps
                if self.events['control_write']!=controls:
                    self.diag_previous=(n,self.pc,steps,self.steps)
                if decision is not None:
                    decision.update(executed=self.steps>steps,
                        control_written=self.events['control_write']!=controls,
                        actual_pc_after=self.diagnostic_ref(self.pc))
                    self.diag_last_decision=decision
                    if decision['pool'] is not None:
                        self.diag_pool_test=decision;observer['allocation_tests']+=1
        def drive(self):
            self.run(CAPS['native_including_installation']-state['installation'])
            assert self.done and not self.refused
        def observe(self,e,label):
            idx=self.boundary[self.pc];observer['field_boundaries']+=1
            self.diag_last_boundary=dict(native_steps=self.steps,pc=self.diagnostic_ref(self.pc),
                                         field_index=idx)
            if label.startswith('is2.q.'):
                self.diag_last_source=dict(self.diag_last_boundary)
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
        def terminal_receipts(self):
            out=[];ids=[];ordinals=[];seen=set();visits=0;status='complete';error=None
            cursor=self.ref('receipts').refs[1]
            while cursor is not None:
                if cursor in seen:status='receipt-cycle';break
                if visits>=CAPS['reserve_cells']:status='truncated';break
                if visits%4096==0:observation_check()
                seen.add(cursor);ids.append(id(cursor));ordinals.append(self.diagnostic_ref(cursor))
                visits+=1;bits=[];cell=cursor.refs[1];bit_seen=set()
                while cell is not None:
                    if cell in bit_seen:status='payload-cycle';break
                    if visits>=CAPS['reserve_cells']:status='truncated';break
                    if visits%4096==0:observation_check()
                    bit_seen.add(cell);visits+=1;bits.append(str(cell.bit(0)));cell=cell.refs[1]
                if status!='complete':break
                try:
                    rd=Reader(''.join(bits));word=rd.v().translate(str.maketrans('01','LR'))
                    code=rd.take(2);rd.end();out.append((word,REPLIES[int(code,2)]))
                except Exception as exc:status='malformed-payload';error=type(exc).__name__+': '+str(exc);break
                cursor=cursor.refs[3]
            observer['terminal_receipt_cell_visits']+=visits
            complete=status=='complete';payload_match=complete and before is not None and out==before[0]
            identity_match=complete and before is not None and ids==before[1]
            return dict(traversal_status=status,error=error,payloads=out,occurrence_process_ids=ids,
                occurrence_source_ordinals=ordinals,cell_visits=visits,cell_visit_limit=CAPS['reserve_cells'],
                comparison_status='not-reached' if before is None else 'incomplete' if not complete else
                    'complete-match' if payload_match and identity_match else 'comparison-false',
                payload_match=payload_match,occurrence_identity_match=identity_match,
                matches_original_calls=complete and out==[(w,r) for w,r,_,_ in port.calls])
    kernel=None;old_inst=im.Installer;e=None
    def timeout(*_):raise TimeoutError('one component batch wall ceiling')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(CAPS['seconds'])
    try:
        result['source_hashes_at_start']={str(p.relative_to(BASE)):sha(p) for p in SOURCE_PATHS}
        fixed,graph,labels=component_program(cases,static['syntactic_calls'])
        result['component_image']=dict(stack_sha256=hashlib.sha256(im.canonical(fixed)).hexdigest(),
            composed_sha256=graph.fingerprint(),stack_states=len(fixed['code']),
            stack_registers=len(fixed['registers']),field_nodes=len(graph.nodes),
            reference_registers=len(graph.refs),bit_registers=len(graph.bits),
            native_code_vertices=len(graph.nodes)+sum(sc.NATIVE_COST[n.op] for n in graph.nodes),
            cut='inject.general',production_fallback_retained=True)
        raw=(BASE/'constructor/positive-full.input.bits').read_text()
        n,theta,cache,cyl=parse_arithmetic(raw)
        rd=Reader(raw);rd.u();rd.rational()
        frames=[raw[:rd.pos]]+['01'+v(word.translate(str.maketrans('LR','01')))
                             for word in acquisition(n,cache)]+['00'+v(raw)]
        result['input_signature']=dict(raw_file_sha256=sha(BASE/'constructor/positive-full.input.bits'),
            raw_text_sha256=hashlib.sha256(raw.encode()).hexdigest(),raw_bits=len(raw),
            cases_sha256=hashlib.sha256(im.canonical(cases)).hexdigest(),cases=cases,
            frame_bits=list(map(len,frames)),frame_sha256=[hashlib.sha256(f.encode()).hexdigest() for f in frames],
            offered_stream_with_End_E_sha256=hashlib.sha256(''.join(f+'E' for f in frames).encode()).hexdigest(),
            schedule_sha256=hashlib.sha256(im.canonical(acquisition(n,cache))).hexdigest())
        source=prepare_source([dict(row).get(h,0) for row in cache for h in range(15)])
        port=OriginalPort(source);conn=OriginalConnection(port)
        initial_source=source_snapshot();result['original_source_initial']=initial_source
        state['phase']='installation'
        im.Installer=Installer;kernel=BitKernel();e=Component(graph,kernel,conn,labels,fixed)
        result['native_installations']=1;result['material']=dict(reserve_cells=CAPS['reserve_cells'],
            source_blocks=len(e.sources),native_census=e.census(),initializer_seed=e.inst.seed_census(),
            all_supplied_grants=True,all_installed_relocations_verified=True,
            scope='one retained component owner, actual source/connection and kernel; no whole-session offer')
        state['phase']='service-initialization';e.drive()
        state['phase']='original-prefix';rd=Reader(raw);rd.u();rd.rational();e.feed(raw[:rd.pos])
        for i,word in enumerate(acquisition(n,cache)):
            state['phase']='original-acquisition';state['acquisition_index']=i
            e.feed('01'+v(word.translate(str.maketrans('LR','01'))))
        paired=e.cache();attained=tuple(tuple((h,paired[15*j+h][1]) for h in range(15)
                    if paired[15*j+h][0]) for j in range(n))
        assert attained==cache
        receipts,ids=e.receipts();assert len(receipts)==74 and len(set(ids))==74
        result['acquisition']=dict(input_sha256=hashlib.sha256(raw.encode()).hexdigest(),
            original_calls=len(port.calls),receipt_occurrences=len(ids),cache_matches_original_frame=True,
            source_census=source.a.census(),connection_census=conn.a.census(),
            native_steps_before_component=e.steps,receipts_preserved_through_component=None,
            receipt_payloads=receipts,receipt_occurrence_process_ids=ids,
            strict_post_feed_comparison='not-reached')
        before=e.receipts()
        result['acquisition']['reconciliation_baseline']=dict(payloads=before[0],occurrence_process_ids=before[1])
        state['phase']='original-Arithmetic-frame';e.feed('00'+v(raw))
        after=e.receipts();matches=after==before
        result['acquisition']['strict_post_feed_comparison']='complete-match' if matches else 'comparison-false'
        result['acquisition']['receipts_preserved_through_component']=matches
        assert matches
        result['status']='bounded-component-complete-with-stopped-general-prefix';result['exit']=0
    except BaseException as exc:
        result['status']='unresolved-bounded-component';result['exit']=1
        result['error']=dict(type=type(exc).__name__,message=str(exc),
            locations=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name)
                       for x in traceback.extract_tb(exc.__traceback__)])
    finally:
        im.Installer=old_inst
        if e is None:e=state.get('component')
        if state['installer'] is not None:state['installation']=state['installer'].steps
        if e is not None:state['native']=e.steps
        result['executed_native_transitions']=state['native']
        result['executed_installation_transitions']=state['installation']
        result['aggregate_native_transitions']=state['native']+state['installation']
        result['general_entry_prefixes']=state['map_entries']
        done={x['index'] for x in result['cases']}
        result['unexecuted']=[dict(index=i,name=c['name']) for i,c in enumerate(cases) if i not in done]
        diagnostic=result['diagnostic']=dict(schema=1,active_phase=state['phase'],
            active_case=state['current'],active_case_name=None if state['current'] is None else cases[state['current']]['name'],
            case_start=state.get('case_start'),branches=list(state['branches']),
            acquisition_index=state.get('acquisition_index'),errors=[])
        def capture(name,read):
            try:diagnostic[name]=read()
            except BaseException as exc:
                diagnostic['errors'].append(dict(observation=name,type=type(exc).__name__,message=str(exc)))
                result['status']='unresolved-bounded-component';result['exit']=1
        inst=state['installer']
        if inst is not None:
            diagnostic['installer']=dict(done=inst.done,phase=inst.phase,steps=inst.steps,
                pc=None if inst.pc is None else dict(index=inst.nodes.index(inst.pc),
                    action=inst.pc.action,literal=inst.pc.literal),charges=dict(inst.charges))
        if e is not None:
            # Capture actual flags/PC before any bounded walks and before kernel.close().
            diagnostic['terminal']=dict(native_steps=e.steps,done=e.done,refused=e.refused,
                ended=e.ended,acked=e.acked,cpu_bit=e.bit,pc_process_id=id(e.pc),
                flags={k:b.bit(0) for k,b in e.bitregs.items()})
            result['native_events']=dict(e.events)
            capture('terminal_pc',lambda:e.pc_identity(e.pc))
            def operand_observations():
                cells=dict(cpu=e.cpu,pc=e.pc,rr=e.rr,saved=e.saved,address=e.address,
                    native_operand=e.pc.refs[0],pool=e.regs['pool'].refs[0])
                for name in ('port.cell','output.return'):
                    cells[name]=e.regs[name];cells[name+'.target']=e.regs[name].refs[0]
                return {k:dict(identity=e.diagnostic_ref(cell),grant=None if cell is None else cell.grant)
                        for k,cell in cells.items()}
            capture('terminal_operands',operand_observations)
            capture('last_attempt',lambda:None if e.diag_attempt is None else
                dict(pc=e.pc_identity(e.diag_attempt[0]),native_steps_before=e.diag_attempt[1]))
            capture('immediate_predecessor',lambda:None if e.diag_previous is None else
                dict(pc=e.pc_identity(e.diag_previous[0]),actual_successor=e.pc_identity(e.diag_previous[1]),
                     native_steps_before=e.diag_previous[2],native_steps_after=e.diag_previous[3]))
            capture('last_boundary',lambda:e.diag_last_boundary)
            capture('last_source_boundary',lambda:None if e.diag_last_source is None else
                dict(e.diag_last_source,field=e.field_identity(e.diag_last_source['field_index'])))
            capture('last_conditional_decision',lambda:e.diag_last_decision)
            capture('last_allocation_test',lambda:e.diag_pool_test)
            capture('input_delivery',lambda:dict(last_offer=e.diag_delivery,
                offered_stream_with_End_E_sha256=e.diag_offer.hexdigest()))
            capture('pool',lambda:remaining_pool(e,observation_check))
            if 'pool' in diagnostic:observer['terminal_pool_cell_visits']+=diagnostic['pool']['visited_cells']
            capture('receipts',e.terminal_receipts)
            capture('cuts',lambda:dict(handoffs=e.handoffs,return_cuts=e.return_cuts,
                cut_checks=e.cut_checks,source_instruction_entries=state['source_instruction_entries'],
                push_entries=state['push_entries'],general_entry_prefixes=state['map_entries']))
            def frozen_observation():
                mismatches=[];scanned=0;changed=0
                for cell,refs in e.frozen.items():
                    if scanned>=CAPS['reserve_cells']:break
                    if scanned%4096==0:observation_check()
                    scanned+=1
                    if tuple(cell.refs)!=refs:
                        changed+=1
                        if len(mismatches)<16:mismatches.append(e.diagnostic_ref(cell))
                observer['frozen_entries_scanned']+=scanned
                return dict(scanned=scanned,total=len(e.frozen),complete=scanned==len(e.frozen),
                    changed=changed,first_changed_source_ordinals=mismatches)
            capture('frozen',frozen_observation)
            capture('native_image_at_installation',lambda:e.diag_initial_native)
            capture('native_image_at_terminal',e.native_identity)
            capture('native_census_at_terminal',e.census)
            capture('observer_mapping_container_bytes',lambda:sum(sys.getsizeof(x) for x in
                (e.diag_ids,e.diag_nodes,e.diag_aliases,e.diag_fields,e.diag_instructions)))
            pool=diagnostic.get('pool',{})
            result['scoped_material_used']=dict(pool_cells_consumed=None if pool.get('status')!='complete' else
                CAPS['reserve_cells']-pool['visited_cells'],frozen_entries=len(e.frozen),output_bits=len(e.output))
        if source is not None:
            capture('original_source_at_terminal',source_snapshot)
            diagnostic['original_source_matches_initial']=None if initial_source is None or 'original_source_at_terminal' not in diagnostic else diagnostic['original_source_at_terminal']==initial_source
            diagnostic['original_connection']=dict(port_source_is_original=port.source is source,
                port_root_is_original=port.root is source.root,connection_port_is_original=conn.port is port,
                binding_root_is_original=conn.binding.fields['source'] is source.root,
                native_connection_is_original=None if e is None else e.connection is conn,
                calls=port.calls,port_observations=dict(port.observations),
                source_census=source.a.census(),connection_census=conn.a.census(),
                source_charges=dict(source.a.charges),connection_charges=dict(conn.a.charges))
        diagnostic['observer']=dict(counts=dict(observer),cost_scope='All observer construction, retention and final scans share batch wall time and monitored resident peaks; host work is additional to native transitions.',
            mapping_bytes_scope='Shallow containers only; not a backing bound.',H_reserved=None,W_backed=None,
            retries=0,budget_increases=0)
        if kernel is not None:kernel.close()
        capture('resident_ceiling_including_observer',observation_check)
        signal.alarm(0);result['elapsed_seconds']=time.monotonic()-start
        result['resident_monitor']=monitor(False);result['owned_jobs_joined']=True
    return result


def remaining_pool(e,check):
    cursor=e.regs['pool'].refs[0];seen=set();status='complete';links_match=True
    base=e.diag_ids[e.pool[0]]
    while cursor is not None:
        if cursor in seen:status='cycle';break
        if len(seen)>=CAPS['reserve_cells']:status='truncated';break
        if len(seen)%4096==0:check()
        seen.add(cursor);ordinal=e.diagnostic_ref(cursor);nxt=cursor.refs[0]
        expected=None if ordinal==base+len(e.pool)-1 else ordinal+1 if isinstance(ordinal,int) else 'unmapped'
        links_match &= e.diagnostic_ref(nxt)==expected
        cursor=nxt
    return dict(status=status,visited_cells=len(seen),visit_limit=CAPS['reserve_cells'],
                initial_suffix_links_match=links_match,stop=e.diagnostic_ref(cursor),head=e.pool_observation())


def fall(c,n):
    z=1
    for i in range(n):z*=max(0,c-i)
    return z


def main():
    assert sys.version_info>=(3,12) and __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--static',action='store_true')
    parser.add_argument('--result',type=Path,default=HERE/'guarded-star-component-diagnostics.json')
    args=parser.parse_args();target=args.result.resolve()
    if not args.static and (target==(HERE/'guarded-star-component-results.json').resolve() or target.exists()):
        parser.error('diagnostic result must be a new path distinct from the retained failed result')
    if not args.static and not target.parent.is_dir():parser.error('diagnostic result parent must exist')
    static,rows=full_source()
    result=dict(schema=1,identity='native-guarded-star-component/v1',static=static,
        complete_session_rows=rows,whole_session_prices_complete=False,
        H_reserved=None,W_backed=None,full_run_released=False,all_six_outcomes_complete=False)
    if args.static:
        cases=test_cases();p,g,labels=component_program(cases,static['syntactic_calls'])
        print(json.dumps(dict(static=static,rows=[{k:r[k] for k in ('case','K','execution_transitions_upper','chosen_reserve_payload_bytes')} for r in rows],
            component=dict(states=len(p['code']),field=len(g.nodes),refs=len(g.refs),bits=len(g.bits),
                           native=len(g.nodes)+sum(sc.NATIVE_COST[x.op] for x in g.nodes))),sort_keys=True))
        return
    batch_start=time.monotonic();result['component']=run_batch(static)
    result['source_hashes']={str(p.relative_to(BASE)):sha(p) for p in SOURCE_PATHS}
    temporary=target.with_name(target.name+'.tmp')
    payload=json.dumps(result,sort_keys=True,indent=2)+'\n'
    # Charge the result's serialization allocation/time as well as the live observer.
    component=result['component'];component['resident_monitor']=monitor(False)
    component['elapsed_through_result_serialization_seconds']=time.monotonic()-batch_start
    component['diagnostic']['serialized_bytes_before_accounting']=len(payload.encode())
    if (component['elapsed_through_result_serialization_seconds']>=CAPS['seconds'] or
        component['resident_monitor']['conservative_sum_of_peaks_bytes']>CAPS['monitored_resident_bytes']):
        component['status']='unresolved-bounded-component';component['exit']=1
        component['diagnostic']['errors'].append(dict(observation='result-serialization',type='ResourceCeiling',
            message='original batch wall or resident ceiling exceeded including diagnostics'))
    temporary.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');temporary.replace(target)
    print(json.dumps(dict(status=result['component']['status'],exit=result['component']['exit'],
        cases=result['component']['cases'],unexecuted=result['component']['unexecuted'],
        aggregate_native_transitions=result['component']['aggregate_native_transitions'],
        elapsed_seconds=result['component']['elapsed_seconds']),sort_keys=True),flush=True)
    raise SystemExit(result['component']['exit'])


if __name__=='__main__':main()
