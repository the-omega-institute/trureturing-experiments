"""IS2 source-cost/layout discrimination, without evaluator execution.

Derive costs from the assembler's actual acyclic call graph, public frames and
bounded source loops. Install the fixed graph with ZERO execution reserve and
never tick it. Retain an independently installed mutable B owner at prefix zero.
The byte census is a scoped reachable Python payload census, not RSS, an
allocator upper bound, a reserve, or a completion promise.
"""
from __future__ import annotations

import ast
from collections import Counter
import ctypes as C
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import struct
import sys
import sysconfig
import time
import types

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path[:0] = [str(HERE), str(BASE)]
import constructor_ir as ci
import evaluator_program as ep
import immutable_spine as im
import is2_service_certificate as old
import native_engine as ne
from field_graph import b_graph, install_b
from native_install import BitKernel, Block, Source, Installer, InstallNode, Wire, Init
from service_native import OriginalConnection
from session import OriginalPort, prepare_source
from framing import parse_arithmetic, Reader, v

assert sys.version_info >= (3, 12) and __debug__
CAPS = dict(wall_seconds=600, installation_events=90_000_000,
            source_blocks=650_000, modeled_objects=2_100_000,
            monitored_RSS_bytes=2 * 1024**3, B_native_steps=2_000_000,
            full_evaluator_ticks=0, evaluator_reserve=0, B_reserve=128,
            source_bits=135, graph_nodes=70_000)
CASES = ('positive-full', 'positive-partial', 'null-nine-distinct-left')
PINS = {
    'positive-full': (4, 387, 74, 2176, '37dfee71eda5a8fbee3d373aab446982f15976785ef965ec2e4a1b6ea98c5402'),
    'positive-partial': (4, 386, 73, 2115, '1a4ae049cf81b1f434398aa05cb8eabcf4d8689f7873f8d75f9ae1e870c34764'),
    'null-nine-distinct-left': (9, 839, 149, 9676, 'f5518d034eef5bd00c87ff9d43c5ca878ef9d012f47e04bedcbe45163ee3e29d')}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rss():
    n = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return n if sys.platform == 'darwin' else 1024*n


def checkpoint(ins):
    assert ins.steps <= CAPS['installation_events'], 'installation event cap'
    assert len(ins.sources) <= CAPS['source_blocks'], 'source block cap'
    if ins.steps % 4096 < 128:
        assert rss() <= CAPS['monitored_RSS_bytes'], 'RSS containment cap'


class BoundedInstaller(Installer):
    def run(self, checkpoint=None):
        assert len(self.sources) <= CAPS['source_blocks']
        while not self.done:
            self.tick(1 << 26)
            globals()['checkpoint'](self)
            if checkpoint:
                checkpoint(self)


class CallAssembler(ci.Assembler):
    def __init__(self):
        super().__init__()
        self.edges = {}

    def call(self, name, outs, *args):
        self.edges.setdefault(self.context, Counter())[name] += 1
        return super().call(name, outs, *args)


def fixed_source():
    with im.shared_assembler():
        a = CallAssembler()
        ci.primitives(a)
        primitive_names = tuple(a.functions)
        ep.define_evaluator(a)
        raw = a.compile(ep.main)
    fixed = im.prune(raw)
    field = im.field_graph(fixed)
    linked = old.compose(field)
    assert len(linked.nodes) <= CAPS['graph_nodes']
    assert (len(fixed['code']), len(fixed['registers']), len(linked.nodes),
            len(linked.refs), len(linked.bits)) == (10983, 360, 62069, 604, 57)
    assert sha(BASE/'constructor/positive-full.input.bits') == PINS['positive-full'][4]
    assert hashlib.sha256(im.canonical(fixed)).hexdigest() == '7257d4ad4b5730860ad811d41c3057ea6accdf1c413dbf7af75be206a5b24595'
    assert linked.fingerprint() == 'db4fb22177743a3aafd286c5431b42d09460ab47cfe4ca488d40dd1f28600b4a'
    # Acyclic routine invocation, distinct from loops within a routine.
    active, memo = set(), {}
    def depth(f):
        assert f not in active, ('recursive call', f)
        if f in memo:
            return memo[f]
        active.add(f)
        d = 1 + max((depth(g) for g in a.edges.get(f, {})), default=0)
        active.remove(f)
        memo[f] = d
        return d
    for f in ('main', *a.functions):
        depth(f)
    return a, primitive_names, raw, fixed, field, linked, memo


def evaluator_bound(n, length, t, a, primitives, raw):
    # IS.9--10, including unreduced represented rational payloads.
    m = 64*(n+length+t+2)
    y = 64*(t+2)*(m+1)
    d0 = length+y+(n+1)**2+(16*n+3)*(1 << min(n, 243))+m+1024
    b = 16**(len(a.functions)+2)*(d0+1)
    partitions = n**n  # n >= 4. Bell(n) <= n^n.
    maps = 36**min(n, 72)*27**min(n, 243)
    subsets = 1 << min(n, 243)
    # Each local state, with calls treated as opaque edges, has at most this
    # many visits per invocation. The ordinary proof lists the actual loops.
    multipliers = {f: (2*b+3)**2 for f in primitives}
    multipliers.update({
        'main': t+2,
        'readcache': 15*n+2, 'conjoin': 15*n+2, 'blocks': n+2,
        'nextpart': (n+2)**2, 'side': 16*(n+2)**2,
        'countlabels': 257*(subsets+1)*(n+257)**3,
        'falling': n+2, 'odometer': 2*n+2, 'incident': 1,
        'inject': (maps+1)*(2*n+2),
        'normalizer': (partitions+1)**2*(n+2),
        'mixture': 1, 'addword': 1})
    assert set(multipliers) == {'main', *a.functions}
    # Use the complete emitted source, including unreachable dispatch arms.
    q0 = len(raw['code'])+1
    memo = {}
    def cost(f):
        if f not in memo:
            memo[f] = q0*multipliers[f]*(1+sum(k*cost(g) for g,k in a.edges.get(f, {}).items()))
        return memo[f]
    steps = cost('main')
    return dict(m=m, Y=y, B=b, partitions_upper=partitions,
                maps_upper=maps, subsets_upper=subsets, source_Q=q0,
                routine_local_multipliers=multipliers, routine_costs=memo,
                stack_step_upper=steps, stack_push_upper=steps)


def frontend_bound(i, s, n, t, q, length, y, field_count):
    # A bound on every scan word, count value and total retained output.
    j = i+s+15*n+t+q+length+y+1024
    f = field_count+1
    c = {'take': f}
    c['readv'] = f*(j+1)*(1+c['take'])
    c['readu'] = f*(1+c['readv'])
    c['unary'] = f*(j+1)**2
    c['shape'] = f*(j+1)*(1+c['unary'])
    c['prior'] = f*(j+1)*(1+2*c['readu'])
    c['mask'] = f*(j+1)*(1+c['take'])
    c['cylinders'] = f*(j+1)*(1+c['readu']+c['unary']+c['readv']+c['mask']+c['take'])
    c['arithmetic'] = f*(1+c['readu']+c['shape']+c['prior']+c['mask']+c['cylinders'])
    c['classify'] = f*(j+1)
    # Includes receive, copy/cat/wrap, prediction construction if used, print,
    # request entry, End/Ack/Return, and the final suspended closed cut.
    top = f*(j+1)**2*(1+sum(c.values()))
    total = (i+q+3)*top
    return dict(J=j, F0=f-1, routines=c, top=top, field_action_upper=total,
                extracted_cells_upper=total)


def source_census(k):
    d = (5*k*k+13*k+6)//2
    vertices, headers = 8*k+8, 3*k+3
    return dict(O=vertices+headers+d, A=16+4*vertices+2*headers+d, C=d)


def port_census(q,s):
    f = int(q>0)
    return dict(O=1+2*s+5*q+f, A=10+2*s+7*q+2*f, C=2*s+3*q)


def row(name, a, primitives, raw, linked, seed):
    bits = (BASE/'constructor'/(name+'.input.bits')).read_text()
    n, theta, cache, cyl = parse_arithmetic(bits)
    p = old.static_demand(bits,10983,360,34)
    assert (n,len(bits),p['acquisition_calls'],p['acquisition_operand_bits'],sha(BASE/'constructor'/(name+'.input.bits'))) == PINS[name]
    eb = evaluator_bound(n,len(bits),len(cyl),a,primitives,raw)
    fb = frontend_bound(p['total_input_bits'],p['complete_operand_bits'],n,len(cyl),
                        p['complete_calls'],len(bits),eb['Y'],len(old.service_graph().nodes))
    k = 5+len(bits)+eb['stack_push_upper']+fb['extracted_cells_upper']
    f = len(linked.nodes)+sum(old.NATIVE_COST[n.op] for n in linked.nodes)
    assert f == 225718
    nn = k+f+len(linked.refs)+len(linked.bits)+1
    rr = 4*(k+f+1)+len(linked.refs)
    cc = 8*(k+1)+21*f+len(linked.bits)
    native = dict(O=3*nn+1+seed['objects'],A=2*rr+5*nn+15+seed['refs'],C=2*cc+2+seed['bits'])
    src = source_census(15*n)
    conn = port_census(p['complete_calls'],p['complete_operand_bits'])
    install = 128*nn+16*rr+4*cc+64
    vc = install+f*(eb['stack_step_upper']+fb['field_action_upper']+2*len(bits)+2*eb['Y']+360+64)
    oo = 7*p['complete_calls']+2*eb['Y']+6
    # Host data facts, independent of a protected offer or allocator overhead:
    # each source/target reserve Block has the actual four-pointer refs list,
    # one-byte data buffer and slotted Block; the descriptor is distinct.
    sizes = host_shapes()
    reserve_payload_slope = 2*(sizes['Block']+sizes['list4']+sizes['bytearray1'])+sizes['Source']
    return dict(case=name,N=n,L=len(bits),t=len(cyl),sha256=PINS[name][4],
                acquisition_reads=p['acquisition_calls'],acquisition_operand_bits=p['acquisition_operand_bits'],
                reads=p['complete_calls'],operand_bits=p['complete_operand_bits'],
                I=p['total_input_bits'],represented_prior=p['represented_prior'],
                evaluator=eb,frontend=fb,K_upper=k,source_blocks_upper=nn,
                native_census_upper=native,original_source_exact=src,connection_exact=conn,
                installation_native_upper=install,composed_native_upper=vc,
                output_bits_upper=oo,
                provider_vertices_upper=p['complete_operand_bits']+p['complete_calls'],
                provider_navigation_upper=p['complete_operand_bits'],
                input_latch_writes=3*(p['total_input_bits']+p['complete_calls']+3),
                receipt_payload_bits=2*p['complete_operand_bits']+3*p['complete_calls'],
                host_reserve_owned_payload_floor_at_chosen_K=reserve_payload_slope*k,
                reserve_payload_bytes_per_K=reserve_payload_slope,
                chosen_K_exceeds_Py_ssize_t=(k>sys.maxsize),
                actual_H_reserved=None,actual_W_backed=None,
                comparison='no same-scope protected byte/work offer; upper demand not a lower bound',
                full_session_executed=False)


def host_shapes():
    b = Block([None]*4,bytearray(1),8)
    s = Source(b)
    return {name:sys.getsizeof(obj) for name,obj in {
        'Block':b,'Source':s,'list0':[],'list4':[None]*4,'list1':[None],
        'tuple4':(None,)*4,'dict0':{},'set0':set(),
        'bytearray0':bytearray(),'bytearray1':bytearray(1),
        'bytearray3':bytearray(3),'int0':0,'int_30bit':(1<<30)-1,
        'int_31bit':1<<30,'str0':'','str1':'0'}.items()}


def payload_census(roots):
    """Identity-union of owned payload; opaque runtime leaves are separated.

    Traverses the actual object dicts/slots and ctypes keepalive maps. Class,
    function, module and code closures are deliberately NOT recursively claimed.
    Work buffers of this measurement are outside its roots and priced apart.
    """
    opaque_types = (type, types.ModuleType, types.FunctionType, types.CodeType,
                    types.MethodType, types.BuiltinFunctionType, types.FrameType)
    seen, todo, counts, bytecounts = set(),list(roots),Counter(),Counter()
    opaque = Counter()
    ctypes_owned_buffer_allowance = 0
    while todo:
        x = todo.pop()
        xid = id(x)
        if xid in seen:
            continue
        seen.add(xid)
        name = type(x).__module__+'.'+type(x).__name__
        counts[name] += 1
        size = sys.getsizeof(x)
        bytecounts[name] += size
        if isinstance(x,opaque_types):
            opaque[name] += 1
            continue
        if isinstance(x,dict):
            todo.extend(x.keys());todo.extend(x.values())
        elif isinstance(x,(list,tuple,set,frozenset)):
            todo.extend(x)
        elif isinstance(x,(str,bytes,bytearray,int,float,type(None),bool)):
            pass
        else:
            d = getattr(x,'__dict__',None)
            if d is not None:
                todo.append(d)
            for cls in type(x).__mro__:
                slots = cls.__dict__.get('__slots__',())
                if isinstance(slots,str):slots=(slots,)
                for slot in slots:
                    if slot not in ('__dict__','__weakref__') and hasattr(x,slot):
                        todo.append(getattr(x,slot))
            # ctypes native backing and its separately retained keepalive maps.
            if isinstance(x,(C.Array,C.Structure,C._Pointer)):
                keep = getattr(x,'_objects',None)
                if keep is not None:todo.append(keep)
                if getattr(x,'_b_needsfree_',False):
                    # getsizeof does not expose separately allocated ctypes
                    # buffers. Small inline buffers can be counted twice here;
                    # this coordinate is an allowance, not an exact heap size.
                    ctypes_owned_buffer_allowance += C.sizeof(x)
    return dict(objects=len(seen),shallow_payload_bytes=sum(bytecounts.values()),
                ctypes_owned_buffer_allowance=ctypes_owned_buffer_allowance,
                payload_plus_ctypes_allowance=sum(bytecounts.values())+ctypes_owned_buffer_allowance,
                types=dict(counts),type_bytes=dict(bytecounts),opaque_runtime_leaves=dict(opaque),
                excludes='allocator arenas/fragmentation, module/code closure, stack, mappings, OS/runtime/C compiler and census work buffers')


def actual_reference_union(direct,b,source,connection):
    descriptors = direct.inst.sources+b.installation.sources
    blocks = [s.block for s in descriptors]+[s.forward for s in descriptors]
    assert all(x is not None for x in blocks)
    objects = descriptors+blocks+direct.inst.nodes+b.installation.nodes
    # One shared immutable wiring, two different mutable Init headers.
    assert direct.inst.kernel.wires is b.kernel.wires
    assert direct.inst.kernel.state is not b.kernel.state
    assert direct.inst.kernel.end is b.kernel.end
    assert direct.inst.kernel.tmp is not None and b.kernel.tmp is None
    wiring_checks=0
    for wires in direct.inst.kernel.wires.values():
        for i,wire in enumerate(wires):
            expected=C.addressof(wires[i+1]) if i+1<len(wires) else C.addressof(direct.inst.kernel.end)
            assert C.cast(wire.next,C.c_void_p).value==expected and wire.grant==1
            wiring_checks+=1
    installer_control_checks=0
    for ins in (direct.inst,b.installation):
        own={id(x) for x in ins.nodes}
        for x in ins.nodes:
            assert x.zero is None or id(x.zero) in own
            assert x.one is None or id(x.one) in own
            if x.action=='bind-cursor':
                from native_install import LAYOUTS
                width=LAYOUTS[x.literal][1]
                assert C.cast(x.binding,C.c_void_p).value==C.addressof(ins.kernel.wires[width][0])
            installer_control_checks+=1
    domain = {id(x) for x in descriptors+blocks}
    refs=0
    for block in blocks:
        for target in (*block.refs,block.next):
            if target is not None:
                assert id(target) in domain
            refs += 1
    for s in descriptors:
        assert id(s.block) in domain and id(s.forward) in domain
        assert s.next is None or id(s.next) in domain
        refs += 3
    arenaobjects = source.a.objects+connection.a.objects
    adomain = {id(x) for x in arenaobjects}
    for x in arenaobjects:
        for key,y in x.fields.items():
            if key != 'bit' and y is not None:
                assert id(y) in adomain
    assert connection.port.root is source.root
    ids = [id(x) for x in objects+arenaobjects]
    assert max(ids).bit_length() <= 8*struct.calcsize('P')
    return dict(source_descriptors=len(descriptors),source_and_target_Blocks=len(blocks),
                checked_Block_and_descriptor_references=refs,
                checked_shared_wire_successors=wiring_checks,
                checked_installer_control_vertices=installer_control_checks,
                original_source_and_connection_objects=len(arenaobjects),
                actual_pointer_bytes=struct.calcsize('P'),
                distinct_counted_object_ids=len(set(ids)),
                pointer_value_max_bits=max(ids).bit_length(),
                targets_are_objects=True,integer_to_reference_casts=0,
                shared_wiring_once=True,independent_mutable_initializer_headers=True,
                modeled_local_widths=[direct.width,b.width,source.a.width,connection.a.width],
                modeling_widths_are_not_host_pointer_widths=True)


def source_locations():
    wanted={
        'constructor_ir.py':('Assembler','primitives'),
        'evaluator_program.py':('define_evaluator','main'),
        'native_install.py':('BitKernel','Installer'),
        'service_native.py':('OriginalConnection',),
        'session.py':('prepare_source','OriginalPort'),
        'field_service.py':('Service',),
        'native-prerequisite/immutable_spine.py':('Direct',),
        'native-prerequisite/is2_service_certificate.py':('ServiceDirect','compose'),
    }
    out={}
    for filename,names in wanted.items():
        tree=ast.parse((BASE/filename).read_bytes())
        out[filename]={}
        for n in ast.walk(tree):
            if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names:
                out[filename][n.name]=[n.lineno,n.end_lineno]
        out[filename]['sha256']=sha(BASE/filename)
    return out


def main():
    start=time.monotonic()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('bounded demand bundle')))
    signal.alarm(CAPS['wall_seconds'])
    result=dict(candidate=im.IDENTITY,caps=CAPS,
                full_evaluator_executions=0,full_service_chronologies=0,
                evaluator_execution_reserve=0,actual_resource_offer=None)
    kernel=None
    try:
        a,primitives,raw,fixed,field,linked,depths=fixed_source()
        result['static']=dict(raw_states=len(raw['code']),reachable_states=len(fixed['code']),
            calls={f:dict(v) for f,v in a.edges.items()},call_depths=depths,
            primitive_names=primitives,field_nodes=len(linked.nodes),
            native_nodes=len(linked.nodes)+sum(old.NATIVE_COST[n.op] for n in linked.nodes))
        kernel=BitKernel()
        im.Installer=BoundedInstaller
        ne.Installer=BoundedInstaller
        first_bits=(BASE/'constructor/positive-full.input.bits').read_text()
        _,_,first_cache,_=parse_arithmetic(first_bits)
        pattern=[dict(row).get(h,0) for row in first_cache for h in range(15)]
        source=prepare_source(pattern)
        port=OriginalPort(source)
        connection=OriginalConnection(port)
        # Complete fixed code relocation, NO evaluator seed/loader/instruction.
        direct=im.Direct(linked,'',0,kernel,connection=connection,service_return=True)
        assert direct.steps==0 and not direct.output and not direct.allocated
        direct_seed=direct.inst.seed_census()
        result['rows']=[row(name,a,primitives,raw,linked,direct_seed) for name in CASES]
        # Separately mutable B; an actual paid initialization to ready at p=0.
        bg=b_graph()
        b,code=install_b(('R','M',('J',2)),'00',CAPS['B_reserve'],bg,
                        lambda *x,**kw:ne.NativeEngine(*x,kernel=kernel,**kw))
        b.seal(checkpoint)
        b.run(CAPS['B_native_steps'])
        assert b.waiting and b.last is bg.entries['ready']
        assert b.native_steps<=CAPS['B_native_steps']
        assert not port.calls and not connection.calls
        assert source_census(60)==dict(O=len(source.a.objects),A=source.a.slots(),C=source.a.bit_fields)
        result['installed']=dict(
            fixed_evaluator=dict(reserve=0,executed_steps=direct.steps,
                source_blocks=len(direct.inst.sources),native_census=direct.census(),
                seed_census=direct_seed,installation_steps=direct.inst.steps,
                installation_charges=dict(direct.inst.charges)),
            separate_B=dict(program=['R','M',['J',2]],input='00',requested_prefix=0,
                completed_B_instructions=0,initial_loader_complete=True,
                reserve=CAPS['B_reserve'],native_steps=b.native_steps,
                source_blocks=len(b.installation.sources),native_census=b.census(),
                seed_census=b.installation.seed_census(),installation_steps=b.installation.steps),
            original_source=dict(k=60,census=source.a.census(),actual_calls=0,
                factory='unchanged prepare_source; original positive-full public-cache fixture convention',
                source_pattern_sha256=hashlib.sha256(bytes(pattern)).hexdigest(),
                preparation_charges=dict(source.a.charges)),
            connection=dict(census=connection.a.census(),actual_calls=0),
            references=actual_reference_union(direct,b,source,connection))
        # Each existing B-pool cell has the same source/descriptor/target
        # increment. Program/input are installed WITHIN that pool.
        bc=b.census()
        b_intercept=dict(O=bc['objects']-3*CAPS['B_reserve'],
                         A=bc['refs']-13*CAPS['B_reserve'],
                         C=bc['bits']-16*CAPS['B_reserve'])
        assert b_intercept==dict(O=1281203,A=5223157,C=18722307)
        result['separate_B_affine_census']=dict(
            intercept=b_intercept,coefficient=dict(O=3,A=13,C=16),
            source_blocks_intercept=len(b.installation.sources)-CAPS['B_reserve'],
            semantic_parameters='finite program a, input ell, requested instruction prefix p',
            reserve='a+ell+6+4*C+D_B; C=ell+p+1; D_B=1+C*(C+1)/2+C*(C+1)',
            scope='ordinary source-layout derivation; one actual p=0 installation, not all-prefix tests')
        wires=1+sum(len(x) for x in kernel.wires.values())
        shared=dict(O=wires,A=wires,C=22*wires)
        result['shared_immutable_seed']=shared
        for row_data in result['rows']:
            components=(row_data['native_census_upper'],row_data['original_source_exact'],
                        row_data['connection_exact'],
                        dict(O=bc['objects'],A=bc['refs'],C=bc['bits']))
            common={k:sum(c[k] for c in components)-shared[k] for k in ('O','A','C')}
            row_data['common_native_with_retained_actual_B_upper']=common
            row_data['common_native_bit_length']={k:v.bit_length() for k,v in common.items()}
            row_data['common_native_at_64bit_carriers_upper']=common['O']+64*common['A']+common['C']
            row_data['chosen_upper_object_count_fits_64bit_names']=common['O']<2**64
            row_data['byte_work_offer_comparison']=dict(host_upper_complete=False,
                H_reserved=None,W_backed=None,demand_fits_offer=None,
                source_native_upper_complete=True,whole_host_demand_complete=False,
                selected_reserve_layout_fits_runtime=False,
                upper_failure_proves_impossibility=False)
        # The module-level label maps are retained preparation too.
        roots=[raw,fixed,field,linked,direct,b,source,port,connection,kernel,im.GRAPH_LABELS]
        result['installed']['owned_reachable_payload']=payload_census(roots)
        model_union_objects=direct.census()['objects']+b.census()['objects']+len(source.a.objects)+len(connection.a.objects)
        # This deliberately sums the two seed censuses, thus overcounts shared
        # wiring; it is a safe containment check, not the common-cut census.
        assert model_union_objects<=CAPS['modeled_objects']
        result['runtime']=dict(python=sys.version,implementation=sys.implementation.name,
            executable=sys.executable,pointer_bytes=struct.calcsize('P'),
            Py_ssize_t_max=sys.maxsize,int_bits_per_digit=sys.int_info.bits_per_digit,
            int_digit_bytes=sys.int_info.sizeof_digit,shallow_sizes=host_shapes(),
            Wire_bytes=C.sizeof(Wire),Init_bytes=C.sizeof(Init),
            physical_memory_bytes=os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES'),
            limits={n:list(resource.getrlimit(getattr(resource,n))) for n in
                    ('RLIMIT_AS','RLIMIT_DATA','RLIMIT_RSS','RLIMIT_CPU','RLIMIT_NOFILE')},
            facts_are_protected_offers=False,
            include_path=sysconfig.get_path('include'))
        result['source_locations']=source_locations()
        result['remaining']=dict(
            author_owned=['complete peak CPython allocation/temporary/allocator and runtime work bounds',
                'runtime-capacity-legal sufficient reserve or a proved sufficient adaptive bound for the full three sessions',
                'complete evaluator session execution observers and their retained/transient copies',
                'host preparation/observer work substitution for each requested B program/input/prefix'],
            external_facts=['no same-scope protected host-byte reserve or backed host-work offer was exposed by inspected facilities',
                'retained shared external owner material and its actual overlap/capacity has no supplied manifest'],
            not_external=['local synchronous OriginalPort/OriginalConnection provider completion',
                'local Direct send/End/Ack/indirect Return completion under paid allocation/rights'],
            downstream=['full original evaluator/session runs','strict exact_checks --require-all',
                'author D reconciliation','fresh native and two CLI acceptance','CI/merge/first dev push'])
        result['status']='source-cost-recurrence-and-full-code-layout-complete; whole-host-demand-and-offer-incomplete'
        result['bundle_exit']=0
    finally:
        if kernel is not None:kernel.close()
        signal.alarm(0)
        result['owned_jobs_joined']=True
        result['elapsed_seconds']=time.monotonic()-start
        result['peak_RSS_bytes']=rss()
        assert result['peak_RSS_bytes']<=CAPS['monitored_RSS_bytes']
        result['program_sha256']=sha(Path(__file__))
        dest=HERE/'is2-complete-demand-results.json'
        tmp=dest.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        tmp.replace(dest)
    print(json.dumps({k:result[k] for k in ('status','bundle_exit','elapsed_seconds','peak_RSS_bytes','owned_jobs_joined')}))


if __name__=='__main__':
    main()
