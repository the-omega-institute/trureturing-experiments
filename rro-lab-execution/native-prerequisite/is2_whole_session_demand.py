"""Source-only certificate for the fixed IS2 consumer, with zero execution ticks.

Assemble public code, count syntactic sites, prove/check public frame guards,
and substitute a two-coordinate recurrence. No native installation, evaluator,
B machine, provider or chronology is invoked. One bounded invocation can wait
for the author to insert its printed numerical table before final source checks.
"""
from __future__ import annotations

import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import select
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
THEORY = Path('/Users/auricstudio/trureturing-cad0ceaa-op')
VOLUME = 'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_IMMUTABLE_SPINE_NATIVE_REALIZATION.md'
PINS = {
    'positive-full': (4, 387, 74, 2176, '37dfee71eda5a8fbee3d373aab446982f15976785ef965ec2e4a1b6ea98c5402'),
    'positive-partial': (4, 386, 73, 2115, '1a4ae049cf81b1f434398aa05cb8eabcf4d8689f7873f8d75f9ae1e870c34764'),
    'null-nine-distinct-left': (9, 839, 149, 9676, 'f5518d034eef5bd00c87ff9d43c5ca878ef9d012f47e04bedcbe45163ee3e29d')}
CAPS = {'seconds': 600, 'monitored_RSS_bytes': 2*1024**3,
        'evaluator_ticks': 0, 'B_ticks': 0, 'installations': 0,
        'provider_calls': 0, 'chronology_runs': 0}
sys.path[:0] = [str(HERE), str(BASE)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checkpoint():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        rss *= 1024
    assert rss <= CAPS['monitored_RSS_bytes']
    return rss


def forbidden(*args, **kwargs):
    raise AssertionError('execution/installation forbidden in source-only bundle')


def source_counts():
    import constructor_ir as ci
    import evaluator_program as ep
    import immutable_spine as im
    import is2_service_certificate as sc
    import native_install as ni
    import native_engine as ne
    # Fail if a static helper accidentally crosses an execution boundary.
    ci.run_graph = forbidden
    im.Direct.__init__ = im.Direct.tick = im.Direct.run = forbidden
    ni.BitKernel.__init__ = forbidden
    ni.Installer.__init__ = ni.Installer.tick = ni.Installer.run = forbidden
    ne.NativeEngine.__init__ = ne.NativeEngine.microtick = forbidden

    class CountAssembler(ci.Assembler):
        def __init__(self):
            super().__init__()
            self.calls_by_body = {}
            self.own = {}

        def call(self, name, outs, *args):
            self.calls_by_body.setdefault(self.context, Counter())[name] += 1
            return super().call(name, outs, *args)

        def counted(self, name, body):
            def wrapped(a):
                start = len(self.g.code)
                body(a)
                emitted = self.g.code[start:]
                self.own[name] = {'instructions': len(emitted)+(1 if name == 'main' else 12),
                                  'pushes': sum(x[0] == 'push' for x in emitted)}
                checkpoint()
            return wrapped

        def define(self, name, args, outs, body):
            return super().define(name, args, outs, self.counted(name, body))

    with im.shared_assembler():
        a = CountAssembler()
        ci.primitives(a)
        ep.define_evaluator(a)
        raw = a.compile(a.counted('main', ep.main))
    fixed = im.prune(raw)
    assert hashlib.sha256(im.canonical(fixed)).hexdigest() == '7257d4ad4b5730860ad811d41c3057ea6accdf1c413dbf7af75be206a5b24595'
    field = im.field_graph(fixed)
    linked = sc.compose(field)
    assert linked.fingerprint() == 'db4fb22177743a3aafd286c5431b42d09460ab47cfe4ca488d40dd1f28600b4a'
    assert (len(fixed['code']), len(fixed['registers']), len(linked.nodes),
            len(linked.refs), len(linked.bits)) == (10983, 360, 62069, 604, 57)
    # Pool advancement is by get, not a reuse/reset of extracted cells.
    assert not any(x.op == 'copy' and x.args[0] == 'pool' for x in linked.nodes)
    checkpoint()
    return a.own, {f: dict(c) for f, c in a.calls_by_body.items()}, len(sc.service_graph().nodes)


def guard_rows(name, cache, cylinders):
    """Only represented masks/values and conflict flags; never evaluate a law."""
    n = len(cache)
    restrictions = [cache]
    for cyl in cylinders:
        if cyl is None:
            restrictions.append(None)  # conflict/empty-event flag, valid=0
            continue
        joined = []
        for left, right in zip(cache, cyl):
            x, y = dict(left), dict(right)
            assert all(h not in x or x[h] == b for h, b in y.items())
            joined.append(tuple(sorted((x | y).items())))
        restrictions.append(tuple(joined))
    facts = []
    for data in restrictions:
        if data is None:
            facts.append({'valid': False, 'normalizer_returns_before_side': True})
            continue
        rows = [dict(x) for x in data]
        left = [tuple(x[h] for h in range(7)) for x in rows]
        assert all(all(h in x for h in range(7)) for x in rows)
        if name == 'null-nine-distinct-left':
            assert len(set(left)) == 9
            assert all(all(x[h] == 0 for h in range(7, 15)) for x in rows)
            facts.append({'valid': True, 'r': 9, 's': 1,
                          'distinct_edges_at_one_right': 9, 'map_settings': 0})
        else:
            assert len(set(left)) == 1
            assert all(all(x[h] == 0 for h in range(14)) for x in rows)
            full = all(14 in x for x in rows)
            if full:
                assert all(x[14] == 0 for x in rows)
            else:
                assert n == 4 and all(14 in x for x in rows[:3]) and 14 not in rows[3]
            facts.append({'valid': True, 'r': 1, 'right_full': full,
                          'right_partitions_visited_upper': 1 if full else 15,
                          'right_blocks_if_count_positive_upper': 1 if full else 2,
                          'countlabels_allowed_labels_upper': 2})
    return facts


def recurrence(name, n, length, t, h, own, calls):
    m0 = 64*n+4*h+64
    w = 4*m0+8
    e = 8*n+2
    partial = name == 'positive-partial'
    null = name == 'null-nine-distinct-left'
    subset = 16 if partial else 1
    records = max(15*n, 36, n*n, subset, 2*n)
    maps = 0 if null else 36*(27**(2 if partial else 1))
    p_right = 15 if partial else 1
    # Reuse the established all-domain output bound; it is not a run count.
    m = 64*(n+length+t+2)
    y = 64*(t+2)*(m+1)
    lam = {f: 1 for f in own}
    for f in ('add','sub','mul','div','encu','encv','readu'):
        lam[f] = w+2
    lam.update(gcd=2*w+3, popword=e+2, prepend=e+2,
               get=records+2, put=records+2, zeros=records+2,
               concat=y+2, readv=length+2,
               main=t+2, readcache=15*n+2, conjoin=15*n+2,
               blocks=n+2, falling=n+2, nextpart=(n+2)**2,
               side=16*(n+2)**2,
               countlabels=max(257*9*(n+1), 2*(subset+1)*(n+1)),
               odometer=2*n+2, inject=(maps+1)*(2*n+2),
               normalizer=max(n+2,p_right+1))
    edges = {f: dict(g) for f,g in calls.items()}
    edges['normalizer']['countlabels'] = 1 if partial else 0
    edges['normalizer']['nextpart'] = 1 if partial else 0
    if null:
        for f in ('incident','eq','odometer'):
            edges['inject'][f] = 0
    memo, active = {}, set()
    def cost(f):
        if f in memo:
            return memo[f]
        assert f not in active
        active.add(f)
        children = [(multiplicity, cost(g)) for g,multiplicity in edges.get(f,{}).items() if multiplicity]
        z = {unit: lam[f]*(own[f][unit]+sum(m*c[unit] for m,c in children))
             for unit in ('instructions','pushes')}
        active.remove(f)
        assert z['pushes'] <= z['instructions']
        memo[f] = z
        return z
    result = cost('main')
    return dict(m0=m0,w=w,entry_bits=e,records=records,subset_entries=subset,
                maps_per_inject_upper=maps,right_partitions_upper=p_right,
                multipliers=lam,calls_after_guard=edges,body_costs=memo,
                evaluator_instructions_upper=result['instructions'],
                evaluator_pushes_upper=result['pushes'],Y=y)


def row(name, own, calls, f0):
    from framing import Reader, parse_arithmetic, v
    from is2_chronology import acquisition
    raw = (BASE/'constructor'/(name+'.input.bits')).read_text()
    n, theta, cache, cylinders = parse_arithmetic(raw)
    words = acquisition(n, cache)
    assert (n,len(raw),len(words),sum(map(len,words)),hashlib.sha256(raw.encode()).hexdigest()) == PINS[name]
    facts = guard_rows(name,cache,cylinders)
    rd = Reader(raw)
    rd.u()
    assert rd.take(1) == '0'
    prior_a, prior_b = rd.u(), rd.u()
    assert (prior_a, prior_b) == (1,3)
    hn, ha, hb, ht = [max(1,x.bit_length()) for x in (n,prior_a,prior_b,len(cylinders))]
    q, s = len(words)+1, sum(map(len,words))+2
    commands = [raw[:rd.pos]]+['01'+v(w.translate(str.maketrans('LR','01'))) for w in words]
    commands += ['00'+v(raw), '01'+v('00'), '11']
    i = sum(map(len,commands))
    assert len(commands) == q+3
    # Compute cylinder lengths directly from masks, without normalizer results.
    lengths = [1 if x is None else 1+15*n+sum(len(r) for r in x) for x in cylinders]
    nonempty = sum(x is not None for x in cylinders)
    initial = 16*n+2*(hn+ha+hb)+5
    validation = (2*len(raw)+31*n+len(cylinders)+2*(hn+ha+hb+ht)+9
                  +sum(2*c+1 for c in lengths)+nonempty*(15*n+1))
    af = i+(q+3)+2+initial+validation+5*s+14*q+4
    assert af == {'positive-full':18713,'positive-partial':18264,'null-nine-distinct-left':74987}[name]
    tf = 16*(f0+1)*(af+q+4)
    ev = recurrence(name,n,len(raw),len(cylinders),max(ha,hb),own,calls)
    k = 5+len(raw)+af+ev['evaluator_pushes_upper']
    nn = k+226380
    rr, cc = 4*k+903480, 8*k+4740143
    native = dict(O=3*k+777713,A=13*k+3037783,C=16*k+11648685)
    source_k = 15*n
    d = (5*source_k**2+13*source_k+6)//2
    src = dict(O=11*source_k+11+d,A=16+4*(8*source_k+8)+2*(3*source_k+3)+d,C=d)
    conn = dict(O=2*s+5*q+2,A=2*s+7*q+12,C=2*s+3*q)
    b = dict(O=1281587,A=5224821,C=18724355)
    shared = dict(O=98404,A=98404,C=22*98404)
    common = {key:native[key]+src[key]+conn[key]+b[key]-shared[key] for key in native}
    install = 128*nn+16*rr+4*cc+64
    execution = 128*(ev['evaluator_instructions_upper']+len(raw)+2*ev['Y']+360+64)+5*tf
    output = 7*q+2*ev['Y']+6
    signals = i+q+3
    return dict(case=name,N=n,L=len(raw),t=len(cylinders),sha256=PINS[name][4],
                acquisition_reads=len(words),acquisition_operand_bits=s-2,
                q=q,S=s,I=i,input_signals=signals,input_latch_writes=3*signals,
                represented_prior=[prior_a,prior_b],prior_widths=[ha,hb],cylinder_lengths=lengths,
                guard_observations=facts,evaluator=ev,
                frontend=dict(initialization_extractions=initial,validation_extractions=validation,
                              exact_extractions=af,F0=f0,field_actions_upper=tf),
                K=k,K_bit_length=k.bit_length(),source_blocks=nn,native_census=native,
                original_source=src,connection=conn,common_with_B=common,
                installation_transitions_upper=install,execution_transitions_upper=execution,
                output_bits_upper=output,receipt_bits=2*s+3*q,
                frozen_entries_upper=2*len(raw)+ev['evaluator_pushes_upper']+5,
                run_join_calls=signals+1,run_join_characters_upper=(signals+1)*output,
                chosen_reserve_payload_bytes=500*k,
                chosen_K_fits_Py_ssize_t=k<=sys.maxsize,
                chosen_K_pointer_array_bytes_fit_Py_ssize_t=8*k<=sys.maxsize,
                chosen_common_objects_fit_64bit_names=common['O']<2**64,
                host_upper_complete=False,H_reserved=None,W_backed=None,
                full_run_release_ready=False,full_session_executed=False)


def final_source_checks():
    old = subprocess.run(['git','show','1c80c55f641d0a536216eae20358cfb6aaddfe64:'+VOLUME],
                         cwd=THEORY,capture_output=True,check=True).stdout
    now = (THEORY/VOLUME).read_bytes()
    assert len(old.splitlines()) == 1744 and now.startswith(old)
    assert len(now.splitlines()) <= 5000 and b'\r' not in now
    assert now.endswith('## 追加锚（本行以下为增补区）\n'.encode())
    assert now[len(old):].lstrip().startswith(b'## 15.')
    expected_main = {VOLUME}
    expected_private = {'rro-lab-execution/native-prerequisite/is2_chronology.py',
                        'rro-lab-execution/native-prerequisite/is2_whole_session_demand.py',
                        'rro-lab-execution/native-prerequisite/is2-whole-session-demand-results.json'}
    checks = []
    for root, expected in ((THEORY,expected_main),(BASE.parent,expected_private)):
        p = subprocess.run(['git','diff','--check'],cwd=root,capture_output=True,text=True)
        assert p.returncode == 0, p.stdout+p.stderr
        status = subprocess.run(['git','status','--porcelain'],cwd=root,capture_output=True,text=True,check=True)
        found = {line[3:] for line in status.stdout.splitlines()}
        assert found <= expected, (root,found)
        checks.append({'command':'git diff --check','cwd':str(root),'exit':p.returncode,
                       'changed_paths':sorted(found)})
    for path in (HERE/'is2_chronology.py', Path(__file__)):
        ast.parse(path.read_bytes())
    return dict(prefix_bytes_preserved=len(old),prefix_sha256=hashlib.sha256(old).hexdigest(),
                theory_LF=len(now.splitlines()),checks=checks)


def main():
    assert sys.version_info >= (3,12) and __debug__
    start = time.monotonic()
    signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('source-only allowance')))
    signal.alarm(CAPS['seconds'])
    own,calls,f0 = source_counts()
    rows = [row(name,own,calls,f0) for name in PINS]
    result = dict(schema=1,identity='IS2-source-specific-whole-session-demand/v1',
                  caps=CAPS,ordinary_proof_location=VOLUME+' §§15–18',
                  own_sites=own,syntactic_calls=calls,rows=rows,
                  source_hashes={str(p.relative_to(BASE)):sha(p) for p in [
                      BASE/'constructor_ir.py',BASE/'evaluator_program.py',BASE/'field_graph.py',
                      BASE/'field_service.py',BASE/'field_io.py',BASE/'native_install.py',
                      BASE/'native_engine.py',BASE/'native_bits.c',BASE/'session.py',BASE/'lab.py',
                      BASE/'service_native.py',HERE/'immutable_spine.py',HERE/'is2_service_certificate.py',
                      HERE/'is2_chronology.py',Path(__file__)]},
                  evaluator_ticks=0,B_ticks=0,native_installations=0,chronology_runs=0,
                  H_reserved=None,W_backed=None,host_upper_complete=False,
                  full_run_release_ready=False,all_six_outcomes_complete=False)
    print(json.dumps({'own_sites':own,'F0':f0,'rows':[
        {key:r[key] for key in ('case','K','K_bit_length','execution_transitions_upper','chosen_reserve_payload_bytes')}
        | {'Sx':r['evaluator']['evaluator_instructions_upper'],
           'Px':r['evaluator']['evaluator_pushes_upper'],'AF':r['frontend']['exact_extractions'],
           'TF':r['frontend']['field_actions_upper']} for r in rows]},sort_keys=True),flush=True)
    if '--finish-on-stdin' in sys.argv:
        print('Awaiting verify on stdin; source-only bundle and 600-second alarm remain owned.',flush=True)
        while not select.select([sys.stdin],[],[],1)[0]:
            checkpoint()
        assert sys.stdin.readline().strip() == 'verify'
    result['source_checks'] = final_source_checks()
    result['elapsed_seconds'] = time.monotonic()-start
    result['monitored_peak_RSS_bytes'] = checkpoint()
    result['owned_jobs_joined'] = True
    result['bundle_exit'] = 0
    target = HERE/'is2-whole-session-demand-results.json'
    temporary = target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    temporary.replace(target)
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ('bundle_exit','elapsed_seconds','monitored_peak_RSS_bytes','owned_jobs_joined')}),flush=True)


if __name__ == '__main__':
    main()
