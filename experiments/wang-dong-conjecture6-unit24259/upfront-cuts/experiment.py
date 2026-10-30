# SPDX-License-Identifier: Apache-2.0
"""One preregistered upfront-cut feasibility check; never a solver campaign."""
import collections, datetime, gzip, hashlib, json, resource, time
from pathlib import Path
import z3
from encoding import Encoding, oracle
ROOT = Path(__file__).resolve().parent

def save(name, data):
    target = ROOT / name
    tmp = target.with_suffix(target.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(target)

def action(table, states, word):
    states = set(states)
    for letter in word:
        assert letter in (0, 1)
        states = {table[letter][q] for q in states}
    return states

def replay(model):
    table = model['table']; Q = set(range(19)); S = set(range(13)); A = {0}
    assert len(table) == 2 and all(len(row) == 19 for row in table)
    assert all(type(v) is int and v in Q for row in table for v in row)
    # Independent set-based transitive closures from every state.
    for root in Q:
        seen = {root}
        while True:
            enlarged = seen | {table[a][q] for a in (0, 1) for q in seen}
            if enlarged == seen: break
            seen = enlarged
        assert seen == Q
    rw, zw = model['reaching_word'], model['reset_word']
    assert len(rw) <= 126 and action(table, Q, rw) == S
    assert len(zw) <= 220 and len(action(table, Q, zw)) == 1
    for w in ([], [0], [1]):
        image = action(table, S, w)
        assert A <= image and len(image) == 13
    return {'total_entries':38,'strong_connectivity_roots':19,'reaching_length':len(rw),
            'reset_length':len(zw),'initial_cuts_replayed':3,'passed':True}

def verify_certificate(table, cert):
    # Separate set representation; independently establish full closure and depths.
    def unpack(mask): return frozenset(q for q in range(19) if mask & (1 << q))
    depth = {unpack(mask): d for mask, d in cert['nodes']}
    assert len(depth) == len(cert['nodes']) and depth[frozenset(range(13))] == 0
    assert all(0 <= d <= 17 and len(s) == 13 and 0 in s for s,d in depth.items())
    edges = {(unpack(src), a): unpack(dst) for src,a,dst in cert['edges']}
    assert len(edges) == len(cert['edges'])
    assert set(edges) == {(s,a) for s,d in depth.items() if d < 17 for a in (0,1)}
    for (s,a),t in edges.items():
        assert frozenset(action(table,s,[a])) == t and t in depth
    todo=collections.deque([frozenset(range(13))]); found={todo[0]:0}
    while todo:
        s=todo.popleft()
        if found[s] == 17: continue
        for a in (0,1):
            t=edges[s,a]
            if t not in found: found[t]=found[s]+1; todo.append(t)
    assert found == depth
    return {'nodes':len(depth),'edges':len(edges),'minimum_depths_and_complete_closure':True}

def main():
    start=time.monotonic(); z3.set_param('memory_max_size',1536)
    result={'status':'building','solver_version':z3.get_version_string(),
            'target':{'n':19,'m':6,'s':13,'r':1,'B':17,'S':list(range(13)),'A':[0],
                      'alphabet':[0,1],'reach_bound':126,'reset_bound':220,'original_letter_cost':1},
            'timeout_ms':60000,'z3_memory_limit_mb':1536,'solver_checks':0,'candidate':None,
            'started_utc':datetime.datetime.now(datetime.UTC).isoformat()}
    save('outcome.json',result)
    t=time.monotonic(); enc=Encoding(19,8191,1,17,126,220,60000)
    enc.add_cut([0]); enc.add_cut([1]); result['build_seconds']=time.monotonic()-t
    result['formula_counts']=enc.counts(); result['initial_cut_words']=enc.cuts
    t=time.monotonic(); smt=enc.solver.sexpr().encode()
    (ROOT/'upfront.smt2.gz').write_bytes(gzip.compress(smt,mtime=0))
    result['smtlib_bytes']=len(smt); result['smtlib_sha256']=hashlib.sha256(smt).hexdigest()
    del smt; result['serialization_seconds']=time.monotonic()-t
    result['status']='checking'; result['solver_checks']=1; save('outcome.json',result)
    t=time.monotonic(); status=enc.solver.check(); result['solve_seconds']=time.monotonic()-t
    result['status']=str(status); result['solver_statistics']={str(k):v for k,v in enc.solver.statistics()}
    if status == z3.unknown: result['reason_unknown']=enc.solver.reason_unknown()
    save('outcome.json',result)
    if status == z3.sat:
        # Persist raw decoded entries and complete symbolic words before any replay/oracle.
        model=enc.solver.model()
        decoded={'table':[[model.eval(t,model_completion=True).as_long() for t in row] for row in enc.table],
                 'reaching_word':[int(z3.is_true(model.eval(a,model_completion=True))) for a in enc.reach_letters],
                 'reset_word':[int(z3.is_true(model.eval(a,model_completion=True))) for a in enc.reset_letters]}
        save('decoded-raw.json',decoded)
        # Select actual reaching prefix independently using set images.
        Q=set(range(19)); states=Q; hit=0 if states==set(range(13)) else None
        for i,a in enumerate(decoded['reaching_word'],1):
            states=action(decoded['table'],states,[a])
            if hit is None and states==set(range(13)): hit=i
        assert hit is not None
        decoded['reaching_word']=decoded['reaching_word'][:hit]
        save('model.json',decoded)  # Before admissibility replay and full oracle.
        t=time.monotonic(); result['admissibility_replay']=replay(decoded)
        result['admissibility_replay_seconds']=time.monotonic()-t; save('outcome.json',result)
        t=time.monotonic(); checked=oracle(decoded['table'],8191,1,17,True)
        result['oracle_seconds']=time.monotonic()-t
        save('oracle.json',checked)
        result['oracle']={k:v for k,v in checked.items() if k not in ('nodes','edges')}
        if checked['survives']:
            t=time.monotonic(); result['certificate_replay']=verify_certificate(decoded['table'],checked)
            result['certificate_replay_seconds']=time.monotonic()-t
            result['candidate']='model.json with oracle.json; pending distinct-seat and kernel verification'
        else:
            word=checked['word']; image=action(decoded['table'],range(13),word)
            assert len(word) <= 17 and (0 not in image or len(image) < 13)
            assert sum(1 << q for q in image)==checked['mask']
            # Independent set BFS verifies shortest failure, without bitset helpers.
            todo=collections.deque([(frozenset(range(13)),0)]); seen={todo[0][0]}; minimum=None
            while todo:
                states,d=todo.popleft()
                if 0 not in states or len(states)!=13: minimum=d; break
                if d==len(word): continue
                for a in (0,1):
                    nxt=frozenset(action(decoded['table'],states,[a]))
                    if nxt not in seen: seen.add(nxt); todo.append((nxt,d+1))
            assert minimum==len(word)
            result['violation_replay']={'word':word,'image':sorted(image),'cardinality':len(image),
                'contains_A':0 in image,'shortest_length':minimum,'passed':True}
            assert len(word)>=2, 'No advance: upfront one-letter cut violated'
            result['engineering_advance']='admissible model beyond both one-letter cuts; no counterexample'
    result['elapsed_seconds']=time.monotonic()-start
    result['peak_process_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['finished_utc']=datetime.datetime.now(datetime.UTC).isoformat()
    save('outcome.json',result)

if __name__=='__main__':
    try: main()
    except BaseException as exc:
        old=json.loads((ROOT/'outcome.json').read_text())
        old['status']='mismatch_or_error'; old['error_type']=type(exc).__name__; old['error']=str(exc)
        old['peak_process_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        save('outcome.json',old)
        raise
