# SPDX-License-Identifier: Apache-2.0
"""Disposable fixed-table audit; never calls the production synthesis entry."""
import collections
import gzip
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import time
import z3

ROOT = Path(__file__).resolve().parent / 'readiness'
spec = importlib.util.spec_from_file_location('audited_readiness', ROOT / 'readiness.py')
production = importlib.util.module_from_spec(spec)
spec.loader.exec_module(production)
z3.set_param('memory_max_size', 512)
start = time.monotonic()
counts = collections.Counter()

def states(mask, n):
    return frozenset(q for q in range(n) if mask >> q & 1)

def mask(subset):
    return sum(2 ** q for q in subset)

def step(table, subset, a):
    return frozenset(table[a][q] for q in subset)

def act(table, subset, word):
    for a in word:
        subset = step(table, subset, a)
    return subset

def literal(table, S, A, bound):
    for length in range(bound + 1):
        for word in itertools.product(range(2), repeat=length):
            target = act(table, S, word)
            if not A <= target or len(target) != len(S):
                return False, word, target
    return True, None, None

def sc(table):
    n = len(table[0])
    for source in range(n):
        reached = {source}
        for _ in range(n - 1):
            reached |= {table[a][q] for q in reached for a in range(2)}
        if len(reached) != n:
            return False
    return True

def image_layers(table, bound):
    layer = {frozenset(range(len(table[0])))}
    layers = [layer]
    for _ in range(bound):
        layer = {step(table, sub, a) for sub in layer for a in range(2)}
        layers.append(layer)
    return layers

def certificate_check(table, S, A, bound, certificate):
    nodes = dict(certificate['nodes'])
    assert nodes[mask(S)] == 0
    assert len(nodes) == certificate['visited'] == certificate['processed']
    expected_edges = set()
    for encoded, depth in nodes.items():
        sub = states(encoded, len(table[0]))
        assert A <= sub and len(sub) == len(S) and 0 <= depth <= bound
        if depth < bound:
            for a in range(2):
                dest = mask(step(table, sub, a))
                assert dest in nodes and nodes[dest] <= depth + 1
                expected_edges.add((encoded, a, dest))
    assert set(map(tuple, certificate['edges'])) == expected_edges
    distance = {mask(S): 0}
    for depth in range(bound):
        for encoded in list(distance):
            if distance[encoded] != depth:
                continue
            for a in range(2):
                dest = mask(step(table, states(encoded, len(table[0])), a))
                distance.setdefault(dest, depth + 1)
    assert nodes == distance

# Exhaust every binary three-state table, every nonempty proper A in every S,
# and every bound 0..3. Literal enumeration uses sets, not production bitsets.
for entries in itertools.product(range(3), repeat=6):
    table = [list(entries[:3]), list(entries[3:])]
    assert production.strong(table) == sc(table)
    counts['strong_connectivity_tables'] += 1
    for S_mask in range(1, 8):
        S = states(S_mask, 3)
        for A_mask in range(1, 8):
            A = states(A_mask, 3)
            if not A < S:
                continue
            for bound in range(4):
                expected, _, _ = literal(table, S, A, bound)
                got = production.oracle(table, S_mask, A_mask, bound, True)
                assert got['survives'] == expected
                counts['oracle_literal_cases'] += 1
                if expected:
                    certificate_check(table, S, A, bound, got)
                    counts['surviving_certificates_checked'] += 1
                else:
                    word = got['word']
                    assert len(word) <= bound
                    target = act(table, S, word)
                    assert mask(target) == got['mask']
                    assert not A <= target or len(target) < len(S)
                    counts['rejection_words_replayed'] += 1

# Independent one-letter controls distinguish the two failure disjuncts.
controls = [([[1,2,0],[0,1,2]], 'avoidance'),
            ([[0,0,2],[0,1,2]], 'compression')]
for table, kind in controls:
    result = production.oracle(table, 3, 1, 1)
    assert result['word'] == [0] and not result['survives']
    target = act(table, frozenset({0,1}), result['word'])
    assert ((0 not in target and len(target) == 2) if kind == 'avoidance'
            else (0 in target and len(target) == 1))
    counts[kind + '_only_control'] += 1

# All table variables are fixed before every solver check: no synthesis.
tables = list(itertools.product(range(3), repeat=6))
for entries in tables[::37]:
    table = [list(entries[:3]), list(entries[3:])]
    for S_mask, reach_bound, reset_bound in ((3,0,0),(3,1,2),(5,2,3),(7,0,3)):
        S, A = states(S_mask,3), frozenset({0})
        enc = production.Encoding(3,S_mask,1,2,reach_bound,reset_bound,2000)
        enc.add_cut((0,)); enc.add_cut((1,0))
        for a in range(2):
            for q in range(3):
                enc.solver.add(enc.table[a][q] == table[a][q])
        layers = image_layers(table, max(reach_bound,reset_bound))
        expected = (sc(table) and any(S in layer for layer in layers[:reach_bound+1])
                    and any(len(sub) == 1 for layer in layers[:reset_bound+1] for sub in layer)
                    and all(A <= act(table,S,w) and len(act(table,S,w)) == len(S)
                            for w in ((),(0,),(1,0))))
        status = enc.solver.check()
        assert status != z3.unknown
        assert (status == z3.sat) == expected
        counts['fixed_table_encoding_cases'] += 1
        if status == z3.sat:
            decoded = enc.decode(enc.solver.model())
            assert decoded['table'] == table
            assert len(decoded['reaching_word']) <= reach_bound
            assert len(decoded['reset_word']) <= reset_bound
            assert act(table,frozenset(range(3)),decoded['reaching_word']) == S
            assert len(act(table,frozenset(range(3)),decoded['reset_word'])) == 1
            counts['decoded_fixed_witnesses_replayed'] += 1

# Explicit nonvacuous admissible control, replayed from fixed original letters.
table = [[1,2,0],[0,1,0]]
control_layers = image_layers(table,4)
assert not any(len(sub) == 1 for layer in control_layers[:4] for sub in layer)
assert any(len(sub) == 1 for sub in control_layers[4])
enc = production.Encoding(3,3,1,0,1,4,2000)
for a in range(2):
    for q in range(3): enc.solver.add(enc.table[a][q] == table[a][q])
assert enc.solver.check() == z3.sat
decoded = enc.decode(enc.solver.model())
assert act(table,frozenset(range(3)),decoded['reaching_word']) == frozenset({0,1})
assert len(act(table,frozenset(range(3)),decoded['reset_word'])) == 1
counts['nonvacuous_fixed_witness_control'] += 1

# Construct/serialize the n19 encoding, never solve it. Verify retained SMT.
raw = gzip.decompress((ROOT / 'initial-admissibility.smt2.gz').read_bytes())
assert hashlib.sha256(raw).hexdigest() == '6efe92bc0fa40f8ec9314c4e04d2b972dfe9e44289113886d50e18b9e792ad30'
enc = production.Encoding(19,8191,1,17,126,220,20000)
generated = enc.solver.sexpr().encode()
assert generated == raw
assert len(z3.parse_smt2_string(raw.decode())) == 427
before = enc.counts()
enc.add_cut((0,))
assert len(enc.solver.assertions()) == 429
# Check the n19 bit width against set images without solving any n19 formula.
for table in ([[ (q+1)%19 for q in range(19)], [18-q for q in range(19)]],
              [[0 for _ in range(19)], list(range(19))],
              [[(q*7+3)%19 for q in range(19)], [(q*q+1)%19 for q in range(19)]]):
    substitutions = [(enc.table[a][q],z3.BitVecVal(table[a][q],5))
                     for a in range(2) for q in range(19)]
    for encoded in [0,8191,(1<<19)-1] + [1<<q for q in range(19)]:
        for a in range(2):
            symbolic = enc.img(z3.BitVecVal(encoded,19),a)
            got = z3.simplify(z3.substitute(symbolic,*substitutions)).as_long()
            assert got == mask(step(table,states(encoded,19),a))
            counts['n19_symbolic_images_no_solver'] += 1
result = {'counts':dict(counts),'n19_initial_counts':before,
          'n19_smt_identical':True,'n19_solver_checks':0,
          'source_sha256':hashlib.sha256((ROOT/'readiness.py').read_bytes()).hexdigest(),
          'smt_sha256':hashlib.sha256(raw).hexdigest(),
          'audit_control_correction':'Initial reset bound 3 was too short; independent enumeration gives minimum 4; fixed audit control to bound 4, no production change',
          'z3_version':z3.get_version_string(), 'elapsed_seconds':time.monotonic()-start,
          'verdict':'all fixed-table audit assertions passed'}
Path('audit-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
