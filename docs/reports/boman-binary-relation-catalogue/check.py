#!/usr/bin/env python3
"""Exact, self-contained check of the proposed finite positive catalogue.

The 51 published rules and eight larger witness matrices are taken from Magnus
Boman's MIT-licensed verification package, revision
75278b521c476b55b39321c5809172dc8c0fc1ab. Full attribution/license: README.md.
The three extra rules are Proposition 6.8 of arXiv:2609.13317v1.
Every property is re-evaluated here; no source proof or solver is trusted.
This program is experimental evidence, not a Lean proof.
"""
from itertools import product

ALL_P = 32767
EMPTY, NO_COMP, FALSE, UNIVERSAL = (1 << i for i in range(15, 19))

# Conjunction of lhs bits implies every rhs bit. Published rules first.
RULES = [(32, 18),
 (128, 2),
 (1, 516),
 (18, 32),
 (66, 32),
 (2049, 1024),
 (72, 4096),
 (4104, 64),
 (257, 1024),
 (24, 64),
 (2112, 256),
 (272, 64),
 (4352, 64),
 (4112, 64),
 (6144, 320),
 (16385, 8192),
 (16393, 4096),
 (76, 1),
 (388, 8),
 (8576, 8),
 (324, 8192),
 (16452, 8192),
 (16516, 8192),
 (16652, 2048),
 (192, 65536),
 (65536, 242),
 (65540, 131072),
 (66048, 32768),
 (73728, 32768),
 (40, 32768),
 (4098, 32768),
 (640, 32768),
 (32768, 29690),
 (32769, 131072),
 (32772, 131072),
 (33792, 131072),
 (3, 131072),
 (33, 131072),
 (129, 131072),
 (1026, 131072),
 (1056, 131072),
 (1152, 131072),
 (131072, 524287),
 (8, 8192),
 (1024, 2049),
 (2048, 16384),
 (4096, 25088),
 (4097, 72),
 (2056, 256),
 (1032, 262144),
 (262144, 32589),
 (70, 131072),
 (578, 32768),
 (8258, 32768)]

# Exact adjacency matrices as integer rows, bit y of row x indicating R(x,y).
# Independent exhaustive representatives for n <= 4; larger rows from source.
WITNESSES = {182: (20, 1, 2, 4, 10),
 498: (0, 0, 3),
 597: (1, 2, 7),
 8374: (2, 4, 9, 2),
 8590: (2, 5, 2),
 8717: (3, 7, 6),
 8758: (148, 137, 226, 101, 14, 19, 19, 120),
 8789: (1, 3, 5, 15),
 8974: (6, 13, 11, 6),
 9044: (1, 1, 1, 15),
 16626: (0, 0, 1),
 18802: (0, 1, 3),
 18930: (0, 1),
 18998: (28, 49, 434, 326, 104, 393, 167, 27, 211),
 19206: (26, 61, 59, 37, 47, 31),
 19280: (0, 3),
 24718: (2, 1, 8, 4),
 24758: (2, 4, 1, 1),
 25102: (10, 9, 48, 3, 36, 20),
 25142: (164, 145, 50, 7, 41, 138, 43, 28),
 25350: (24, 25, 25, 17, 9),
 26806: (2, 4, 1),
 26964: (1, 1, 3),
 27022: (2, 1),
 27190: (76, 41, 98, 52, 7, 81, 26),
 27398: (6, 5, 3, 7),
 27406: (6, 5, 3),
 28181: (3, 6, 5),
 28421: (3, 7, 7),
 28501: (1, 3),
 29188: (2, 6, 6),
 29277: (1, 2),
 29524: (1, 1, 1),
 29690: (0, 0),
 31572: (1, 1),
 31738: (0,),
 32589: (3, 3),
 32605: (1,)}

# Inclusion-maximal consistent closed positive masks omitting each property.
MAXIMAL_OMISSION = {1: (26806, 27022, 27190, 27406, 31572, 31738),
 2: (32605,),
 3: (31738,),
 4: (18930, 26806, 27190, 27398, 28501, 31572),
 5: (27022, 27406, 32589),
 6: (27022, 27406, 32605),
 7: (8717, 26806, 27022, 27190, 27406, 28181, 28421, 29188),
 8: (18802, 27190, 27406, 32605),
 9: (16626, 24718, 25102, 26806, 27190, 28181, 29277),
 10: (18930, 26806, 26964, 27022),
 11: (26806, 27022, 27190, 27406, 29277, 31572, 31738),
 12: (8590, 8974, 24718, 24758, 25102, 25142, 25350, 29277, 29524, 29690),
 13: (8717, 18930, 26806, 27022, 27190, 27406, 28501),
 14: (182, 597, 18930, 18998, 19206, 19280),
 15: (498, 8374, 8590, 8717, 8758, 8789, 8974, 9044)}

def profile(rows):
 n=len(rows);D=range(n)
 def r(x,y):return (rows[x]>>y)&1
 p=[
 all(r(x,x) for x in D),
 all(not r(x,x) for x in D),
 all(any(r(x,y) for y in D) for x in D),
 all(not r(x,y) or r(y,x) for x in D for y in D),
 all(not(r(x,y) and r(y,x)) or x==y for x in D for y in D),
 all(not r(x,y) or not r(y,x) for x in D for y in D),
 all(not(r(x,y) and r(y,z)) or r(x,z) for x in D for y in D for z in D),
 all(not(r(x,y) and r(y,z)) or not r(x,z) for x in D for y in D for z in D),
 all(r(x,y) or r(y,z) or not r(x,z) for x in D for y in D for z in D),
 all(not r(x,y) or any(r(x,z) and r(z,y) for z in D) for x in D for y in D),
 all(r(x,y) or r(y,x) for x in D for y in D),
 all(r(x,y) or r(y,x) or x==y for x in D for y in D),
 all(not(r(x,y) and r(x,z)) or r(y,z) for x in D for y in D for z in D),
 all(not(r(x,y) and r(x,z)) or any(r(y,w) and r(z,w) for w in D) for x in D for y in D for z in D),
 all(not(r(x,y) and r(x,z)) or r(y,z) or r(z,y) or y==z for x in D for y in D for z in D)]
 return sum((1<<i) for i,b in enumerate(p) if b)


def bits(mask):
    return tuple(i + 1 for i in range(15) if mask & (1 << i))


def horn(seed, rules=RULES):
    while True:
        before = seed
        for lhs, rhs in rules:
            if seed & lhs == lhs:
                seed |= rhs
        if seed == before:
            return seed


def auxiliary(rows, p):
    n = len(rows)
    if not any(rows):
        p |= EMPTY
    if not any(rows[y] for x in range(n) for y in range(n)
               if rows[x] & (1 << y)):
        p |= NO_COMP
    if all(row == (1 << n) - 1 for row in rows):
        p |= UNIVERSAL
    return p


def graph_profile(rows):
    """Independent bit-set/composition implementation for matrix cross-checks."""
    n = len(rows)
    full = (1 << n) - 1
    cols = [sum(1 << x for x in range(n) if rows[x] & (1 << y))
            for y in range(n)]
    comp = [0] * n
    for x in range(n):
        for y in range(n):
            if rows[x] & (1 << y):
                comp[x] |= rows[y]
    complement_comp = [0] * n
    for x in range(n):
        for y in range(n):
            if not rows[x] & (1 << y):
                complement_comp[x] |= full ^ rows[y]
    values = [
        all(rows[x] & (1 << x) for x in range(n)),
        not any(rows[x] & (1 << x) for x in range(n)),
        all(rows),
        rows == tuple(cols),
        all((rows[x] & cols[x]) & ~(1 << x) == 0 for x in range(n)),
        all(rows[x] & cols[x] == 0 for x in range(n)),
        all(comp[x] & ~rows[x] == 0 for x in range(n)),
        all(comp[x] & rows[x] == 0 for x in range(n)),
        all(complement_comp[x] & rows[x] == 0 for x in range(n)),
        all(rows[x] & ~comp[x] == 0 for x in range(n)),
        all(rows[x] | cols[x] == full for x in range(n)),
        all(rows[x] | cols[x] | (1 << x) == full for x in range(n)),
        all(rows[x] & ~rows[y] == 0
            for x in range(n) for y in range(n) if rows[x] & (1 << y)),
        all(rows[y] & rows[z]
            for x in range(n) for y in range(n) for z in range(n)
            if rows[x] & (1 << y) and rows[x] & (1 << z)),
        all(rows[x] & ~(rows[y] | cols[y] | (1 << y)) == 0
            for x in range(n) for y in range(n) if rows[x] & (1 << y)),
    ]
    return sum(1 << i for i, value in enumerate(values) if value)


def main():
    assert len(RULES) == 54
    assert len(WITNESSES) == 38
    # Source Section 12, five finite MCS model profiles (M4 is infinite).
    source_models = {
        'M1': ((1,), 32605),
        'M2': ((2, 1), 27022),
        'M3': ((6, 5, 3), 27406),
        'M5': ((2, 4, 1), 26806),
        'M6': ((0,), 31738),
    }
    for label, (rows, expected) in source_models.items():
        assert profile(rows) == graph_profile(rows) == expected
        print('FIDELITY', label, 'profile', expected)

    # Exhaustively reproduce all three assertions of source Proposition 6.8
    # on every labelled relation with n <= 4; also check every rule on them.
    counts = []
    small_profiles = set()
    for n in range(1, 5):
        count = 0
        for rows in product(range(1 << n), repeat=n):
            p = profile(rows)
            assert p == graph_profile(rows), (n, rows, p)
            small_profiles.add(p)
            sem = auxiliary(rows, p)
            for lhs, rhs in RULES:
                assert sem & lhs != lhs or sem & rhs == rhs, (rows, lhs, rhs)
            assert p & 70 != 70  # P2+P3+P7 inconsistent
            assert (p & 578 == 578) == (not any(rows))
            assert (p & 8258 == 8258) == (not any(rows))
            count += 1
        counts.append(count)
        print('FIDELITY_PROP_6_8', 'n', n, 'relations', count,
              '54_rules_sound', 'two_semantic_evaluators_agree')
    assert counts == [2, 16, 512, 65536]
    assert len(small_profiles) == 188
    print("FIDELITY_SECTION_11 distinct_small_profiles=188 relations=66066")

    for p, rows in WITNESSES.items():
        assert 1 <= len(rows) <= 9
        assert all(0 <= row < (1 << len(rows)) for row in rows)
        assert profile(rows) == graph_profile(rows) == p
        sem = auxiliary(rows, p)
        assert not horn(p) & FALSE
        assert all(sem & l != l or sem & r == r for l, r in RULES)
        print('MATRIX', p, 'n', len(rows), 'rows', ','.join(map(str, rows)),
              'properties', ','.join(map(str, bits(p))))

    closures = {a: horn(a) for a in range(1 << 15)}
    consistent = {a: c for a, c in closures.items() if not c & FALSE}
    closed = sorted(set(consistent.values()))
    assert len(consistent) == 7936
    assert len(closed) == 330
    assert len(closures) - len(consistent) == 24832
    maximal_union = set()
    for i in range(15):
        omitted = sorted({c & ALL_P for c in closed if not c & (1 << i)})
        maximal = tuple(c for c in omitted if not any(
            c != d and c & d == c for d in omitted))
        assert maximal == MAXIMAL_OMISSION[i + 1]
        maximal_union.update(maximal)
        print('MAXIMAL_OMISSION', i + 1, ','.join(map(str, maximal)))
    assert maximal_union == set(WITNESSES)

    # For every closure record a real model, and a real refuter of each
    # omitted property. Each referenced row was evaluated literally above.
    for c in closed:
        mask = c & ALL_P
        model = next(p for p in WITNESSES if p & mask == mask)
        refs = {}
        for i in range(15):
            if not mask & (1 << i):
                refs[i + 1] = next(p for p in MAXIMAL_OMISSION[i + 1]
                                   if p & mask == mask)
                assert WITNESSES[refs[i + 1]]
                assert not profile(WITNESSES[refs[i + 1]]) & (1 << i)
        print('CLOSED_MASK', mask, 'model', model, 'refuters',
              ','.join(str(i) + ':' + str(p) for i, p in refs.items()))

    obligations = 0
    for a, c in closures.items():
        assert a & c == a
        assert horn(c) == c
        for i in range(15):
            if not c & (1 << i):
                p = next(p for p in MAXIMAL_OMISSION[i + 1]
                         if p & (c & ALL_P) == c & ALL_P)
                assert p & a == a and not p & (1 << i)
            obligations += 1
    assert obligations == 491520
    print('SUMMARY antecedents=32768 consistent=7936 inconsistent=24832',
          'consistent_closures=330 obligations=491520 matrices=38 bound=9')
    print('SCOPE universal rule soundness remains a proof obligation;',
          'the exact finite calculation is not kernel verification')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
