"""Exact incidence-only recovery and literal/indexed inverse-label checks.

Python >=3.12; standard library; run with -B. All inputs are pinned Git objects.
The geometry constructor receives transported literal edges, never source w/v/z.
Coordinate access occurs only in the separate validation functions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb
import json
from pathlib import Path
import subprocess
import sys
import types

assert sys.version_info >= (3, 12)
assert not sys.flags.optimize, 'Correctness assertions must be enabled.'
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
EXPERIMENT_BASE = 'b490f4f17c97645ce7c6495ca8635bce3d3a1af7'
THEORY_BASE = 'd6e6417f0e055ac0701e5853b9eaf70a316ef709'
THEORY_PRE_REPAIR = 'd07bda9538cfff024c409a259a73559b9b2edd40'
SUPPLIER_BASE = '865f38257ead1b4e857ff915c23b0b04276bd7f3'
BLOBS = {}


def load(path, revision, name):
    raw = subprocess.check_output(['git', '-C', str(ROOT), 'show', revision+':'+path])
    BLOBS[revision+':'+path] = sha256(raw).hexdigest()
    module = types.ModuleType(name)
    module.__file__ = str(ROOT/path)
    exec(compile(raw, 'pinned:'+revision+':'+path, 'exec'), module.__dict__)
    return module


source = load('rro-structural-recovery/exact_checks.py', SUPPLIER_BASE, 'source')
positive = load('rro-positive-certificates/exact_checks.py', EXPERIMENT_BASE, 'positive')
assert source.bilinear.__code__.co_code == positive.source.bilinear.__code__.co_code
A, T, M = source.LEFT, source.RIGHT, source.M
VECTORS = tuple(product(range(3), repeat=3))
ZERO = ((0, 0), (0, 0, 0), (0, 0), 0)


def add(a, b):
    return tuple((x+y) % 3 for x, y in zip(a, b))


def scale(n, a):
    return tuple(n*x % 3 for x in a)


def heap(a, b, c):
    return add(add(a, scale(2, b)), c)


def mul(i, x, y):
    v, w, z, e = x
    vv, ww, zz, ee = y
    ww, zz = scale((-1)**e, ww), scale((-1)**e, zz)
    return add(v, vv), add(w, ww), add(add(z, zz), source.bilinear(i, v, ww)), (e+ee) % 2


def power(i, x, n):
    ans = ZERO
    for _ in range(n):
        ans = mul(i, ans, x)
    return ans


def inverse(i, x):
    v, w, z, e = x
    return scale(2, v), scale(-(-1)**e, w), scale((-1)**e, add(scale(2, z), source.bilinear(i, v, w))), e


def group_source(i):
    """Rebuild actual S_i power pairs, then compare the literal supplier."""
    ai, ti = {a: j for j, a in enumerate(A)}, {t: j for j, t in enumerate(T)}
    all_elements = tuple(product(source.Z, source.W, source.Z, (0, 1)))
    edges, inverse_edges, orders = {}, {}, Counter()
    for x in all_elements:
        xi = inverse(i, x)
        assert mul(i, x, xi) == mul(i, xi, x) == ZERO
        assert inverse(i, xi) == x
        o = next(n for n in (1, 2, 3, 6) if power(i, x, n) == ZERO)
        orders[o] += 1
        if o != 6:
            continue
        a, t = power(i, x, 2), power(i, x, 3)
        assert a[1] == (0, 0, 0) and a[3] == 0 and t[0] == (0, 0) and t[3] == 1
        pair = ai[a[0], a[2]], ti[t[1], t[2]]
        assert pair not in edges
        edges[pair] = x
        aa, tt = power(i, xi, 2), power(i, xi, 3)
        inverse_edges[pair] = ai[aa[0], aa[2]], ti[tt[1], tt[2]]
        aa0, tt0 = A[pair[0]], T[pair[1]]
        assert x == (scale(2, aa0[0]), tt0[0], add(tt0[1], aa0[1]), 1)
        assert inverse_edges[pair] == (ai[scale(2, aa0[0]), scale(2, aa0[1])], pair[1])
    E = set(edges)
    assert E == positive.literal_edges(i) and len(E) == M == 1944
    assert orders == {1: 1, 2: 243, 3: 2186, 6: 1944}
    assert Counter(a for a, _ in E) == {a: 27 for a in range(72)}
    assert Counter(t for _, t in E) == {t: 8 for t in range(243)}
    assert Counter(inverse_edges.values()) == {e: 1 for e in E}
    # Thus a uniform unread edge has a uniform inverse edge in BOTH sources.
    return E, inverse_edges, {'group_elements': len(all_elements), 'orders': dict(sorted(orders.items())),
                             'paired_power_bijection': len(edges), 'uniform_inverse_edge_counts': [1, M]}


def transport(E, lm, rm):
    return {(lm[a], rm[t]) for a, t in E}


def quotient(E):
    rows, cols = positive.index(E)
    l, r = defaultdict(set), defaultdict(set)
    for a in sorted(rows):
        l[frozenset(rows[a])].add(a)
    for t in sorted(cols):
        r[frozenset(cols[t])].add(t)
    left = tuple(sorted((frozenset(c) for c in l.values()), key=lambda c: tuple(sorted(c))))
    right = tuple(sorted((frozenset(c) for c in r.values()), key=lambda c: tuple(sorted(c))))
    li = {a: j for j, c in enumerate(left) for a in c}
    ri = {t: j for j, c in enumerate(right) for t in c}
    lines = tuple(frozenset(ri[t] for t in rows[min(c)]) for c in left)
    incidences = {(li[a], ri[t]) for a, t in E}
    lifted = {(a, t) for l, r in incidences for a in left[l] for t in right[r]}
    assert lifted == E
    return left, right, lines


class Geometry:
    """Consumes only acquired labels and certified literal incidence."""
    def __init__(self, E):
        self.E = E
        self.left, self.right, self.lines = quotient(E)
        if not (len(self.left) == 36 and len(self.right) == 27 and
                all(len(c) == 2 for c in self.left) and
                all(len(c) == 9 for c in self.right) and
                all(len(l) == 3 for l in self.lines)):
            raise ValueError('G0 literal twin capacities/incidence required')
        graph = {j: {j} for j in range(36)}
        planes, histogram = set(), Counter()
        for j, k in combinations(range(36), 2):
            a, b = self.lines[j], self.lines[k]
            if a & b:
                continue
            transversals = [l for l in self.lines if l & a and l & b]
            histogram[len(transversals)] += 1
            if len(transversals) == 3:
                assert all(not x & y for x, y in combinations(transversals, 2))
                plane = frozenset().union(*transversals)
                assert len(plane) == 9 and a | b <= plane
                planes.add(plane)
                graph[j].add(k)
                graph[k].add(j)
        remaining, directions = set(graph), []
        while remaining:
            c = {min(remaining)}
            while True:
                expanded = c | set().union(*(graph[j] for j in c))
                if expanded == c:
                    break
                c = expanded
            directions.append(frozenset(c))
            remaining -= c
        self.directions = tuple(directions)
        self.planes = tuple(sorted(planes, key=lambda p: tuple(sorted(p))))
        self.histogram = dict(sorted(histogram.items()))
        assert len(directions) == 4 and all(len(d) == 9 for d in directions)
        self.through = {}
        for d, cls in enumerate(directions):
            for j in cls:
                for p in self.lines[j]:
                    assert (p, d) not in self.through
                    self.through[p, d] = self.lines[j]
        assert len(self.through) == 27*4 and len(planes) == 18
        self.plane_through = {}
        for p in self.planes:
            dirs = tuple(d for d, cls in enumerate(directions)
                         if sum(self.lines[j] <= p for j in cls) == 3)
            assert len(dirs) == 2
            for point in p:
                assert (point, dirs) not in self.plane_through
                self.plane_through[point, dirs] = p
        assert len(self.plane_through) == 27*6
        # These are acquired attachments, not the hidden maps A->X, T->Y.
        self.left_of = {a: j for j, cls in enumerate(self.left) for a in cls}
        self.point_of = {t: j for j, cls in enumerate(self.right) for t in cls}
        self.chart = next(self.charts())[3]
        self.coords = {p: c for c, p in self.chart.items()}
        self.H = {(p, q, r): self.chart[heap(self.coords[p], self.coords[q], self.coords[r])]
                  for p, q, r in product(range(27), repeat=3)}

    def lift(self):
        """Reverse the quotient using only retained acquired fibres."""
        return {(a, t) for j, line in enumerate(self.lines) for p in line
                for a in self.left[j] for t in self.right[p]}

    def H_labels(self, x, y, z):
        """A point is represented by its whole acquired right-label fibre."""
        try:
            p, q, r = (self.point_of[t] for t in (x, y, z))
        except KeyError as exc:
            raise ValueError('Acquired right labels required') from exc
        return self.right[self.H[p, q, r]]

    def fixed_frame_chart(self):
        """An incidence-computed chart with marks in the declared public P0.

        The selected-axis chart sends intrinsic classes 0,1,2 to coordinate
        axes. Rescale their images in P0 so the fourth class maps to P0[3].
        This is a public representation choice, not the actual hidden chart.
        """
        P0 = ((0, 0, 1), (1, 0, 0), (1, 2, 1), (1, 1, 1))
        p, q, _ = sorted(self.lines[min(self.directions[3])])
        fourth = add(self.coords[q], scale(2, self.coords[p]))
        assert all(fourth)

        def image(v, factors):
            out = (0, 0, 0)
            for x, n, d in zip(v, factors, P0[:3]):
                out = add(out, scale(x*n, d))
            return out

        factors = next(f for f in product((1, 2), repeat=3)
                       if image(fourth, f) in (P0[3], scale(2, P0[3])))
        coords = {p: image(v, factors) for p, v in self.coords.items()}
        assert len(set(coords.values())) == 27
        return coords, P0

    def full_incidence(self):
        """Full lines, planes and parallel classes from H alone."""
        points = range(27)
        lines = {frozenset((p, q, self.H[q, p, q])) for p, q in combinations(points, 2)}
        planes = set()
        for p, q, r in combinations(points, 3):
            if r in {p, q, self.H[q, p, q]}:
                continue
            u, v = (p, q, self.H[q, p, q]), (p, r, self.H[r, p, r])
            plane = frozenset(self.H[a, p, b] for a, b in product(u, v))
            assert len(plane) == 9
            planes.add(plane)
        # Translate each line through every point, using an origin on that line.
        parallel = {frozenset(frozenset(self.H[p, min(line), a] for a in line)
                              for p in points) for line in lines}
        assert set().union(*parallel) == lines
        return lines, planes, parallel

    @staticmethod
    def meet(a, b):
        c = a & b
        assert len(c) == 1
        return next(iter(c))

    def charts(self):
        for o in range(27):
            for dirs in permutations(range(4), 3):
                axis_choices = [tuple(sorted(self.through[o, d]-{o})) for d in dirs]
                for signs in product((0, 1), repeat=3):
                    axes = [(o, choices[s], choices[1-s]) for choices, s in zip(axis_choices, signs)]
                    xy = {(x, y): self.meet(self.through[axes[0][x], dirs[1]],
                                            self.through[axes[1][y], dirs[0]])
                          for x, y in product(range(3), repeat=2)}
                    plane_dirs = tuple(sorted(dirs[:2]))
                    chart = {(x, y, z): self.meet(self.through[xy[x, y], dirs[2]],
                                                  self.plane_through[axes[2][z], plane_dirs])
                             for x, y, z in VECTORS}
                    assert len(set(chart.values())) == 27
                    yield o, dirs, signs, chart

    def inverse_label(self, y):
        a, t = y
        cls = next((c for c in self.left if a in c), None)
        if cls is None or y not in self.E:
            raise ValueError('An acquired actual input edge is required')
        return next(iter(cls-{a})), t


def relabelings():
    yield 'identity', tuple(range(72)), tuple(range(243))
    yield 'reverse_ambient', tuple(127-a for a in range(72)), tuple(255-t for t in range(243))
    yield 'odd_ambient_permutations', tuple((37*a+11) % 128 for a in range(72)), tuple((73*t+19) % 256 for t in range(243))
    l = sorted(range(128), key=lambda x: sha256(('left:'+str(x)).encode()).digest())[:72]
    r = sorted(range(256), key=lambda x: sha256(('right:'+str(x)).encode()).digest())[:243]
    yield 'sha256_order', tuple(l), tuple(r)


def disjoint_line_covers(lines):
    """Literal disjoint nine-line cover predicate on the actual 27 points."""
    through = {p: [j for j, line in enumerate(lines) if p in line] for p in range(27)}
    covers = set()

    def visit(uncovered, chosen):
        if not uncovered:
            assert len(chosen) == 9
            covers.add(frozenset(chosen))
            return
        p = min(uncovered)
        for j in through[p]:
            if lines[j] <= uncovered:
                visit(uncovered-lines[j], chosen+(j,))

    visit(frozenset(range(27)), ())
    return covers


def exact_cover_check(geom):
    """Enumerate the disjoint-cover object, rather than assuming parallelism."""
    covers = disjoint_line_covers(geom.lines)
    selected = set(geom.directions)
    assert selected <= covers and len(covers) == 40 and len(selected) == 4
    mixed = covers-selected
    assert len(mixed) == 36
    profiles = Counter()
    for cover in covers:
        assert len(cover) == 9 and set().union(*(geom.lines[j] for j in cover)) == set(range(27))
        assert all(not geom.lines[j] & geom.lines[k] for j, k in combinations(cover, 2))
        profiles[tuple(sorted(len(cover & d) for d in geom.directions if cover & d))] += 1
    witness = min(mixed, key=lambda x: tuple(sorted(x)))
    assert len({d for d, cls in enumerate(geom.directions) if cls & witness}) > 1
    return {'all_disjoint_nine_line_covers': len(covers), 'three_transversal_components': len(selected),
            'extra_mixed_covers': len(mixed),
            'direction_occupancy_profiles': {str(k): v for k, v in sorted(profiles.items())},
            'mixed_cover_line_indices': sorted(witness),
            'mixed_cover_point_blocks': [sorted(geom.lines[j]) for j in sorted(witness)],
            'old_cover_definition_equals_selected_classes': covers == selected}


def direction_check(geom, coordinates):
    def projective(v):
        pivot = next(x for x in v if x)
        return scale(pivot, v)

    def marks(coords):
        ans = []
        for cls in geom.directions:
            directions = set()
            for j in cls:
                p, q, r = sorted(geom.lines[j])
                d = projective(add(coords[q], scale(2, coords[p])))
                assert add(coords[p], scale(2, add(coords[q], scale(2, coords[p])))) == coords[r]
                directions.add(d)
            assert len(directions) == 1
            ans.append(next(iter(directions)))
        return ans

    source_marks, chart_marks = marks(coordinates), marks(geom.coords)
    declared_coords, P0 = geom.fixed_frame_chart()
    declared_marks = marks(declared_coords)
    assert tuple(declared_marks) == P0
    assert all(declared_coords[h] == heap(declared_coords[p], declared_coords[q], declared_coords[r])
               for (p, q, r), h in geom.H.items())
    assert set(source_marks) == {(0, 0, 1), (1, 0, 0), (1, 2, 1), (1, 1, 1)}
    assert {(1, 0, 0), (0, 1, 0), (0, 0, 1)} <= set(chart_marks)
    fourth = next(d for d in chart_marks if sum(x != 0 for x in d) == 3)
    assert all(fourth)
    determinants = []
    for a, b, c in combinations(source_marks, 3):
        det = (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
               +a[2]*(b[0]*c[1]-b[1]*c[0])) % 3
        assert det
        determinants.append(det)
    return {'source_P0_marks_by_intrinsic_class': source_marks,
            'canonical_chart_marks_by_intrinsic_class': chart_marks,
            'source_frame_triple_determinants': determinants,
            'source_marks_are_chart_marks': set(source_marks) == set(chart_marks),
            'declared_P0_chart_marks_by_intrinsic_class': declared_marks,
            'declared_P0_chart_points': [declared_coords[p] for p in range(27)],
            'declared_P0_chart_H_triples_checked': len(geom.H),
            'chart_contract': 'three selected axes, either orientation; fourth has three nonzero components'}


def geometry_check(E, lm, rm):
    labelled = transport(E, lm, rm)
    D = transport(positive.g0_witness(E), lm, rm)
    completed, rounds = positive.g0_complete(D)
    assert completed == labelled
    geom = Geometry(completed)
    invr = {label: T[t][0] for t, label in enumerate(rm)}
    coordinates = {j: invr[min(c)] for j, c in enumerate(geom.right)}
    assert all(len({invr[label] for label in c}) == 1 for c in geom.right)
    chart_count, point_checks = 0, 0
    e0, e1, e2 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    for origin, dirs, signs, chart in geom.charts():
        offset = coordinates[origin]
        basis = [add(coordinates[chart[e]], scale(2, offset)) for e in (e0, e1, e2)]
        assert len({add(add(scale(x, basis[0]), scale(y, basis[1])), scale(z, basis[2]))
                    for x, y, z in VECTORS}) == 27
        for c, p in chart.items():
            predicted = offset
            for x, b in zip(c, basis):
                predicted = add(predicted, scale(x, b))
            assert coordinates[p] == predicted
            point_checks += 1
        chart_count += 1
    assert chart_count == 27*24*8 == 5184
    # Every chart is affine, so its H is the same by cancellation of offset
    # in (b+Au)-(b+Av)+(b+Aw). Enumerate all triples in one chart as well.
    for (p, q, r), result in geom.H.items():
        assert coordinates[result] == heap(coordinates[p], coordinates[q], coordinates[r])
    full_lines, full_planes, intrinsic_parallel = geom.full_incidence()
    assert len(full_lines) == 117 and set(geom.lines) <= full_lines
    # Independent chart-space construction checks the H-only plane definition.
    coordinate_planes = set()
    for p, q, r in combinations(range(27), 3):
        if r in {p, q, geom.H[q, p, q]}:
            continue
        generated = frozenset(geom.chart[add(geom.coords[p], add(
            scale(x, add(geom.coords[q], scale(2, geom.coords[p]))),
            scale(y, add(geom.coords[r], scale(2, geom.coords[p])))))]
            for x, y in product(range(3), repeat=2))
        coordinate_planes.add(generated)
    assert full_planes == coordinate_planes
    assert len(full_planes) == 39 and set(geom.planes) <= full_planes
    assert len(intrinsic_parallel) == 13 and set(map(len, intrinsic_parallel)) == {9}
    assert all(set().union(*cls) == set(range(27)) for cls in intrinsic_parallel)
    assert all(sum({p, q} <= line for line in full_lines) == 1
               for p, q in combinations(range(27), 2))
    assert all(sum(p in plane for plane in full_planes) == 13 for p in range(27))
    full_directions = Counter()
    for line in full_lines:
        p, q, _ = sorted(line)
        delta = add(geom.coords[q], scale(2, geom.coords[p]))
        pivot = next(x for x in delta if x)
        full_directions[scale(pivot, delta)] += 1
    assert len(full_directions) == 13 and set(full_directions.values()) == {9}
    assert all(sum(p in line for line in full_lines) == 13 for p in range(27))
    assert all(sum(line <= plane for line in full_lines) == 12 for plane in full_planes)
    assert all(sum(line <= plane for plane in full_planes) == 4 for line in full_lines)
    # Exact reverse, including literal labels/capacities and the four marks.
    reverse = geom.lift()
    assert all(geom.H_labels(min(geom.right[p]), min(geom.right[q]), min(geom.right[r]))
               == geom.right[h] for (p, q, r), h in geom.H.items())
    assert reverse == labelled
    assert sum((lm[a], rm[t]) in reverse for a in range(72) for t in range(243)) == M
    for a, t in E:
        v, s = A[a]
        ai = A.index((scale(2, v), scale(2, s)))
        assert geom.inverse_label((lm[a], rm[t])) == (lm[ai], rm[t])
        assert geom.inverse_label(geom.inverse_label((lm[a], rm[t]))) == (lm[a], rm[t])
    for bad in ((next(x for x in range(128) if x not in lm), rm[0]),
                next((a, t) for a in lm for t in rm if (a, t) not in labelled)):
        try:
            geom.inverse_label(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid input accepted')
    return geom, {'certificate_pairs': len(D), 'completion_round_sizes': rounds,
                  'literal_endpoint_pairs_checked': 72*243, 'reverse_edges': len(reverse),
                  'left_twins': [36, 2], 'right_twins': [27, 9], 'marked_incidences': 108,
                  'disjoint_line_transversal_counts': geom.histogram,
                  'parallel_classes': [4, 9], 'marked_grid_planes': len(geom.planes),
                  'charts': chart_count, 'chart_point_comparisons': point_checks,
                  'H_triples_enumerated': len(geom.H), 'H_chart_triples_implied_by_affinity': chart_count*len(geom.H),
                  'full_affine_lines': len(full_lines), 'full_affine_planes': len(full_planes),
                  'full_affine_parallel_classes': len(full_directions),
                  'full_lines_per_point': 13, 'full_lines_per_plane': 12, 'full_planes_per_line': 4,
                  'literal_inverse_edges': len(E),
                  'observable_left_attachments': len(geom.left_of),
                  'observable_right_attachments': len(geom.point_of),
                  'H_fibre_outputs_checked': len(geom.H),
                  'finite_structure_sizes': {
                      'acquired_relation_edges': len(geom.E),
                      'left_label_class_attachments': len(geom.left_of),
                      'right_label_class_attachments': len(geom.point_of),
                      'quotient_incidences': sum(map(len, geom.lines)),
                      'one_chart_points': len(geom.chart), 'stored_H_entries': len(geom.H),
                      'full_line_point_memberships': sum(map(len, full_lines)),
                      'full_plane_point_memberships': sum(map(len, full_planes)),
                      'full_parallel_line_memberships': sum(map(len, intrinsic_parallel)),
                      'charts_enumerated_streaming': chart_count},
                  'exact_cover_discriminator': exact_cover_check(geom),
                  'chart_direction_correspondence': direction_check(geom, coordinates)}


def perfect_matchings(items):
    if not items:
        yield {}
        return
    a, *rest = items
    for b in rest:
        for tail in perfect_matchings([x for x in rest if x != b]):
            yield tail | {a: b, b: a}


def g1_check(E, inverses, lm, rm):
    labelled = transport(E, lm, rm)
    D = transport(positive.g1_witness(E), lm, rm)
    assert positive.g1_complete(D) == labelled
    left, right, lines = quotient(labelled)
    assert len(left) == len(right) == 9
    assert all(len(c) == 8 for c in left) and all(len(c) == 27 for c in right)
    assert all(len(l) == 1 for l in lines)
    try:
        Geometry(labelled)
    except ValueError:
        pass
    else:
        raise AssertionError('G1 admitted to G0 recovery')
    # All possible swaps of the inverse partner with another twin, fixing a.
    # Full prepared array contains each actual edge once and repeats the input.
    # Each transformed hidden array uses the same once-changed full map.
    literal_classes, literal_right, _ = quotient(E)
    for cls in literal_right:
        assert len({T[t][0][:2] for t in cls}) == 1
        assert len({(T[t][0][2], T[t][1]) for t in cls}) == 27
    assert all(len({A[a][0] for a in cls}) == 8 for cls in literal_classes)
    witnesses, array_checks, example = 0, 0, None
    edge_order = sorted(E)
    for cls in literal_classes:
        for a in sorted(cls):
            inva = next(inverses[a, t][0] for t in range(243) if (a, t) in E)
            for b in sorted(cls-{a, inva}):
                sigma = list(range(72))
                sigma[inva], sigma[b] = b, inva
                altered_lm = tuple(lm[sigma[x]] for x in range(72))
                assert transport(E, altered_lm, rm) == labelled
                t = min(t for aa, t in E if aa == a)
                prepared = edge_order+[(a, t)]
                for aa, tt in prepared:
                    assert (sigma[aa], tt) in E
                    assert (altered_lm[sigma[aa]], rm[tt]) == (lm[aa], rm[tt])
                assert (altered_lm[a], rm[t]) == (lm[a], rm[t])
                first, second = (lm[inva], rm[t]), (altered_lm[inva], rm[t])
                assert first != second
                if example is None:
                    example = {'input_literal': [A[a], T[t]], 'input_labels': [lm[a], rm[t]],
                               'first_inverse_labels': first, 'second_inverse_labels': second,
                               'fixed_input_left_index': a, 'swapped_left_indices': [inva, b],
                               'swapped_literal_endpoints': [A[inva], A[b]],
                               'prepared_N': len(prepared), 'prepared_edge_order': 'sorted literal edge indices, then input',
                               'same_array_sha256': sha256(bytes(bit for aa, tt in prepared
                                     for bit in positive.bits_of((lm[aa], rm[tt])))).hexdigest()}
                witnesses += 1
                array_checks += len(prepared)
    assert witnesses == 72*6 == 432
    # Conditional on certified E and any supported full prepared label array,
    # within-component relabeling leaves likelihood M^-N fixed. The inherited
    # source inversion is a uniformly transported perfect matching on 8 labels.
    matchings = list(perfect_matchings(list(sorted(left[0]))))
    assert len(matchings) == 105
    x, y = sorted(left[0])[:2]
    one = Counter(p[x] for p in matchings)
    assert set(one.values()) == {15} and len(one) == 7
    assert sum(p[x] == y and p[y] == x for p in matchings) == 15
    assert sum(p[x] == p[x] for p in matchings) == 105
    fixed_input_partner = [p for p in matchings if p[x] == y]
    assert len(fixed_input_partner) == 15 and {p[x] for p in fixed_input_partner} == {y}
    assert len({tuple(sorted(p.items())) for p in fixed_input_partner}) > 1
    right_component = next(l for cls, l in zip(left, lines) if x in cls)
    right_labels = right[next(iter(right_component))]
    partial_inverse = Counter((p[x], t) for p in matchings for t in right_labels)
    assert len(partial_inverse) == 7*27 and set(partial_inverse.values()) == {15}
    assert F(15, len(matchings)*len(right_labels)) == F(1, 189)
    # Count all local bijections as an independent weighting check, not 72!/243!.
    actual = sorted(literal_classes[0])
    named_partner = {a: next(inverses[a, t][0] for t in range(243) if (a, t) in E) for a in actual}
    weights = Counter()
    for perm in permutations(actual):
        to_label = dict(zip(perm, range(8)))
        weights[tuple(to_label[named_partner[a]] for a in perm)] += 1
    assert len(weights) == 105 and set(weights.values()) == {384} and sum(weights.values()) == 40320
    # Anonymous G1 component equality cannot choose an affine F3 plane: even
    # swapping two components preserves E but fails every assigned affine heap.
    points = tuple(product(range(3), repeat=2))
    nonaffine = 0
    for j, k in combinations(range(9), 2):
        perm = list(points)
        perm[j], perm[k] = perm[k], perm[j]
        f = dict(zip(points, perm))
        assert any(f[heap(a, b, c)] != heap(f[a], f[b], f[c]) for a, b, c in product(points, repeat=3))
        nonaffine += 1
    assert len(positive.components(labelled)) == 9
    assert len({frozenset(line) for line in lines}) == 9
    assert sum(map(len, lines)) == 9
    return {'certificate_pairs': len(D), 'literal_components_K8_27': 9,
            'twin_quotient_components_K1_1': 9, 'twin_quotient_edges': sum(map(len, lines)), 'components': 9, 'left_twins': [9, 8], 'right_twins': [9, 27],
            'all_literal_endpoint_pairs': 72*243, 'free_w3_z_choices_per_right_class': 27,
            'same_observation_witnesses': witnesses,
            'same_prepared_slot_comparisons': array_checks, 'example': example,
            'local_bijections_enumerated': 40320, 'inverse_matchings': 105, 'bijections_per_matching': 384,
            'point_input_partner_probability': '1/7',
            'fixed_input_partner_datum': {'input': x, 'retained_partner': y,
                                         'compatible_whole_class_matchings': len(fixed_input_partner),
                                         'fixed_input_output_support': len({p[x] for p in fixed_input_partner})},
            'left_label_only_inverse_edge_support': len(partial_inverse),
            'left_label_only_inverse_edge_mass': '1/189', 'mutual_partner_joint_probability': '1/7',
            'independently_redrawn_mutual_probability': '1/49', 'repeated_input_agreement': '1',
            'nonaffine_component_transpositions': nonaffine}


def distribution(values):
    c = Counter(values)
    return {x: F(n, c.total()) for x, n in c.items()}


def conditional_check(geom, inverses, lm, rm):
    E = geom.E
    reference = {(lm[a], rm[t]): (lm[b], rm[u]) for (a, t), (b, u) in inverses.items()}
    complete = [dict(enumerate(positive.bits_of(y))) for y in sorted(E)]
    masks = [{}]+[{h: bit} for h in range(15) for bit in (0, 1)]
    y = min(E)
    for fields in (range(7), range(7, 15), (0, 3, 7, 11), (1, 2, 8, 14)):
        masks.append({h: positive.bits_of(y)[h] for h in fields})
    masks += complete
    checked, domain_mass_sum = 0, F(0)
    for mask in masks:
        if len(mask) == 15:
            pair = (sum(mask[h] << (6-h) for h in range(7)),
                    sum(mask[h+7] << (7-h) for h in range(8)))
            support = {pair} if pair in E else set()
        else:
            support = positive.conditioned(E, mask)
        assert support
        recovered = distribution(geom.inverse_label(y) for y in support)
        actual = distribution(reference[y] for y in support)
        assert recovered == actual and sum(recovered.values()) == 1
        assert len(recovered) == len(support)
        # Direct literal-edge joint event mass equals conditional*likelihood.
        likelihood = F(len(support), M)
        assert all(prob*likelihood == F(1, M) for prob in actual.values())
        domain_mass_sum += likelihood
        checked += 1
    assert checked == 1+30+4+1944
    # A repeated slot is one hidden edge, while distinct original slots are iid
    # conditional on E. Check exact joint cylinders on partial/complete/unread.
    a = min(E)[0]
    ma = {h: positive.bits_of((a, 0))[h] for h in range(7)}
    support = positive.conditioned(E, ma)
    assert len(support) == 27
    same_slot = distribution((geom.inverse_label(y), geom.inverse_label(y)) for y in support)
    assert all(x == y for x, y in same_slot) and len(same_slot) == 27
    distinct = {(geom.inverse_label(x), geom.inverse_label(y)): F(1, 27*27)
                for x, y in product(support, repeat=2)}
    assert sum(distinct.values()) == 1 and sum(p for (x, y), p in distinct.items() if x == y) == F(1, 27)
    completed = min(support)
    assert positive.conditioned(E, dict(enumerate(positive.bits_of(completed)))) == {completed}
    # Sequentially acquiring a partial old slot telescopes, including rereads.
    old_mask = ma.copy()
    telescoped = F(1)
    steps = []
    for h in (7, 7, 8, 9, 10, 11, 12, 13, 14):
        before = positive.conditioned(E, old_mask)
        bit = positive.bits_of(completed)[h]
        after = positive.conditioned(E, old_mask | {h: bit})
        p = F(len(after), len(before))
        steps.append(str(p))
        telescoped *= p
        old_mask[h] = bit
    assert telescoped == F(1, 27) and len(positive.conditioned(E, old_mask)) == 1
    return {'positive_masks': checked, 'single_bit_rectangles': 30, 'representative_mixed_or_field_rectangles': 4,
            'complete_edge_rectangles': 1944, 'sum_likelihoods_over_declared_rectangles': str(domain_mass_sum),
            'same_old_slot_inverse_pairs': 27, 'distinct_partial_slot_inverse_pairs': 729,
            'same_slot_agreement': '1', 'distinct_partial_slot_agreement': '1/27',
            'partial_completion_and_repeat_probabilities': steps, 'telescoped_completion_probability': str(telescoped)}


def geometry_law_check(geom):
    """Push the existing indexed edge law through intrinsic point/H readouts.

    Count literal lifts per right twin class, so no uniform quotient sampling
    assumption replaces the actual source law. Repeated indices use one draw.
    """
    E = geom.E
    point_of = {t: p for p, cls in enumerate(geom.right) for t in cls}
    coords = geom.coords
    y = min(E)
    partial_mask = {h: positive.bits_of(y)[h] for h in range(7)}
    partial = positive.conditioned(E, partial_mask)
    complete = {max(E)}
    scenarios = [('three_unread', [E, E, E], (0, 1, 2)),
                 ('partial_complete_unread', [partial, complete, E], (0, 1, 2)),
                 ('same_partial_complete_same_partial', [partial, complete], (0, 1, 0)),
                 ('same_partial_three_times', [partial], (0, 0, 0))]
    result = {}
    for name, domains, indices in scenarios:
        masses = [Counter(point_of[t] for a, t in domain) for domain in domains]
        denominator = 1
        for domain in domains:
            denominator *= len(domain)
        output, total, assignments = Counter(), 0, 0
        for points in product(*(sorted(m) for m in masses)):
            weight = 1
            for p, m in zip(points, masses):
                weight *= m[p]
            p, q, r = (points[j] for j in indices)
            h = geom.H[p, q, r]
            assert coords[h] == heap(coords[p], coords[q], coords[r])
            output[h] += weight
            total += weight
            assignments += 1
        assert total == denominator
        if name in ('three_unread', 'partial_complete_unread'):
            assert len(output) == 27 and set(output.values()) == {denominator//27}
        else:
            assert len(output) == 3 and set(output.values()) == {denominator//3}
        result[name] = {'original_index_pattern': indices, 'literal_domain_sizes': list(map(len, domains)),
                        'point_assignments_enumerated': assignments, 'literal_edge_array_count': denominator,
                        'output_points': len(output), 'equal_output_mass': str(F(next(iter(output.values())), denominator))}
    masses = Counter(point_of[t] for a, t in partial)
    assert len(masses) == 3 and set(masses.values()) == {9}
    correct_joint = {(p, p): F(n, len(partial)) for p, n in masses.items()}
    independent_joint = {(p, q): F(n*m, len(partial)**2) for p, n in masses.items() for q, m in masses.items()}
    assert correct_joint != independent_joint
    assert sum(v for (p, q), v in correct_joint.items() if p == q) == 1
    assert sum(v for (p, q), v in independent_joint.items() if p == q) == F(1, 3)
    result['repeated_point_input_agreement'] = {'indexed': '1', 'incorrect_independent_redraw': '1/3'}
    return result


def native_acquisition(E, witness):
    """Reuse the existing tree/read interpreter; no new source port or draw."""
    q = 576  # Public scan budget, identical for both source branches.
    prepared = sorted(witness)+[min(witness)]*(q-len(witness))+[min(E), max(E), min(E), max(E)]
    N = len(prepared)
    bits = tuple(b for y in prepared for b in positive.bits_of(y))
    term = positive.tree(bits)
    cache = [{} for _ in prepared]
    operands = []

    def query(j, h):
        b = 15*j+h
        word = (1,)*b+(0, 0)
        reply = positive.read(term, word)
        assert reply in (positive.BETA, positive.BRANCH)
        value = int(reply == positive.BRANCH)
        assert h not in cache[j] or cache[j][h] == value
        cache[j][h] = value
        operands.append(word)
        return value

    for j in range(q):
        for h in range(15):
            query(j, h)
    acquired = {(sum(cache[j][h] << (6-h) for h in range(7)),
                 sum(cache[j][h+7] << (7-h) for h in range(8))) for j in range(q)}
    assert acquired == witness
    candidate0, candidate1 = positive.g0_complete(acquired), positive.g1_complete(acquired)
    assert (candidate0 is None) != (candidate1 is None)
    completed_E = candidate0[0] if candidate0 else candidate1
    assert completed_E == E
    recover = Geometry(completed_E) if candidate0 else None
    partial, complete, unread, repeated = q, q+1, q+2, q+3
    for h in range(7):
        query(partial, h)
    for h in range(15):
        query(complete, h)
    partial_initial_cache = [sum(1 << h for h in cache[partial]),
                             sum(bit << h for h, bit in cache[partial].items())]
    before = positive.conditioned(E, cache[partial])
    assert len(before) == 27 and cache[unread] == cache[repeated] == {}
    for h in (0, 0):
        query(complete, h)
    if recover:
        output_law = distribution(recover.inverse_label(y) for y in before)
        assert len(output_law) == 27
        for h in range(7, 15):
            query(partial, h)
        pair = (sum(cache[partial][h] << (6-h) for h in range(7)),
                sum(cache[partial][h+7] << (7-h) for h in range(8)))
        output = recover.inverse_label(pair)
        assert pair == prepared[partial] and output in output_law
        point_of = {t: j for j, cls in enumerate(recover.right) for t in cls}
        p, qpoint = point_of[pair[1]], point_of[prepared[complete][1]]
        point_output = sorted(recover.H_labels(pair[1], prepared[complete][1], pair[1]))
        assert recover.lift() == completed_E
    else:
        output, output_law, point_output = None, {}, None
    # Fixed root and deep absence are invocations without cache information.
    for word, expected in (((), positive.BRANCH), ((0, 0, 0, 0, 0), positive.ABSENT)):
        assert positive.read(term, word) == expected
        operands.append(word)
    assert cache[unread] == cache[repeated] == {}
    return {'N': N, 'scan_slots': q, 'scan_Read_invocations': 15*q,
            'scan_literal_operand_sum': (15*q)*(15*q+3)//2, 'scan_literal_operand_max': 15*q+1,
            'all_Read_invocations': len(operands), 'all_distinct_Read_addresses': len(set(operands)),
            'all_literal_operand_sum': sum(map(len, operands)), 'all_literal_operand_max': max(map(len, operands)),
            'preparation_events': 75*N+4, 'paid_preparer_draws': N,
            'partial_original_index_1based': partial+1, 'complete_original_index_1based': complete+1,
            'partial_initial_mask_value': partial_initial_cache,
            'indexed_cache_encoding': 'original sample j is row j-1; [mask,value], bit h is local coordinate h (0..14)',
            'indexed_cache_mask_value_rows': [[sum(1 << h for h in row),
                                               sum(bit << h for h, bit in row.items())] for row in cache],
            'unread_original_indices_1based': [unread+1, repeated+1],
            'partial_output_support_before_completion': len(output_law) if recover else None,
            'pointwise_inverse_recoverable': bool(recover), 'realized_inverse_output_labels': output,
            'H_partial_complete_same_partial_output_right_labels': point_output,
            'original_array_sha256': sha256(bytes(bits)).hexdigest()}


def ambiguity_check(E, lm, rm, geom):
    labelled = transport(E, lm, rm)
    edge_order = sorted(E)
    translations = 0
    for delta in VECTORS[1:]:
        aa = {a: A.index((v, add(s, source.bilinear(0, v, delta)))) for a, (v, s) in enumerate(A)}
        tt = {t: T.index((add(w, delta), z)) for t, (w, z) in enumerate(T)}
        assert {(aa[a], tt[t]) for a, t in E} == E
        newlm, newrm = [None]*72, [None]*243
        for a in range(72):
            newlm[aa[a]] = lm[a]
        for t in range(243):
            newrm[tt[t]] = rm[t]
        assert transport(E, newlm, newrm) == labelled
        assert all((newlm[aa[a]], newrm[tt[t]]) == (lm[a], rm[t]) for a, t in edge_order)
        assert add((0, 0, 0), delta) != (0, 0, 0)
        translations += 1
    # Same quotient point, nine different literal right labels: retaining
    # associations matters for derived point output and literal edge reversal.
    cls = geom.right[0]
    assert len(cls) == 9
    swapz = list(rm)
    x, y = sorted(cls)[:2]
    ix, iy = rm.index(x), rm.index(y)
    assert T[ix][0] == T[iy][0] and T[ix][1] != T[iy][1]
    swapz[ix], swapz[iy] = swapz[iy], swapz[ix]
    assert transport(E, lm, swapz) == labelled
    # E/closure forgets old index assignment. The same two completed labels
    # have identical M^-2 likelihood at every compatible g, but replay differs.
    y0, y1 = min(labelled), max(labelled)
    h, hh = [y0, y1, y0, y0], [y1, y0, y0, y0]
    assert set(h[:2]) == set(hh[:2]) and h[0] != hh[0]
    likelihoods = []
    for history in (h, hh):
        likelihood = F(1)
        for pair in history[:2]:
            assert pair in labelled
            likelihood *= F(1, M)
        likelihoods.append(likelihood)
    assert likelihoods[0] == likelihoods[1] == F(1, M)**2
    return {'nonzero_translation_same_observation_witnesses': translations,
            'right_twin_reference_swap': {'literal_right_indices': [ix, iy], 'labels': [x, y],
                                          'source_endpoints': [T[ix], T[iy]]},
            'old_index_memory_witness': {'first_two_completed_pairs': h[:2],
                                         'permuted_first_two_completed_pairs': hh[:2],
                                         'common_completed_likelihood': str(F(1, M)**2),
                                         'replayed_first_labels_differ': True}}


def falling(n, k):
    ans = 1
    for j in range(k):
        ans *= n-j
    return ans if k <= n else 0


def extension_mass(n, bucket_sizes, requirements):
    """Exact marginal of a UNIFORM FULL injection, by equal extension counts.

    Each distinct named endpoint has a bucket constraint. Unmentioned endpoints
    retain (n-d)_(domain_size-d) full extensions, cancelling the full prior.
    Bucket assignment is never a replacement source or an independent redraw.
    """
    counts = Counter(requirements.values())
    numerator = 1
    for b, size in enumerate(bucket_sizes):
        numerator *= falling(size, counts[b])
    return F(numerator, falling(n, len(requirements)))


def left_inverse_projection(E, J, cache, targets, N=4):
    """Full-prior pushforward of declared left-high-bit projections.

    Enumerate actual left endpoint assignments once per DISTINCT original index.
    All actual right choices integrate to degree 27, and all unused slots to 1.
    Full injections are integrated by exact extension counts, including inverse
    endpoints. This is a latent law, not an acquired inversion operation.
    """
    assert N >= 4 and all(0 <= j < N for j in set(cache) | set(targets))
    rows, _ = positive.index(E)
    assert len(rows) == 72 and set(map(len, rows.values())) == {27}
    indices = sorted(set(cache) | set(targets))
    counts = Counter()
    for values in product(range(72), repeat=len(indices)):
        sampled = dict(zip(indices, values))
        for output in product((0, 1), repeat=len(targets)):
            requirements = {}
            restrictions = [(sampled[j], b) for j, b in cache.items()]
            restrictions += [(J[sampled[j]], b) for j, b in zip(targets, output)]
            for a, b in restrictions:
                if a in requirements and requirements[a] != b:
                    break
                requirements[a] = b
            else:
                ones = sum(requirements.values())
                counts[output, len(requirements)-ones, ones] += 1
    joint = Counter()
    for (out, zeros, ones), count in counts.items():
        probability = extension_mass(128, (64, 64), dict(enumerate((0,)*zeros+(1,)*ones)))
        joint[out] += F(count, 72**len(indices))*probability
    evidence = sum(joint.values())
    assert evidence > 0
    law = {out: p/evidence for out, p in joint.items() if p}
    assert sum(law.values()) == 1
    return evidence, law, len(indices), 72**len(indices)


def hypergeom(population, successes, draws):
    masses = {k: F(comb(successes, k)*comb(population-successes, draws-k), comb(population, draws))
              for k in range(max(0, draws-(population-successes)), min(successes, draws)+1)}
    assert sum(masses.values()) == 1
    return masses


def weighted_mean(masses, observable, likelihood=lambda k: 1):
    denominator = sum(p*likelihood(k) for k, p in masses.items())
    assert denominator > 0
    return sum(p*likelihood(k)*observable(k) for k, p in masses.items())/denominator


def anchored_inverse_check(E, J):
    # Condition only on acquired Y_0=(0,0), retaining all possible actual anchor
    # edges and full injections. Y_1 has only its left high bit acquired as 0.
    # Unread slots 2,3 integrate out; right map extensions integrate to 1.
    joint, evidence = F(0), F(0)
    degree = Counter(a for a, t in E)
    for a in range(72):
        for b in range(72):
            # Given lambda(a)=0, the remaining 71 endpoints use 127 labels,
            # with bucket capacities 63/64. J(a) is a distinct named endpoint.
            requirements = {} if b == a else {b: 0}
            evidence += F(degree[a]*degree[b], M*M)*extension_mass(127, (63, 64), requirements)
            requirements[J[a]] = 0
            joint += F(degree[a]*degree[b], M*M)*extension_mass(127, (63, 64), requirements)
    answer = joint/evidence
    return answer, evidence/F(128*256), evidence


def mixed_inverse_H_check(E, J):
    """N=4, full anchor + TWO partial fields, joint inverse/H output event.

    H's quotient output is represented by its whole right-label fibre. The event
    checked is whether that fibre contains the acquired anchor label 0. This is
    a latent G0 pushforward on an uncertified cache, not a controller that can
    execute H without the certified relation. All M^2 literal pairs are counted
    by actual-neighbour intersections, not by a reduced substitute graph.
    """
    rows, _ = positive.index(E)
    wcounts = {a: Counter(T[t][0] for t in rows[a]) for a in rows}
    classes = Counter()
    for a in range(72):
        for b in range(72):
            category = 0 if b == a else 1 if b == J[a] else 2
            equal_t = len(rows[a] & rows[b])
            equal_w = sum(n*wcounts[b][w] for w, n in wcounts[a].items())
            classes[category, 0] += equal_t
            classes[category, 1] += equal_w-equal_t
            classes[category, 2] += len(rows[a])*len(rows[b])-equal_w
    assert sum(classes.values()) == M*M
    assert all((heap(w, u, w) == w) == (u == w) for w, u in product(VECTORS, repeat=2))
    mass = Counter()
    # lambda(a)=0 and rho(t)=0 each have the SAME 1/(128*256) mass
    # at every actual anchor edge. Neither likelihood is inserted a second time.
    for (category, right_category), count in classes.items():
        for inverse_bit in (0, 1):
            # Named endpoints 0=anchor a, 1=J(a), 2=other b.
            req = {1: inverse_bit}
            if category == 1 and inverse_bit != 0:
                continue
            if category != 0:
                req[1 if category == 1 else 2] = 0
            left_p = extension_mass(127, (63, 64), req)
            right_p = F(1) if right_category == 0 else F(127, 255)
            contains_anchor = right_category != 2
            mass[inverse_bit, contains_anchor] += F(count, M*M*128*256)*left_p*right_p
    evidence = sum(mass.values())
    assert evidence > 0
    law = {out: p/evidence for out, p in mass.items() if p}
    assert len(law) == 4 and sum(law.values()) == 1
    inv0 = sum(p for (bit, event), p in law.items() if bit == 0)
    event = sum(p for (bit, event), p in law.items() if event)
    assert law[0, True] != inv0*event
    # Repeated original index means cancellation H(p,q,q)=p. A different
    # wholly unread original slot has 72 literal lifts at each of 27 points.
    point_lifts = Counter(T[t][0] for a, t in E)
    assert set(point_lifts.values()) == {72} and len(point_lifts) == 27
    independent_slot_agreement = {w: F(point_lifts[w], M) for w in VECTORS}
    assert set(independent_slot_agreement.values()) == {F(1, 27)}
    return {'N': 4, 'indexed_cache': {'0': 'complete label pair (0,0)',
                                    '1': 'left coordinate 0 = 0; right coordinate 0 = 0',
                                    '2': 'unread', '3': 'unread'},
            'actual_literal_pair_assignments_counted': sum(classes.values()),
            'actual_literal_pair_assignments_enumerated_individually': False,
            'pair_count_method': 'all 72^2 left pairs; exact actual-right-neighbour intersections and w multiplicities',
            'positive_evidence_mass': str(evidence),
            'joint_output_event': 'inverse-left high bit at index 0; acquired right label 0 belongs to H(point_0,point_1,point_0)',
            'joint_distribution': {str(out): str(p) for out, p in sorted(law.items())},
            'inverse_zero_marginal': str(inv0), 'H_contains_anchor_marginal': str(event),
            'product_of_marginals_at_0_True': str(inv0*event),
            'repeated_H_indices_0_1_1_anchor_agreement': '1',
            'distinct_H_indices_0_1_2_anchor_agreement': str(next(iter(independent_slot_agreement.values())))}


def partial_cache_native_fixture(E):
    # Explicit supported realization, common root Read schedule, no certification.
    assert (0, 0) in E
    prepared = [(0, 0), min(E), max(E), min(E)]
    term = positive.tree(tuple(b for y in prepared for b in positive.bits_of(y)))
    requested = [(0, h) for h in range(15)]+[(1, 0), (1, 0)]
    cache = [{} for _ in prepared]
    lengths = []
    for j, h in requested:
        word = (1,)*(15*j+h)+(0, 0)
        reply = positive.read(term, word)
        assert reply == (positive.BRANCH if positive.bits_of(prepared[j])[h] else positive.BETA)
        cache[j][h] = int(reply == positive.BRANCH)
        lengths.append(len(word))
    assert cache[1] == {0: 0} and not cache[2] and not cache[3]
    return {'N': len(prepared), 'Read_invocations': len(requested),
            'distinct_addresses': len(set(requested)), 'literal_operand_sum': sum(lengths),
            'literal_operand_max': max(lengths),
            'indexed_cache': [sorted(row.items()) for row in cache],
            'prepared_label_pairs': prepared, 'certified_E_assumed': False}


def full_injection_check(E0, inv0, E1, inv1):
    sources, raw = {}, {}
    scenarios = [
        ('unread_single', {}, (2,)),
        ('partial_single', {1: 0}, (1,)),
        ('partial_and_unread', {0: 0}, (0, 2)),
        ('two_partial_zero', {0: 0, 1: 0}, (0, 1)),
        ('two_partial_opposite', {0: 0, 1: 1}, (0, 1)),
        ('repeated_partial', {1: 0}, (1, 1)),
        ('repeated_unread', {}, (2, 2)),
        ('distinct_unread', {}, (2, 3)),
    ]
    counts = hypergeom(128, 64, 72)
    anchored_counts = {k+1: p for k, p in hypergeom(127, 63, 71).items()}
    for i, (E, inv) in enumerate(((E0, inv0), (E1, inv1))):
        J = {a: inv[a, t][0] for a, t in E}
        assert len(J) == 72 and all(J[J[a]] == a and J[a] != a for a in J)
        cases = {}
        for name, cache, targets in scenarios:
            z, law, distinct, assignments = left_inverse_projection(E, J, cache, targets)
            zeros = sum(b == 0 for b in cache.values())
            ones = len(cache)-zeros
            integral_z = sum(p*F(k, 72)**zeros*F(72-k, 72)**ones for k, p in counts.items())
            assert z == integral_z
            cases[name] = {'N': 4, 'left_high_bit_cache_by_original_index': cache,
                           'inverse_target_original_indices': targets,
                           'distinct_original_indices_marginalized': distinct,
                           'actual_left_assignments_enumerated': assignments,
                           'positive_evidence_mass': str(z),
                           'conditional_law': {str(out): str(p) for out, p in sorted(law.items())}}
            raw[i, name] = z, law
        z, law = raw[i, 'partial_single']
        correct = law[0,]
        direct = extension_mass(128, (64, 64), {0: 0, 1: 0})/extension_mass(128, (64, 64), {0: 0})
        normalized_Q = weighted_mean(counts, lambda k: F(k-1, 71), lambda k: F(k, 72))
        doubled = weighted_mean(counts, lambda k: F(k-1, 71), lambda k: F(k, 72)**2)
        assert correct == direct == normalized_Q == F(63, 127)
        assert doubled == F(287, 575) and doubled != correct
        anchored, anchored_z, anchored_cond_z = anchored_inverse_check(E, J)
        posterior = {k: p*F(k, 72)/anchored_cond_z for k, p in anchored_counts.items()}
        assert sum(posterior.values()) == 1
        anchored_Q = sum(p*F(k-1, 71) for k, p in posterior.items())
        anchored_mutant = weighted_mean(anchored_counts, lambda k: F(k-1, 71), lambda k: F(k, 72)**2)
        assert anchored == anchored_Q == F(287, 575)
        assert anchored_mutant == F(32879, 65475) and anchored_mutant != anchored
        repeat = raw[i, 'repeated_partial'][1]
        distinct_law = raw[i, 'two_partial_zero'][1]
        assert all(a == b for a, b in repeat)
        assert sum(p for (a, b), p in distinct_law.items() if a == b) < 1
        # A single shared injection also couples DISTINCT unread slots.
        unread = raw[i, 'distinct_unread'][1]
        assert unread[0, 0] != F(1, 4)
        assert unread[0, 0] == sum(p*F(k, 72)**2 for k, p in counts.items())
        for name in ('partial_and_unread', 'two_partial_zero', 'two_partial_opposite'):
            joint_law = raw[i, name][1]
            p = sum(v for (a, b), v in joint_law.items() if a == 0)
            q = sum(v for (a, b), v in joint_law.items() if b == 0)
            assert joint_law[0, 0] != p*q
        sources['G'+str(i)] = {
            'left_inverse_projections': cases,
            'partial_slot_weighting_discriminator': {
                'indexed_cache': 'N=4; index 1 left high bit 0; all other bits/indices unread',
                'target': 'inverse-left high bit at the SAME partial index',
                'actual_full_prior': str(correct), 'posterior_with_normalized_Q': str(normalized_Q),
                'incorrect_extra_support_or_likelihood': str(doubled), 'difference': str(doubled-correct),
                'occupancy_strata': len(counts), 'actual_occupancy_distribution': 'Hypergeom(128,64,72)'},
            'completed_anchor_weighting_discriminator': {
                'indexed_cache': 'N=4; index 0 complete (0,0); index 1 left high bit 0; indices 2,3 unread',
                'target': 'inverse-left high bit of acquired label 0 at index 0',
                'actual_full_prior': str(anchored), 'posterior_with_normalized_Q': str(anchored_Q),
                'incorrect_doubled_likelihood': str(anchored_mutant),
                'difference': str(anchored_mutant-anchored), 'positive_evidence_mass': str(anchored_z),
                'occupancy_strata': len(anchored_counts), 'actual_occupancy_distribution': '1+Hypergeom(127,63,71)',
                'actual_anchor_and_second_left_assignments_enumerated': 72**2},
            'attained_partial_cache_native_fixture': partial_cache_native_fixture(E)}
        if i == 0:
            sources['G0']['mixed_cache_joint_inverse_H'] = mixed_inverse_H_check(E, J)
    # A common source, not a fresh source per target. Here both endpoint-source
    # laws coincide, so the exact theta=1/3 mixture equals them; no all-theta
    # computation claim is needed.
    for name, _, _ in scenarios:
        z0, law0 = raw[0, name]
        z1, law1 = raw[1, name]
        zmix = F(2, 3)*z0+F(1, 3)*z1
        mixed = {out: (F(2, 3)*z0*law0[out]+F(1, 3)*z1*law1[out])/zmix for out in law0}
        assert law0 == law1 == mixed
    return {'prior': 'UJB.9001 full injections (128)_72 and (256)_243; one shared source and maps',
            'full_injections_enumerated': False, 'all_caches_enumerated': False,
            'marginalization': 'exact uniform full-injection extension counts; actual 72/243 endpoints and M=1944',
            'sources': sources, 'source_prior_values_checked': ['0', '1', '1/3'],
            'inverse_projection_scenarios_per_source': len(scenarios),
            'domain': 'named N=4 attained caches; inverse-bit and G0 quotient-fibre-event pushforwards only',
            'repeated_reads': 'same indexed restriction once; repeated targets use one endpoint assignment',
            'unread_slots': 'original prepared slots only; unused coordinates and injection values integrated exactly'}


def theory_identities():
    path = 'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_POSITIVE_RELATION_CERTIFICATES.md'
    versions = {}
    contents = []
    for rev in (THEORY_BASE, THEORY_PRE_REPAIR):
        raw = subprocess.check_output(['git', '-C', str(ROOT.parent/'trureturing'), 'show', rev+':'+path])
        versions[rev+':'+path] = {'sha256': sha256(raw).hexdigest(), 'lines': len(raw.splitlines())}
        contents.append(raw)
    assert contents[1].startswith(contents[0])
    return versions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', type=Path, help='Compare deterministic JSON without writing.')
    args = parser.parse_args()
    assert not (args.output and args.check)
    E0, inv0, g0 = group_source(0)
    E1, inv1, g1 = group_source(1)
    cases, geometries = {}, []
    for name, lm, rm in relabelings():
        assert len(set(lm)) == 72 and len(set(rm)) == 243
        geom, geo = geometry_check(E0, lm, rm)
        cases[name] = {'left_full_injection': lm, 'right_full_injection': rm, 'G0_geometry': geo,
                       'G0_conditional': conditional_check(geom, inv0, lm, rm),
                       'G0_H_conditional': geometry_law_check(geom),
                       'G1_domain_and_inverse': g1_check(E1, inv1, lm, rm)}
        geometries.append((name, lm, rm, geom))
    name, lm, rm, geom = geometries[2]
    native0 = native_acquisition(geom.E, transport(positive.g0_witness(E0), lm, rm))
    native1 = native_acquisition(transport(E1, lm, rm), transport(positive.g1_witness(E1), lm, rm))
    full_injection = full_injection_check(E0, inv0, E1, inv1)
    result = {'schema': 'rro-intrinsic-recovery-2',
              'immutable_theory_sources': theory_identities(),
              'full_injection_partial_cache': full_injection, 'theory_revision': THEORY_BASE,
              'experiment_base': EXPERIMENT_BASE, 'supplier_sha256': BLOBS,
              'arithmetic': 'integers, F3 and fractions.Fraction; MSB-first 7/8-bit labels; zero-based code indices',
              'group_sources': {'G0': g0, 'G1': g1}, 'relabelings': cases,
              'native_representative_relabeling': name, 'native_G0': native0, 'native_G1': native1,
              'G0_reference_memory': ambiguity_check(E0, lm, rm, geom),
              'coverage': {'full_injections_enumerated': False,
                           'all_caches_enumerated': False,
                           'exact_covers_per_relabeling': 40,
                           'full_injection_cache_domain': 'N=4; declared partial/complete/unread and repeated-index projections; exact extension-count marginalization', 'full_relabeling_examples': 4,
                           'charts_per_relabeling': 5184, 'total_charts': 20736,
                           'total_chart_point_comparisons': 559872, 'total_H_triples_enumerated': 78732,
                           'all_source_elements_per_source': 4374, 'group_adapter_edges_per_source': 1944,
                           'conditional_masks_per_G0_relabeling': 1979,
                           'same_observation_G1_witnesses_per_relabeling': 432,
                           'local_G1_bijections_per_relabeling': 40320,
                           'native_scans': 'one G0 and one G1 at fixed N; no fresh observer samples'},
              'limitations': ['Finite author validation, not independent review or Lean certification.',
                              'No enumeration of all full injections, all masks, arbitrary controllers, or all native words.',
                              'Full-injection tests integrate the original prior exactly for the named projections; they do not test every cache or prove the general law.',
                              'Theory identities are immutable pre-repair inputs; correspondence to terminal corrected definitions remains for the bounded alignment assignment.',
                              'Geometry consumes certified acquired E and its actual label associations; acquisition/certification suppliers are reused.',
                              'Affine chart invariance for all triples in all charts uses verified affine point maps and symbolic offset cancellation.',
                              'G1 inverse ambiguity requires positive G1 prior support; G0 recovery requires certified G0 E and an acquired actual input edge.',
                              'No hidden origin/injection, source mutation, inverse/group-operation port, fresh draw, actual serial history, quantum/random port, optimality, or resource improvement.']}
    rendered = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.check:
        assert args.check.read_text() == rendered, 'deterministic result mismatch'
    if args.output:
        args.output.write_text(rendered)
    print(json.dumps({'assertions': 'passed', 'relabelings': len(cases), 'charts': 20736,
                      'group_adapter_edges': 3888, 'result_sha256': sha256(rendered.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
