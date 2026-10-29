#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""One fixed, bounded source-component probe; computed evidence, not a proof."""
import datetime as dt
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import signal
import sys
import time
from fractions import Fraction as Q
from functools import lru_cache

PREFIX = 'Linusson-Restadh/2305.05275v1/Conj3.4/exact-refutation/v1/'
BODY_SHA = '57ac5cb0f76b61501769cb8e63abe5ace64f611de77299d3c841fe091592c4da'
CAPS = dict(n=7, raw_proposals=65536, fully_processed_eligible_pairs=256,
            components_min=3, components_max=12, selections_per_pair=4096,
            execution_seconds=3600, preparation_seconds=1200, overall_seconds=4800)
ROOT = Path(__file__).resolve().parent


def utcnow():
    return dt.datetime.now(dt.timezone.utc)


def atomic(path, obj):
    temp = path.with_name(path.name + '.tmp')
    with temp.open('w') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


class CapReached(Exception):
    pass


class Probe:
    def __init__(self, registration):
        self.registration = registration
        self.started = utcnow()
        self.start = time.monotonic()
        self.stop = self.start + CAPS['execution_seconds']
        self.overall = dt.datetime.fromisoformat(registration['operator_entry_at']) + dt.timedelta(seconds=4800)
        self.public_deadline = dt.datetime.fromisoformat(
            registration['strict_preregistration']['deadline'].replace('Z', '+00:00'))
        self.counts = {k: 0 for k in (
            'raw_proposals', 'equal_pairs', 'duplicate_pairs', 'unique_pairs',
            'nonchordal_pairs', 'nested_pairs', 'both_split_pairs',
            'cross_complete_neighborhood_pairs', 'components_below_3',
            'components_above_12', 'eligible_pairs_started',
            'fully_processed_eligible_pairs', 'component_selections',
            'compatible_selections', 'incompatible_cycle_selections',
            'incompatible_clique_selections', 'pairs_with_rhombus_witnesses',
            'pairs_without_new_complementary_pair', 'pair_lp_calls',
            'validation_lp_calls', 'lp_pivots', 'lp_feasible_results',
            'lp_no_feasible_result', 'graph_analyses', 'cycle_template_checks')}
        self.files = []
        self.current = None
        self.outcome = 'running'
        self.reason = None
        self.certificate = None
        self.error = None

    def tick(self):
        if time.monotonic() >= self.stop:
            raise CapReached('actual_execution_3600_seconds')
        if utcnow() >= self.overall:
            raise CapReached('overall_4800_seconds')
        if utcnow() >= self.public_deadline:
            raise CapReached('registered_absolute_deadline')

    def save(self, final=False):
        record = dict(outcome=self.outcome, stop_condition=self.reason,
                      caps=CAPS, execution_started_at=self.started.isoformat(),
                      execution_deadline_at=(self.started + dt.timedelta(seconds=3600)).isoformat(),
                      overall_deadline_at=self.overall.isoformat(),
                      preparation_elapsed_seconds=(self.started - dt.datetime.fromisoformat(
                          self.registration['operator_entry_at'])).total_seconds(),
                      last_raw_index=self.counts['raw_proposals'] - 1 if self.counts['raw_proposals'] else None,
                      actual_elapsed_seconds=time.monotonic() - self.start,
                      execution_ended_at=utcnow().isoformat() if final else None,
                      counters=self.counts, current_pair=self.current,
                      processed_pair_files=self.files, hit_certificate=self.certificate,
                      capability_error=self.error,
                      registration=dict(comment_id=self.registration['comment']['id'],
                          comment_url=self.registration['comment']['html_url'],
                          created_at=self.registration['comment']['created_at'],
                          updated_at=self.registration['comment']['updated_at'],
                          verified_at=self.registration['verified_at'],
                          body_sha256=self.registration['body_sha256']),
                      dependencies=dict(python=sys.version, executable=sys.executable,
                          platform=platform.platform(), third_party_packages=[]),
                      all_own_children_joined=True, child_processes_started=0,
                      formal_status='unresolved; no Lean or independent evaluation',
                      graph_cache= self.analyze.cache_info()._asdict()
                          if hasattr(self, 'analyze') else None)
        atomic(ROOT / 'results.json', record)
        return record

    def prepare(self):
        self.tick()
        self.edges = list(itertools.combinations(range(7), 2))
        self.edge_index = {e: i for i, e in enumerate(self.edges)}
        self.subsets = [s for s in range(128) if bin(s).count('1') >= 2]
        self.coord = {s: i for i, s in enumerate(self.subsets)}
        self.require = [sum(1 << i for i, (u, v) in enumerate(self.edges)
                            if s >> u & 1 and s >> v & 1) for s in range(128)]
        self.edge_coordinates = [(i, self.coord[(1 << u) | (1 << v)])
                                 for i, (u, v) in enumerate(self.edges)]
        self.cycles = []
        for size in range(4, 8):
            for vertices in itertools.combinations(range(7), size):
                for tail in itertools.permutations(vertices[1:]):
                    if tail[0] > tail[-1]:
                        continue
                    cycle = (vertices[0],) + tail
                    mask = 0
                    for u, v in zip(cycle, cycle[1:] + cycle[:1]):
                        mask |= 1 << self.edge_index[tuple(sorted((u, v)))]
                    subset = sum(1 << u for u in cycle)
                    self.cycles.append((mask, self.require[subset], cycle))
        assert len(self.edges) == 21 and len(self.subsets) == 120
        assert len(self.cycles) == 1137
        self.analyze = lru_cache(maxsize=16384)(self.graph_data)

    def vector(self, mask):
        return sum(1 << i for i, s in enumerate(self.subsets)
                   if mask & self.require[s] == self.require[s])

    def cycle(self, mask):
        for i, (edges, induced, vertices) in enumerate(self.cycles):
            if i % 256 == 0:
                self.tick()
            self.counts['cycle_template_checks'] += 1
            if mask & induced == edges:
                return vertices
        return None

    def graph_data(self, mask):
        self.tick()
        self.counts['graph_analyses'] += 1
        neighbors = [0] * 7
        for i, (u, v) in enumerate(self.edges):
            if mask >> i & 1:
                neighbors[u] |= 1 << v
                neighbors[v] |= 1 << u
        vector = self.vector(mask)
        cycle = self.cycle(mask)
        split = None
        if cycle is None:
            for clique in range(128):
                if (mask & self.require[clique] == self.require[clique]
                        and not mask & self.require[127 ^ clique]):
                    split = clique
                    break
        return vector, cycle, split, neighbors

    def components(self, cg, ch):
        different = [i for i in range(120) if (cg ^ ch) >> i & 1]
        remaining = set(different)
        difference_subsets = {self.subsets[i] for i in different}
        components = []
        while remaining:
            first = min(remaining)
            remaining.remove(first)
            todo = [first]
            comp = []
            while todo:
                a = todo.pop()
                comp.append(a)
                adjacent = [b for b in sorted(remaining)
                            if self.subsets[a] & self.subsets[b] in difference_subsets]
                for b in adjacent:
                    remaining.remove(b)
                    todo.append(b)
            components.append(sorted(comp))
        return components

    def selection_vector(self, base, comp_vectors, a):
        v = base
        for i, component in enumerate(comp_vectors):
            if a >> i & 1:
                v |= component
        return v

    def graph_from_vector(self, vector):
        return sum(1 << i for i, c in self.edge_coordinates if vector >> c & 1)

    def simplex_feasible(self, columns, rhs):
        """Rational phase I simplex with artificial basis and Bland's rule.

        This searches A*x=b, x>=0. A returned point is checked independently.
        Failure to find one is not exported as an edge/nonedge certificate.
        """
        self.tick()
        rows, original = len(rhs), len(columns)
        assert all(b >= 0 for b in rhs)
        table = [[Q(columns[j][i]) for j in range(original)]
                 + [Q(i == j) for j in range(rows)] + [Q(rhs[i])]
                 for i in range(rows)]
        basis = list(range(original, original + rows))
        costs = [Q(0)] * original + [Q(1)] * rows
        reduced = [costs[j] - sum(table[i][j] for i in range(rows))
                   for j in range(original + rows)]
        while True:
            self.tick()
            entering = next((j for j, c in enumerate(reduced) if c < 0), None)
            if entering is None:
                break
            possible = [(table[i][-1] / table[i][entering], basis[i], i)
                        for i in range(rows) if table[i][entering] > 0]
            if not possible:
                raise ArithmeticError('phase I unexpectedly unbounded')
            _, _, leaving = min(possible)
            pivot = table[leaving][entering]
            table[leaving] = [v / pivot for v in table[leaving]]
            for i in range(rows):
                self.tick()
                if i != leaving and table[i][entering]:
                    factor = table[i][entering]
                    table[i] = [a - factor * b
                                for a, b in zip(table[i], table[leaving])]
            factor = reduced[entering]
            reduced = [a - factor * b for a, b in zip(reduced, table[leaving][:-1])]
            basis[leaving] = entering
            self.counts['lp_pivots'] += 1
        point = [Q(0)] * (original + rows)
        for i, variable in enumerate(basis):
            point[variable] = table[i][-1]
        assert all(v >= 0 for v in point)
        if sum(point[original:]) != 0:
            return None
        point = point[:original]
        assert all(sum(columns[j][i] * point[j] for j in range(original)) == rhs[i]
                   for i in range(rows))
        return point

    def validate_solver(self):
        cases = [([[1, 0], [0, 1]], [1, 1], True),
                 ([[1, 1]], [1, 0], False),
                 ([[1, -1], [1, 1]], [1, 0], True)]
        for columns, rhs, feasible in cases:
            self.counts['validation_lp_calls'] += 1
            assert (self.simplex_feasible(columns, rhs) is not None) == feasible

    def certificate_for(self, g, h, components, a_g, compatible):
        other = [(a, d) for a, d in sorted(compatible.items()) if d not in (g, h)]
        self.counts['pair_lp_calls'] += 1
        k = len(components)
        g_bits = [(a_g >> i) & 1 for i in range(k)]
        h_bits = [1 - v for v in g_bits]
        columns = [[1] + [(a >> i) & 1 for i in range(k)] + [0] for a, d in other]
        columns += [[0] + [-(g_bits[i] - h_bits[i]) for i in range(k)] + [1],
                    [0] + [0] * k + [1]]
        point = self.simplex_feasible(columns, [1] + h_bits + [1])
        if point is None:
            self.counts['lp_no_feasible_result'] += 1
            return None
        self.counts['lp_feasible_results'] += 1
        t = point[-2]
        if not 0 < t < 1:
            raise ArithmeticError('endpoint-excluded binary hull produced noninterior point')
        terms = [(d, point[i]) for i, (a, d) in enumerate(other) if point[i]]
        assert terms and all(d not in (g, h) and w >= 0 for d, w in terms)
        assert sum(w for d, w in terms) == 1
        cg, ch = self.vector(g), self.vector(h)
        for d, weight in terms:
            assert self.cycle(d) is None
        for i in range(120):
            self.tick()
            assert t * ((cg >> i) & 1) + (1 - t) * ((ch >> i) & 1) == sum(
                weight * ((self.vector(d) >> i) & 1) for d, weight in terms)
        return dict(n=7, G_edge_mask=g, H_edge_mask=h, t=str(t),
                    terms=[dict(edge_mask=d, weight=str(w)) for d, w in terms],
                    exact_checks=dict(all_120_coordinates=True, normalized=True,
                        nonnegative=True, strict_interior=True, chordal=True,
                        endpoints_excluded=True), status='computed exact candidate; unformalized')

    def process_pair(self, j, g, h, cg, ch, components):
        k = len(components)
        self.counts['eligible_pairs_started'] += 1
        self.current = dict(raw_index=j, G_edge_mask=g, H_edge_mask=h,
                            components=k, selections_completed=0)
        self.save()
        comp_vectors = [sum(1 << i for i in comp) for comp in components]
        a_g = sum(1 << a for a, comp in enumerate(components) if cg >> comp[0] & 1)
        all_components = (1 << k) - 1
        a_h = a_g ^ all_components
        for comp in components:
            assert len({(cg >> i) & 1 for i in comp}) == 1
            assert len({(ch >> i) & 1 for i in comp}) == 1
        compatible, incompatible = {}, {}
        base = cg & ch
        assert self.selection_vector(base, comp_vectors, a_g) == cg
        assert self.selection_vector(base, comp_vectors, a_h) == ch
        for a in range(1 << k):
            self.tick()
            expected = self.selection_vector(base, comp_vectors, a)
            mask = self.graph_from_vector(expected)
            actual, cycle, split, neighbors = self.analyze(mask)
            mismatch = actual ^ expected
            if cycle is not None:
                incompatible[a] = dict(kind='chordless_cycle', vertices=list(cycle))
                self.counts['incompatible_cycle_selections'] += 1
            elif mismatch:
                first = (mismatch & -mismatch).bit_length() - 1
                assert bin(self.subsets[first]).count('1') >= 3
                incompatible[a] = dict(kind='clique_mismatch', subset_mask=self.subsets[first],
                                      expected=(expected >> first) & 1, actual=(actual >> first) & 1)
                self.counts['incompatible_clique_selections'] += 1
            else:
                compatible[a] = mask
                self.counts['compatible_selections'] += 1
            self.counts['component_selections'] += 1
            self.current['selections_completed'] += 1
        assert compatible[a_g] == g and compatible[a_h] == h
        assert len(compatible) + len(incompatible) == 1 << k <= 4096
        witness = next((a for a in sorted(compatible)
                        if a not in (a_g, a_h) and a ^ all_components in compatible), None)
        record = dict(raw_index=j, G_edge_mask=g, H_edge_mask=h,
                      G_clique_vector_hex=hex(cg), H_clique_vector_hex=hex(ch),
                      components_as_subset_masks=[[self.subsets[i] for i in comp] for comp in components],
                      G_selection=a_g, H_selection=a_h, selections_enumerated=1 << k,
                      compatible_family=[dict(selection=a, edge_mask=d)
                                         for a, d in sorted(compatible.items())],
                      incompatible_selections=incompatible, rhombus_witness=None,
                      exact_nonedge_candidate=None, lp_result=None)
        pair_path = ROOT / 'pairs' / ('pair-%05d.json' % j)
        atomic(pair_path, record)
        self.save()
        if witness is not None:
            d, e = compatible[witness], compatible[witness ^ all_components]
            assert d not in (g, h) and e not in (g, h)
            cd, ce = self.vector(d), self.vector(e)
            assert all(((cg >> i) & 1) + ((ch >> i) & 1)
                       == ((cd >> i) & 1) + ((ce >> i) & 1) for i in range(120))
            record['rhombus_witness'] = dict(selection=witness, G_prime_edge_mask=d, H_prime_edge_mask=e)
            self.counts['pairs_with_rhombus_witnesses'] += 1
        else:
            self.counts['pairs_without_new_complementary_pair'] += 1
            certificate = self.certificate_for(g, h, components, a_g, compatible)
            record['lp_result'] = ('exact_feasible_point_checked' if certificate else
                                   'no_feasible_point_found; no edge certificate claimed')
            if certificate:
                exclusions = []
                for a in range(1 << k):
                    b = a ^ all_components
                    if a in (a_g, a_h) or a > b:
                        continue
                    failed = a if a in incompatible else b
                    assert failed in incompatible
                    exclusions.append(dict(selection=a, complement=b, excluded_side=failed,
                                           obstruction=incompatible[failed]))
                certificate['complementary_exclusions'] = exclusions
                certificate['pair_data_path'] = str(pair_path.relative_to(ROOT))
                certificate['bridge_obligations'] = [
                    'native Fin 7 graphs and full source imsets',
                    '0/1 vertex correspondence in actual real convex hull',
                    'finite injective-cycle enumeration completeness',
                    'source Delta and complementary-selection witness completeness',
                    'strict interior convex equality implies ordinary nonface',
                    'exact rational arithmetic and all exclusion witnesses in Lean']
                record['exact_nonedge_candidate'] = certificate
                self.certificate = certificate
        atomic(pair_path, record)
        self.files.append(dict(path=str(pair_path.relative_to(ROOT)),
                               sha256=hashlib.sha256(pair_path.read_bytes()).hexdigest()))
        self.counts['fully_processed_eligible_pairs'] += 1
        self.current = None
        self.save()
        return self.certificate is not None

    def execute(self):
        self.prepare()
        self.validate_solver()
        seen = set()
        for j in range(65536):
            self.tick()
            self.counts['raw_proposals'] += 1
            digest = hashlib.sha256((PREFIX + str(j)).encode('utf-8')).digest()
            bits = int.from_bytes(digest[:6], 'big') % (1 << 42)
            g, h = sorted((bits & ((1 << 21) - 1), bits >> 21))
            if g == h:
                self.counts['equal_pairs'] += 1
                continue
            if (g, h) in seen:
                self.counts['duplicate_pairs'] += 1
                continue
            seen.add((g, h))
            self.counts['unique_pairs'] += 1
            cg, cycle_g, split_g, ng = self.analyze(g)
            ch, cycle_h, split_h, nh = self.analyze(h)
            if cycle_g is not None or cycle_h is not None:
                self.counts['nonchordal_pairs'] += 1
                continue
            if g & h in (g, h):
                self.counts['nested_pairs'] += 1
                continue
            if split_g is not None and split_h is not None:
                self.counts['both_split_pairs'] += 1
                continue
            if any(ng[x] != nh[x] and all(mask & self.require[nb] == self.require[nb]
                    for mask in (g, h) for nb in (ng[x], nh[x])) for x in range(7)):
                self.counts['cross_complete_neighborhood_pairs'] += 1
                continue
            components = self.components(cg, ch)
            if len(components) < 3:
                self.counts['components_below_3'] += 1
                continue
            if len(components) > 12:
                self.counts['components_above_12'] += 1
                continue
            if self.process_pair(j, g, h, cg, ch, components):
                return 'computed_exact_candidate'
            if self.counts['fully_processed_eligible_pairs'] == 256:
                return 'fully_processed_eligible_pairs_256'
        return 'raw_proposals_65536'


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: probe.py REGISTRATION_JSON')
    registration = json.loads(Path(sys.argv[1]).read_text())
    comment = registration['comment']
    assert registration['gate_passed'] and registration['body_sha256'] == BODY_SHA
    assert hashlib.sha256(comment['body'].encode()).hexdigest() == BODY_SHA
    assert registration['issue']['number'] == 23531
    assert comment['created_at'] == comment['updated_at']
    assert dt.datetime.fromisoformat(comment['created_at'].replace('Z', '+00:00')) < utcnow()
    assert utcnow() < dt.datetime.fromisoformat(
        registration['strict_preregistration']['deadline'].replace('Z', '+00:00'))
    assert (utcnow() - dt.datetime.fromisoformat(registration['operator_entry_at'])).total_seconds() < 1200
    assert not (ROOT / 'execution-start.json').exists(), 'single attempt already started; no resume'
    probe = Probe(registration)
    with (ROOT / 'execution-start.json').open('x') as f:
        json.dump(dict(started_at=probe.started.isoformat(),
                       execution_deadline_at=(probe.started + dt.timedelta(seconds=3600)).isoformat(),
                       overall_deadline_at=probe.overall.isoformat(), caps=CAPS,
                       registration_comment_id=comment['id']), f, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    (ROOT / 'pairs').mkdir(exist_ok=False)
    probe.save()
    def alarm(signum, frame):
        raise CapReached('actual_execution_3600_seconds')
    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, max(0.001, probe.stop - time.monotonic()))
    exit_code = 0
    try:
        probe.reason = probe.execute()
        probe.outcome = ('computed_exact_candidate' if probe.certificate else 'unresolved/defer')
    except CapReached as e:
        probe.reason = str(e)
        probe.outcome = 'unresolved/defer'
    except Exception as e:
        probe.outcome = 'unresolved/defer'
        probe.reason = 'capability_error_no_retry'
        probe.error = dict(type=type(e).__name__, message=str(e))
        exit_code = 2
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        result = probe.save(final=True)
    print(json.dumps(dict(outcome=result['outcome'], stop_condition=result['stop_condition'],
                         actual_elapsed_seconds=result['actual_elapsed_seconds'],
                         counters=result['counters'], capability_error=result['capability_error'])))
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
