"""Acquired-input attempts on the unchanged complete production evaluator.

A finite supply refusal is an incomplete run. No fixture or reference answer
is an execution input. Inspection below cannot select a machine successor.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction
import argparse
import gc
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
import time

assert sys.version_info >= (3, 12) and __debug__
HERE = Path(__file__).resolve().parent
from constructor_pack import unpack
from framing import (Reader, REPLIES, parse_arithmetic, parse_reply, reply,
                     u, v, rational, words, cylinder_word)
from field_graph import OPS
from field_service import service_graph
from field_io import observe_buffer
from native_install import BitKernel
from paid_code import CodeStore
from service_native import OriginalConnection, ServiceEngine
from service_checks import receipts, cache
from session import OriginalPort, prepare_source
from lab import Arena
from pins import source_manifest

CASES = ('positive-full', 'positive-partial', 'null-nine-distinct-left')
STRICT = ['ML-N-evaluator-execution', 'ML-Cost-full-evaluator-session']


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def file_identity(path):
    raw = Path(path).read_bytes()
    return {'bytes': len(raw), 'sha256': sha(raw)}


def atomic(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, sort_keys=True, indent=2) + '\n')
    tmp.replace(path)


class Chronology:
    """Independent material deltas and common-cut peaks, without event logs.

    Arena.sample is the original material-change cut, including real width
    overlap. Installer callbacks are its actual post-primitive cuts. Since
    runtime native transitions do not allocate, they cannot change its census.
    """
    def __init__(self):
        self.arenas = []
        self.delta = {}
        self.engine = None
        self.peak = None
        self.cut = 0
        self.checks = 0
        self.points = {}
        self.original = {n: getattr(Arena, n) for n in ('__init__', 'alloc', 'write', 'sample')}
        trace = self
        def init(a, *args, **kwargs):
            trace.original['__init__'](a, *args, **kwargs)
            trace.arenas.append(a)
            trace.delta[id(a)] = Counter(objects=0, refs=a.registers, bits=0)
        def alloc(a, kind):
            out = trace.original['alloc'](a, kind)
            trace.delta[id(a)]['objects'] += 1
            return out
        def write(a, obj, key, value, located=False):
            first = key not in obj.fields
            out = trace.original['write'](a, obj, key, value, located)
            if first:
                trace.delta[id(a)]['bits' if key == 'bit' else 'refs'] += 1
            return out
        def sample(a, overlap=0):
            trace.original['sample'](a, overlap)
            trace.sample('arena-material', a, overlap)
        Arena.__init__, Arena.alloc, Arena.write, Arena.sample = init, alloc, write, sample

    def close(self):
        for name, original in self.original.items():
            setattr(Arena, name, original)

    def installation(self):
        e = self.engine
        if e is None or e.installation is None:
            return None
        i = e.installation
        # Seed/arity formula is checked against actual traversal at named cuts.
        seed = self.seed if hasattr(self, 'seed') else i.seed_census()
        self.seed = seed
        c = {k: i.base[k] + i.accounted[k] + seed[k] for k in ('objects', 'refs', 'bits')}
        return {'objects': c['objects'], 'reference_slots': c['refs'],
                'bit_fields': c['bits'], 'grant_bits': c['objects'],
                'reference_width': i.width,
                'bits': c['bits'] + i.width*c['refs'] + c['objects']}

    def sample(self, phase, overlap_owner=None, overlap=0):
        self.cut += 1
        components = {('source' if at == 0 else 'connection'): a.census()
                      for at, a in enumerate(self.arenas)}
        installation = self.installation()
        if installation is not None:
            components['installation'] = installation
        total = {k: sum(c[k] for c in components.values())
                 for k in ('objects', 'reference_slots', 'bit_fields', 'grant_bits', 'bits')}
        total['bits'] += overlap
        if self.peak is None or total['bits'] > self.peak['total']['bits']:
            self.peak = {'cut': self.cut, 'phase': phase, 'components': components,
                         'overlap_bits': overlap,
                         'overlap_owner': ('source' if overlap_owner is self.arenas[0] else 'connection') if overlap else None,
                         'total': total}
        return {'cut': self.cut, 'phase': phase, 'components': components, 'total': total}

    def checkpoint(self, installer):
        self.sample('installation:' + installer.phase)

    def point(self, name):
        for a in self.arenas:
            d = self.delta[id(a)]
            direct_refs = a.registers + sum(k != 'bit' for o in a.objects for k in o.fields)
            direct_bits = sum(k == 'bit' for o in a.objects for k in o.fields)
            assert (d['objects'], d['refs'], d['bits']) == (len(a.objects), direct_refs, direct_bits)
            self.checks += 1
        if self.engine is not None and self.engine.installation is not None:
            direct = self.engine.installation.census()
            inc = self.installation()
            assert (direct['objects'], direct['refs'], direct['bits'], direct['total_bits']) == (
                inc['objects'], inc['reference_slots'], inc['bit_fields'], inc['bits'])
            self.checks += 1
        p = self.sample(name)
        self.points[name] = p
        return p


class MeasuredEngine(ServiceEngine):
    """Inspection hooks only; the inherited primitive transition is unchanged."""
    def __init__(self, *args, **kwargs):
        self.field_operations = Counter()
        self.boundaries = Counter()
        super().__init__(*args, **kwargs)
        self.labels_by_node = {}
        for name, node in self.program.entries.items():
            self.labels_by_node.setdefault(node, []).append(name)

    def tick(self):
        super().tick()
        self.field_operations[self.last.op] += 1
        for label in self.labels_by_node.get(self.last, ()):
            self.boundaries[label] += 1


def cache_tuple(e, n):
    paired = cache(e)
    assert len(paired) == 15*n
    return tuple(tuple((h, paired[15*j+h][1]) for h in range(15) if paired[15*j+h][0])
                 for j in range(n))


def chain_length(cell):
    n = 0
    while cell is not None:
        n += 1
        cell = cell.refs[1]
    return n


def work_view(e):
    root = e.registers['root']
    if root is None:
        return None
    arm_values = []
    vertices, headers, route = set(), set(), set()
    todo = [root]
    while todo:
        p = todo.pop()
        assert p not in vertices
        vertices.add(p)
        todo.extend(c for c in p.refs[1:3] if c is not None)
        if p.refs[3] is not None:
            h = p.refs[3]
            headers.add(h)
            c = h.refs[1]
            while c is not None:
                assert c not in route
                route.add(c)
                c = c.refs[1]
    for key in (1, 2):
        p, values = root.refs[key], []
        while p is not None and p.refs[2] is not None:
            values.append(int(p.refs[2].refs[1] is not None))
            p = p.refs[1]
        arm_values.append(values)
    return {'source_vertices': len(vertices), 'routing_headers': len(headers),
            'routing_cells': len(route), 'objects': len(vertices)+len(headers)+len(route),
            'minus_bits': arm_values[0], 'plus_bits': arm_values[1],
            'growth_target_present': e.registers['target'] is not None,
            'pending_ancestor_present': e.registers['growth.parent'] is not None,
            'pending_append_cell_present': e.registers['growth.child'] is not None}


def required_loader(converted):
    m = len(converted)
    one_positions = [j for j, b in enumerate(converted) if b == '1']
    b = len(one_positions)
    d = 1+m*(m+1)//2+sum(j+2 for j in one_positions)
    return {'B_input_bits': m, 'marked_input_bits': b, 'P': 2+b, 'S': 1+m,
            'V': 5+2*m+b, 'H': 1+m, 'D_route': d,
            'work_objects_at_complete_load': 6+3*m+b+d,
            'derivation': 'ML.10 at n_minus=0,n_plus=m, including every loaded one; not an executed allocation or time estimate'}


def run_case(name, capacity, theory_root, output):
    start = time.perf_counter()
    input_path = HERE/'constructor'/(name+'.input.bits')
    requested = input_path.read_text().strip()
    n, theta, expected_cache, cylinders = parse_arithmetic(requested)
    bits = [dict(row).get(h, 0) for row in expected_cache for h in range(15)]
    trace = Chronology()
    kernel = None
    try:
        source = prepare_source(bits)
        source_state = (source.logical(), tuple(source.events), source.routing())
        preparation = Counter(source.a.charges)
        port = OriginalPort(source)
        connection = OriginalConnection(port)
        trace.point('source-prepared')
        with tempfile.TemporaryDirectory(prefix='rro-full-execution-') as td:
            packed = Path(td)/'W_U.block.gz'
            unpack(HERE/'constructor-code', 'W_U.block.gz', packed)
            raw = gzip.decompress(packed.read_bytes())
            interface = json.loads((HERE/'constructor/INTERFACE.json').read_text())
            assert len(raw) == interface['W_U']['uncompressed_image_bytes']
            assert sha(raw) == interface['W_U']['uncompressed_image_sha256']
            store = CodeStore(raw)
            p = service_graph(store.height, 361)
            kernel = BitKernel()
            e = MeasuredEngine(p, capacity, bindings={'image.root': store.root},
                               store=store, kernel=kernel, connection=connection)
            trace.engine = e
            e.seal(trace.checkpoint)
            trace.point('installed')
            print(json.dumps({'case': name, 'installed_steps': e.installation.steps}), flush=True)
            e.run()
            e.feed(u(n)+rational(theta))
            assert e.bits['initialized'] and not e.output
            acquisition = ['', 'LL', 'LRR', 'LLR', 'LL', 'LR', 'LLL', 'LRL']
            acquisition += ['R'*b+'LL' for b in range(15*n) if (b%15, bits[b]) in expected_cache[b//15]]
            acquisition += ['LL', 'R'*15+'LL', 'R'*45+'LL', 'L'*130, 'R'*(15*n+10)]
            read_outputs = []
            for index, operand in enumerate(acquisition):
                before = len(e.output)
                e.feed('01'+v(operand.translate(str.maketrans('LR', '01'))))
                value = ''.join(map(str, e.output[before:]))
                assert parse_reply(value) == (2, format(REPLIES.index(port.calls[-1][1]), '02b'))
                read_outputs.append(value)
                assert receipts(e) == [(a, b) for a, b, _, _ in port.calls]
                trace.point('read-return:'+str(index))
            assert cache_tuple(e, n) == expected_cache
            assert set(v for _, v, _, _ in port.calls) == set(REPLIES)
            calls_before = len(port.calls)
            native_before = e.native_steps
            output_before = len(e.output)
            if name == 'positive-full':
                # An actual caller supplies the acquired restriction verbatim;
                # the production service validates these received bytes.
                command = '00'+v(requested)
                mode = 'Evaluate'
            else:
                command = '101'+words([cylinder_word(c) for c in cylinders])
                mode = 'Predict-cylinder'
            # Supplying the last End separately permits measurement of the
            # actual autonomous continuation, without replacing its PC.
            e.feed(command, seal=False)
            e.bw('in.bit', 0)
            e.bw('in.end', 1)
            e.bw('in.ready', 1)
            e.resume()
            stages = {}
            obstruction = None
            try:
                while not e.waiting:
                    if e.pc is p.entries['converted.done'] and 'converted' not in stages:
                        converted = observe_buffer(e.registers['converted'])
                        actual_input = observe_buffer(e.registers['eval.input'])
                        assert actual_input == requested
                        # Independent conversion is inspection, never loaded.
                        from constructor_compiler import initial_binary
                        from copier import initial_word
                        from nwmg import encode
                        assert converted == encode(initial_word(initial_binary(requested, 361)))
                        stages['converted'] = {'input': actual_input,
                            'input_identity': {'bits': len(actual_input), 'sha256': sha(actual_input.encode())},
                            'B_input_bits': len(converted), 'B_input_sha256': sha(converted.encode()),
                            'field_actions': e.ticks, 'native_steps': e.native_steps,
                            'required_complete_loader': required_loader(converted),
                            'cut': trace.point('converted')}
                        print(json.dumps({'case': name, 'converted_bits': len(converted)}), flush=True)
                    e.tick()
            except RuntimeError as exc:
                obstruction = str(exc)
                assert e.pc is p.entries['exhausted'], ('unexpected refusal', str(exc))
                assert e.registers['pool'] is None
            actual_arithmetic_output = ''.join(map(str, e.output[output_before:]))
            completed = obstruction is None and e.waiting and not e.bits['eval.busy']
            if completed:
                from framing import parse_arithmetic_output
                kind, payload = parse_reply(actual_arithmetic_output)
                assert kind == 1
                parsed = parse_arithmetic_output(payload)
                expected_output = (HERE/'constructor'/(name+'.output.bits')).read_text().strip()
                assert payload == expected_output
                before = len(port.calls)
                e.feed('01'+v('00'))
                assert len(port.calls) == before+1
                continued = True
            else:
                parsed = None
                continued = False
                assert not actual_arithmetic_output and e.bits['eval.busy']
            assert len(port.calls) == calls_before + int(continued)
            assert source_state == (source.logical(), tuple(source.events), source.routing())
            assert e.native_steps == sum(e.native_events.values())
            assert e.native_steps == e.charges['native_controller_selections']
            assert e.charges['code_bit_reads'] == 21*(e.ticks + int(obstruction is not None))
            assert e.charges['code_reference_reads'] == 4*(e.ticks + int(obstruction is not None))
            assert e.charges['original_calls'] == connection.calls == len(port.calls)
            final = trace.point('terminal')
            pending = {k: chain_length(e.registers[k]) for k in ('input',)}
            occupied = sum(c.bit(7) for c in e.cells)
            assert occupied == capacity if obstruction else occupied <= capacity
            result = {
                'schema': 1, 'case': name, 'scope': 'author finite production attempt',
                'C_complete': completed, 'D_complete': False, 'independent_acceptance': False,
                'whole_goal_complete': False, 'strict_missing_objects': [] if completed else STRICT,
                'input_fixture': {'path': str(input_path.relative_to(HERE)), **file_identity(input_path)},
                'received_command': command, 'received_command_bits': len(command),
                'received_command_sha256': sha(command.encode()), 'mode': mode,
                'actual_original_source': {'N': n, 'z': bits, 'z_sha256': sha(bytes(bits)),
                    'same_bound_root': port.root is source.root, 'current_grant': port.grant,
                    'preparation_events': len(source.events), 'persistent': source.persistent(),
                    'preparation_operations': dict(preparation), 'source_unchanged': True},
                'source_calls': [dict(operand=w, reply=v, examined_vertices=x, edges=y) for w,v,x,y in port.calls],
                'receipt_occurrences': receipts(e), 'cache': cache(e),
                'physical_read_service_outputs': read_outputs,
                'physical_arithmetic_service_output': actual_arithmetic_output,
                'decoded_ML5': parsed, 'continued_original_Read_after_evaluation': continued,
                'predict_or_evaluate_source_calls': len(port.calls)-calls_before-int(continued),
                'counts': {'field_actions': e.ticks, 'native_steps': e.native_steps,
                    'installation_steps': e.installation.steps, 'code_graph_nodes': len(p.nodes),
                    'installed_native_nodes': sum(1 for s in e.installation.sources if s.block.width==21)-len(p.nodes),
                    'B_instruction_commits': e.boundaries['commit'],
                    'copier_completed_steps': 0 if not e.boundaries['commit'] else None,
                    'binary_completed_steps': 0 if not e.boundaries['commit'] else None,
                    'stack_completed_steps': 0 if not e.boundaries['commit'] else None,
                    'source_calls': len(port.calls), 'occupied_reserve_cells': occupied,
                    'reserve_capacity': capacity, 'field_operations': dict(e.field_operations),
                    'graph_boundaries': dict(e.boundaries)},
                'operations': {'installation': dict(e.installation.charges),
                    'native': dict(e.charges), 'native_primitives': dict(e.native_events),
                    'source_preparation': dict(preparation),
                    'source_provider': dict(source.a.charges-preparation),
                    'connection_and_provider_buffers': dict(connection.a.charges)},
                'resource_account': {'scope': 'same actual finite prefix; full ML.16 unfinished',
                    'complete_full_session_vector': completed,
                    'retained_common_cut': final, 'peak_common_cut': trace.peak,
                    'direct_incremental_equal_named_cuts': trace.checks,
                    'named_cuts': trace.points,
                    'unrelated_component_peak_sum_used': False,
                    'legacy_session_records_and_Print_emitted_consumed': False,
                    'runtime_installer_fixed_width': e.installation.width},
                'stages': stages,
                'obstruction': None if completed else {'kind': 'actual finite reserve exhausted',
                    'message': obstruction, 'interpreter_pc': 'exhausted',
                    'native_refusal_occurrences': e.charges['native_refusals'],
                    'eval_busy': e.bits['eval.busy'], 'pool_is_null': e.registers['pool'] is None,
                    'work': work_view(e), 'remaining_input_chain_bits': pending['input'],
                    'remaining_object': 'Complete original ML-L load, then fixed-image copier/B/LAB evaluator, actual returned-head decoding, Arithmetic Print/End/Ack/Return and subsequent same-port Read; full ML.16 of that completed chronology.'},
                'fixed_code': {'graph_sha256': p.fingerprint(), 'image': file_identity(packed),
                    'uncompressed_image_bytes': len(raw), 'uncompressed_image_sha256': sha(raw)},
                'theory_binding': {'root': str(theory_root), 'manifest': file_identity(HERE/'theory-candidate.json'),
                    'volume': file_identity(theory_root/'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md')},
                'source_hashes': {p.name: file_identity(p) for p in sorted(HERE.glob('*.py'))},
                'runtime': {'python': sys.version, 'assertions_enabled': __debug__,
                    'scope_descriptors_inherited': len(os.environ.get('WORKTREE_SCOPE_FDS','').split(','))},
                'host_peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                'seconds': time.perf_counter()-start,
            }
            atomic(output, result)
            print(json.dumps({'case': name, 'complete': completed, 'obstruction': obstruction,
                              'native_steps': e.native_steps, 'result': str(output)}), flush=True)
            return 0 if completed else 2
    finally:
        trace.close()
        if kernel is not None:
            kernel.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=CASES, required=True)
    ap.add_argument('--capacity', type=int, default=200000)
    ap.add_argument('--theory-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    # pins reads actual candidate bytes and all immutable original suppliers.
    source_manifest(args.theory_root, HERE/'theory-candidate.json')
    raise SystemExit(run_case(args.case, args.capacity, args.theory_root, args.output))


if __name__ == '__main__':
    main()
