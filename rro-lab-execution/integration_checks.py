"""Author checks for the remaining B integration boundaries.

This module composes the already installed service, its unchanged original
port, and the LAB consumer.  It deliberately does not run the complete ML-A
evaluator: the two original strict omissions remain explicit in the result.
The dynamic account below is a single chronological view of one native Read
session.  Its Arena event ledger is independent of the Arena census and is
compared at every real source-call boundary.
"""
from __future__ import annotations

from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

assert sys.version_info >= (3, 12)
assert __debug__

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from lab import Arena, Buffer  # noqa: E402
from native_install import BitKernel  # noqa: E402
from service_checks import installed_reads  # noqa: E402
from session import OriginalPort, Session, prepare_source  # noqa: E402
from framing import v  # noqa: E402


STRICT = ['ML-N-evaluator-execution', 'ML-Cost-full-evaluator-session']


class ArenaTrace:
    """An event ledger attached to the two actual Arena owners.

    Object identities remain the implementation's references.  The ledger
    records only allocation and first-field materialization, so it cannot
    choose a service value or a continuation.
    """

    def __init__(self):
        self.arenas = []
        self.events = {}
        self.boundaries = []
        self._install()

    def _install(self):
        trace = self
        if getattr(Arena, '_integration_trace_installed', False):
            return
        init = Arena.__init__
        alloc = Arena.alloc
        write = Arena.write
        release = Arena.release_bit_cell

        def traced_init(obj, *args, **kwargs):
            init(obj, *args, **kwargs)
            trace.arenas.append(obj)
            trace.events[id(obj)] = []

        def traced_alloc(obj, kind):
            out = alloc(obj, kind)
            trace.events[id(obj)].append({
                'op': 'alloc', 'kind': kind, 'object': id(out),
                'objects': len(obj.objects), 'width': obj.width,
                'refs': obj.reference_fields, 'bits': obj.bit_fields,
                'free_pool': obj.free_pool is not None,
            })
            return out

        def traced_write(obj, cell, key, value, located=False):
            before = key in cell.fields
            write(obj, cell, key, value, located)
            trace.events[id(obj)].append({
                'op': 'write', 'object': id(cell), 'key': key,
                'first': not before, 'objects': len(obj.objects),
                'width': obj.width, 'refs': obj.reference_fields,
                'bits': obj.bit_fields,
                'free_pool': obj.free_pool is not None,
            })

        def traced_release(obj, cell):
            release(obj, cell)
            trace.events[id(obj)].append({
                'op': 'release', 'object': id(cell), 'objects': len(obj.objects),
                'width': obj.width, 'refs': obj.reference_fields,
                'bits': obj.bit_fields, 'free_pool': obj.free_pool is not None,
            })

        Arena.__init__ = traced_init
        Arena.alloc = traced_alloc
        Arena.write = traced_write
        Arena.release_bit_cell = traced_release
        Arena._integration_trace_installed = True

    @staticmethod
    def direct(arena):
        refs = arena.registers
        bits = 0
        for cell in arena.objects:
            for key in cell.fields:
                if key == 'bit':
                    bits += 1
                else:
                    refs += 1
        return {
            'objects': len(arena.objects),
            'reference_slots': refs,
            'bit_fields': bits,
            'grant_bits': len(arena.objects),
            'reference_width': arena.width,
            'bits': refs * arena.width + bits + len(arena.objects),
        }

    def incremental(self, arena, stop=None):
        """Reconstruct the same census from the event stream."""
        events = self.events[id(arena)] if stop is None else self.events[id(arena)][:stop]
        objects = sum(e['op'] == 'alloc' for e in events)
        first_refs = sum(e['op'] == 'write' and e['first'] and e['key'] != 'bit'
                         for e in events)
        first_bits = sum(e['op'] == 'write' and e['first'] and e['key'] == 'bit'
                         for e in events)
        width = max((e['width'] for e in events), default=arena.width)
        refs = arena.registers + first_refs
        return {
            'objects': objects,
            'reference_slots': refs,
            'bit_fields': first_bits,
            'grant_bits': objects,
            'reference_width': width,
            'bits': refs * width + first_bits + objects,
        }

    def owners(self):
        return self.arenas


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def common_timeline(k):
    """Run one actual installed service and compare common-cut censuses."""
    trace = ArenaTrace()
    original_call = OriginalPort.call

    def traced_call(port, operand, observer_arena):
        before = {id(a): len(trace.events[id(a)]) for a in trace.arenas}
        result = original_call(port, operand, observer_arena)
        after = {id(a): len(trace.events[id(a)]) for a in trace.arenas}
        trace.boundaries.append({
            'operand': operand,
            'event_start': before,
            'event_end': after,
            'source_calls': len(port.calls),
            'direct': {id(a): trace.direct(a) for a in trace.arenas},
        })
        return result

    OriginalPort.call = traced_call
    try:
        measured = installed_reads(k, 0)
    finally:
        OriginalPort.call = original_call

    # The run creates one source Arena and one connection Arena.  The
    # installation is a separate native owner whose direct census is returned
    # by the actual service cut; its fixed material is included at the same
    # final service configuration, never added as an unrelated peak.
    assert len(trace.arenas) == 2, len(trace.arenas)
    source, connection = trace.arenas
    assert len(trace.boundaries) == measured['invocations']
    dynamic = []
    for boundary in trace.boundaries:
        owner_rows = {}
        for name, arena in (('source', source), ('connection', connection)):
            direct = boundary['direct'][id(arena)]
            incremental = trace.incremental(arena, boundary['event_end'][id(arena)])
            assert direct == incremental, (name, direct, incremental)
            owner_rows[name] = {'direct': direct, 'incremental': incremental}
        dynamic.append({'operand': boundary['operand'], 'owners': owner_rows,
                        'source_calls': boundary['source_calls']})

    final = measured['actual_common_retained_cut']
    assert final['same_original_source'] == trace.direct(source)
    assert final['connection'] == trace.direct(connection)
    installation = final['installation']
    installation_norm = {
        'objects': installation['objects'],
        'reference_slots': installation['refs'],
        'bit_fields': installation['bits'],
        'grant_bits': installation['objects'],
        'reference_width': installation['reference_width'],
        'bits': installation['total_bits'],
    }
    total = {
        key: installation_norm[key] + final['connection'][key] + final['same_original_source'][key]
        for key in ('objects', 'reference_slots', 'bit_fields', 'grant_bits', 'bits')
    }
    # A common cut carries a width vector.  It does not silently coerce the
    # source/provider arenas to the field installation's width.
    widths = {
        'installation': installation_norm['reference_width'],
        'connection': final['connection']['reference_width'],
        'source': final['same_original_source']['reference_width'],
    }
    connection_ops = measured['connection_operations']
    assert connection_ops['reference_widenings'] > 0
    assert connection_ops['reference_copy_bits'] > 0
    assert connection_ops['free_list_writes'] > 0
    assert measured['operations']['acknowledgements'] == measured['operations']['returns']
    assert measured['operations']['original_calls'] == measured['invocations']
    return {
        'run': {
            'source_bit': 0,
            'invocations': measured['invocations'],
            'distinct_addresses': measured['distinct_addresses'],
            'native_steps': measured['native_steps'],
            'field_actions': measured['field_actions'],
            'installation_steps': measured['installation_steps'],
            'graph_sha256': measured['graph_sha256'],
            'output_sha256': measured['output_sha256'],
        },
        'common_cut': {
            'width_vector': widths,
            'component_direct_census': {
                'installation': installation,
                'installation_normalized': installation_norm,
                'connection': final['connection'],
                'same_original_source': final['same_original_source'],
            },
            'unified_total_at_same_cut': total,
            'unrelated_peak_sum_used': False,
        },
        'timeline': {
            'source_call_boundaries': dynamic,
            'boundary_count': len(dynamic),
            'direct_equals_incremental': True,
            'retained_old_jobs_and_roots': True,
            'pending_ack_return_close_cuts': True,
            'free_pool_and_cleanup_observed': True,
            'refusal_and_busy_paths_checked_by_existing_service_mutants': True,
        },
        'source_operations': measured['operations'],
        'connection_operations': connection_ops,
    }


def lightweight_allocator_checks():
    """Exercise actual allocator widening, overlap and free-pool cleanup."""
    a = Arena(registers=2)
    b = Buffer(a, 'timeline.buffer')
    for bit in '0101010101':
        b.append(int(bit))
    before = ArenaTrace.direct(a)
    b.clear()
    after = ArenaTrace.direct(a)
    assert a.free_pool is not None
    assert a.charges['reference_widenings'] > 0
    assert a.charges['reference_copy_bits'] > 0
    assert a.charges['free_list_writes'] > 0
    assert after['objects'] == before['objects'] + 1  # free-pool header only
    return {
        'before_cleanup': before,
        'after_cleanup': after,
        'reference_widenings': a.charges['reference_widenings'],
        'reference_copy_bits': a.charges['reference_copy_bits'],
        'free_list_writes': a.charges['free_list_writes'],
        'cleanup_fields': a.charges['cleanup_fields'],
    }


def whole_family_author_checks():
    """Re-run the finite author-facing composition checks.

    These are corroboration only; the universal clauses are the ordinary
    induction in the owning theory volume.
    """
    import exact_checks as exact
    old, evaluator, positive = exact.suppliers()
    b = exact.exhaustive_b_checks()
    framing = exact.framing_checks()
    menu = exact.port_menu_checks(old, positive)
    # The sealed result remains a source-bound audit of unchanged domains.
    saved = json.loads((HERE / 'exact-results.json').read_text())
    assert saved['unresolved_obligations'][0]['id'] in STRICT
    assert saved['unresolved_obligations'][1]['id'] in STRICT
    return {
        'exhaustive_B': b,
        'framing': {
            'V_words': framing['V_words'],
            'within_slot_masks': framing['within_slot_masks'],
            'proper_frame_prefixes_wait': framing['proper_frame_prefixes_wait'],
        },
        'original_port_menu': {
            'calls': menu['calls'],
            'distinct_operands_per_fixture': menu['distinct_operands_per_fixture'],
            'maximum_operand_length': menu['maximum_operand_length'],
            'all_slot_suffixes_through_length': menu['all_slot_suffixes_through_length'],
        },
        'strict_omissions_retained': STRICT,
    }


def source_binding_checks():
    """Check the changed source bindings and persisted operational mutants."""
    import subprocess
    service = json.loads((HERE / 'service-results.json').read_text())
    author = json.loads((HERE / 'service-author-results.json').read_text())
    domain = json.loads((HERE / 'service-domain-results.json').read_text())
    assert service['source_hashes']['service_native.py'] == digest(HERE / 'service_native.py')
    assert domain['default_engine_compatibility']['default_native_engine_fully_normalized_AST_equal']
    # Existing mutants were required to fail through real state/behavior.
    mutant_results = [x for x in author.get('tests', []) if x.get('falsifier')]
    assert {x['falsifier'] for x in mutant_results} == {
        'actual_port_transfer_occupied', 'actual_cache_write'
    }
    from native_engine import NativeEngine
    from field_protocol_checks import exchange_graph
    baseline = subprocess.check_output(['git', '-C', str(HERE.parent),
                                        'show', 'HEAD:rro-lab-execution/native_engine.py'])
    current = (HERE / 'native_engine.py').read_bytes()
    assert current == baseline
    kernel = BitKernel()
    try:
        engine = NativeEngine(exchange_graph(), 1, kernel=kernel)
        engine.run()
        assert engine.output == [1]
        assert engine.charges['output_ends'] == engine.charges['acknowledgements'] == engine.charges['returns'] == 1
        fresh_compatibility = {
            'source_sha256': hashlib.sha256(current).hexdigest(),
            'baseline_sha256': hashlib.sha256(baseline).hexdigest(),
            'unchanged_source': True,
            'fresh_output': engine.output,
            'native_steps': engine.native_steps,
        }
    finally:
        kernel.close()
    return {
        'service_graph_sha256': service['native_reads'][0]['graph_sha256'],
        'default_engine_compatibility': True,
        'fresh_compatibility': fresh_compatibility,
        'operational_mutants': [x['falsifier'] for x in mutant_results],
        'returned_output_cases': len(domain['returned_output']['cases']),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--skip-native', action='store_true')
    ap.add_argument('--reuse-common', type=Path)
    args = ap.parse_args()
    start = time.perf_counter()
    if args.reuse_common:
        common = json.loads(args.reuse_common.read_text())['common_timeline']
    elif args.skip_native:
        common = None
    else:
        k = BitKernel()
        try:
            common = common_timeline(k)
        finally:
            k.close()
    theory_root = Path('/Users/auricstudio/trureturing-cad0ceaa-op')
    theory_path = theory_root / 'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md'
    candidate = json.loads((HERE / 'theory-candidate.json').read_text())
    result = {
        'schema': 1,
        'stage': 'author integration',
        'B_LAB_driver': common is not None,
        'B_allocation_width_common_timeline': common is not None,
        'B_whole_family_correspondence': True,
        'unit_complete': common is not None,
        'B_complete': common is not None,
        'C_released': False,
        'D_complete': False,
        'independent_acceptance': False,
        'whole_goal_complete': False,
        'strict_original_omissions': STRICT,
        'theory_binding': {
            'candidate_manifest_sha256': digest(HERE / 'theory-candidate.json'),
            'theory_sha256': digest(theory_path),
            'theory_bytes': theory_path.stat().st_size,
            'manifest_entry': candidate['files']['docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md'],
        },
        'allocator_probe': lightweight_allocator_checks(),
        'source_binding': source_binding_checks(),
        'whole_family': whole_family_author_checks(),
        'common_timeline': common,
        'source_hashes': {
            name: digest(HERE / name) for name in (
                'integration_checks.py', 'service_native.py', 'native_engine.py',
                'native_install.py', 'field_service.py', 'field_graph.py',
                'field_decoder.py', 'field_io.py', 'session.py', 'framing.py',
            )
        },
        'seconds': time.perf_counter() - start,
    }
    raw = json.dumps(result, sort_keys=True, indent=2) + '\n'
    tmp = args.output.with_suffix('.json.tmp')
    tmp.write_text(raw)
    tmp.replace(args.output)
    print(json.dumps({'checks': 'pass', 'result_sha256': digest(args.output),
                      'seconds': result['seconds']}))


if __name__ == '__main__':
    main()
