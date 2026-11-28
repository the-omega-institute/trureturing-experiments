"""Author validation of incomplete original production attempts and falsifiers.

Exit zero validates the reported finite prefixes, never original C acceptance.
"""
from __future__ import annotations
import argparse
from collections import Counter
import copy
import hashlib
import json
import os
from pathlib import Path
import sys

assert sys.version_info >= (3, 12) and __debug__
ROOT = Path(__file__).resolve().parent
HERE = ROOT/'rro-lab-execution'
sys.path.insert(0, str(HERE))
from full_execution import atomic, file_identity, required_loader, STRICT
from framing import REPLIES, Reader, parse_arithmetic, parse_reply, u, v, rational
from field_graph import OPS, install_b
from field_io import output_graph
from field_service import service_graph
from service_native import OriginalConnection
from full_execution import MeasuredEngine
from native_engine import NativeEngine
from service_checks import cache, receipts
from session import OriginalPort, prepare_source
from native_install import BitKernel
from lab import Arena, Buffer
from pins import source_manifest


class MeasuredNativeEngine(NativeEngine):
    """Inspect generic component boundaries without adding a service endpoint."""
    def __init__(self, *args, **kwargs):
        self.boundaries = Counter()
        super().__init__(*args, **kwargs)
        self.labels_by_node = {}
        for name, node in self.program.entries.items():
            self.labels_by_node.setdefault(node, []).append(name)

    def tick(self):
        super().tick()
        for label in self.labels_by_node.get(self.last, ()):
            self.boundaries[label] += 1


def rejected(fn):
    try:
        fn()
    except (AssertionError, RuntimeError):
        return True
    raise AssertionError('falsifier was not detected')


def validate_native(ops, events, steps, actions, refused=False):
    assert sum(events.values()) == steps
    assert ops['native_controller_selections'] == steps
    assert ops['native_successor_reference_reads'] == steps-int(refused)
    assert ops['native_control_reference_writes'] == steps-int(refused)
    assert ops['field_commits'] == actions == events['commit']
    assert ops['code_bit_reads'] == events['codebit'] == 21*(actions+int(refused))
    assert ops['code_reference_reads'] == events['coderef'] == 4*(actions+int(refused))
    assert ops['bit_field_reads'] == events['codebit']+events['bfield']
    assert ops['reference_field_reads'] == events['coderef']+events['rfield']
    assert ops['bit_field_writes'] == events['vfield']
    assert ops['reference_field_writes'] == events['wfield']


def validate_attempt(d):
    assert not d['C_complete'] and not d['D_complete'] and not d['independent_acceptance']
    assert d['strict_missing_objects'] == STRICT
    count, ops = d['counts'], d['operations']
    validate_native(ops['native'], ops['native_primitives'], count['native_steps'], count['field_actions'], True)
    assert count['occupied_reserve_cells'] == count['reserve_capacity']
    assert count['B_instruction_commits'] == count['copier_completed_steps'] == count['binary_completed_steps'] == count['stack_completed_steps'] == 0
    assert d['obstruction']['pool_is_null'] and d['obstruction']['eval_busy']
    assert d['obstruction']['native_refusal_occurrences'] == 1
    assert not d['physical_arithmetic_service_output'] and d['decoded_ML5'] is None
    assert not d['continued_original_Read_after_evaluation']
    assert d['predict_or_evaluate_source_calls'] == 0
    actual = d['stages']['converted']['input']
    n, theta, data, cylinders = parse_arithmetic(actual)
    assert actual == (HERE/d['input_fixture']['path']).read_text().strip()
    paired = [(1, dict(row)[h]) if h in dict(row) else (0, 0) for row in data for h in range(15)]
    assert d['cache'] == [list(x) for x in paired]
    assert d['receipt_occurrences'] == [[c['operand'], c['reply']] for c in d['source_calls']]
    assert len(d['physical_read_service_outputs']) == count['source_calls'] == len(d['source_calls'])
    assert set(c['reply'] for c in d['source_calls']) == set(REPLIES)
    for c, frame in zip(d['source_calls'], d['physical_read_service_outputs']):
        assert parse_reply(frame) == (2, format(REPLIES.index(c['reply']), '02b'))
    calls = d['source_calls']
    co = ops['connection_and_provider_buffers']
    assert co['original_calls'] == count['source_calls'] == ops['native']['original_calls']
    assert co['port_logical_operand_bits'] == sum(len(c['operand']) for c in calls)
    assert co['port_reply_bits'] == 2*len(calls)
    assert co['port_returns'] == co['port_acknowledgements'] == len(calls)
    assert ops['source_provider']['field_reads'] == 2*sum(c['examined_vertices'] for c in calls)
    assert ops['source_provider']['navigation_links'] == sum(c['edges'] for c in calls)
    assert ops['native']['output_bits'] == sum(map(len, d['physical_read_service_outputs']))
    assert ops['native']['output_ends'] == ops['native']['acknowledgements'] == ops['native']['returns'] == len(calls)
    source = d['actual_original_source']
    assert source['preparation_events'] == 75*n+4
    assert source['persistent']['V'] == 120*n+8
    assert source['same_bound_root'] and source['current_grant'] and source['source_unchanged']
    account = d['resource_account']
    assert not account['complete_full_session_vector'] and not account['unrelated_component_peak_sum_used']
    for cut in list(account['named_cuts'].values())+[account['peak_common_cut']]:
        for key in ('objects','reference_slots','bit_fields','grant_bits','bits'):
            expected = sum(c[key] for c in cut['components'].values())
            if key == 'bits':
                expected += cut.get('overlap_bits', 0)
            assert cut['total'][key] == expected
        assert cut['total']['bits'] <= account['peak_common_cut']['total']['bits']
    expected = (HERE/'constructor'/(d['case']+'.B-input.bits')).read_text()
    assert required_loader(expected) == d['stages']['converted']['required_complete_loader']
    assert len(expected) == d['stages']['converted']['B_input_bits']
    physical = d['obstruction']['work']['plus_bits']
    consumed = len(expected)-d['obstruction']['remaining_input_chain_bits']
    assert len(physical) in (consumed, consumed+1)
    # One last split/mark may still be pending at the refusal. Every earlier
    # cell is already the input bit that the real loader consumed.
    assert physical[:-1] == [int(b) for b in expected[:len(physical)-1]]
    for name, identity in d['source_hashes'].items():
        assert file_identity(HERE/name) == identity
    return {'case': d['case'], 'source_calls': len(calls), 'receipt_payload_bits': sum(2*len(c['operand'])+3 for c in calls),
            'B_input_bits': len(expected), 'actual_remaining_input_bits': d['obstruction']['remaining_input_chain_bits'],
            'actual_plus_cells': len(d['obstruction']['work']['plus_bits']),
            'verified_committed_physical_loader_bits': len(physical)-1,
            'consumed_converted_bits': consumed,
            'actual_native_steps': count['native_steps'],
            'retained_common_cut': account['retained_common_cut']['total'],
            'peak_common_cut': account['peak_common_cut']['total'],
            'full_execution_complete': False}


def result_mutants(d):
    rejected_names = []
    for key in ('native_controller_selections', 'bit_field_reads'):
        mutant = copy.deepcopy(d)
        mutant['operations']['native'][key] -= 1
        rejected(lambda: validate_attempt(mutant))
        rejected_names.append(key)
    mutant = copy.deepcopy(d)
    mutant['receipt_occurrences'].pop()
    rejected(lambda: validate_attempt(mutant))
    rejected_names.append('missing_receipt_occurrence')
    mutant = copy.deepcopy(d)
    mutant['cache'][15][0] = 0
    rejected(lambda: validate_attempt(mutant))
    rejected_names.append('missing_original_index')
    mutant = copy.deepcopy(d)
    mutant['resource_account']['complete_full_session_vector'] = True
    rejected(lambda: validate_attempt(mutant))
    rejected_names.append('premature_full_session_completion')
    mutant = copy.deepcopy(d)
    mutant['obstruction']['work']['plus_bits'][0] ^= 1
    rejected(lambda: validate_attempt(mutant))
    rejected_names.append('wrong_loaded_physical_input_bit')
    return rejected_names


def native_mutants(kernel):
    results = []
    for kind in ('missing-index', 'missing-receipt'):
        p = service_graph()
        e = MeasuredEngine(p, 4000, kernel=kernel,
            connection=OriginalConnection(OriginalPort(prepare_source([0]*60))))
        selected = ('class.slot', 'class.slot', 1) if kind == 'missing-index' else ('receipts', 1, 'receipt.record')
        op = 'get' if kind == 'missing-index' else 'set'
        at = next(i for i,n in enumerate(p.nodes) if n.op == op and n.args == selected)
        e.code[at].opcode = OPS.index('goto')
        e.run()
        e.feed(u(4)+rational(0))
        operand = 'R'*15+'LL' if kind == 'missing-index' else 'LL'
        e.feed('01'+v(operand.translate(str.maketrans('LR','01'))))
        actual_calls = e.connection.port.calls
        if kind == 'missing-index':
            assert cache(e)[0] == (1,0) and cache(e)[15] == (0,0)
            detected = rejected(lambda: assert_index(cache(e)))
        else:
            assert not receipts(e) and len(actual_calls) == 1
            detected = rejected(lambda: assert_receipt(receipts(e), actual_calls))
        results.append({'mutation': kind, 'detected': detected, 'source_calls': len(actual_calls),
                        'actual_receipts': receipts(e), 'actual_cache': cache(e),
                        'actual_output': e.output, 'native_steps': e.native_steps,
                        'field_actions': e.ticks, 'installation_steps': e.installation.steps})
        del e
        import gc
        gc.collect()
    return results


def assert_index(paired):
    assert paired[15] == (1,0) and paired[0] == (0,0)


def assert_receipt(found, calls):
    assert found == [(a,b) for a,b,_,_ in calls]


def decoder_charge_falsifier(kernel):
    from copier import initial_word
    from nwmg import encode
    p = output_graph()
    e, _ = install_b((), encode(initial_word('0111')), 24000, p,
                     lambda *a, **kw: MeasuredNativeEngine(*a, kernel=kernel, **kw))
    e.run()
    assert ''.join(map(str,e.output)) == '0'
    assert e.boundaries['decode.output'] == 1
    validate_native(e.charges, e.native_events, e.native_steps, e.ticks)
    assert e.native_events['bfield'] > 0
    mutant = Counter(e.charges)
    mutant['bit_field_reads'] -= 1
    rejected(lambda: validate_native(mutant, e.native_events, e.native_steps, e.ticks))
    return {'scope': 'actual returned-head decoder component fixture; not evaluator output',
            'actual_output': e.output, 'field_actions': e.ticks, 'native_steps': e.native_steps,
            'actual_bit_field_reads': e.charges['bit_field_reads'],
            'one_omitted_read_rejected': True}


def premature_phase_falsifier(kernel):
    from field_protocol_checks import exchange_graph
    p = exchange_graph()
    e = MeasuredNativeEngine(p, 1, kernel=kernel)
    at = next(i for i,n in enumerate(p.nodes) if n.op == 'end')
    e.code[at].opcode = OPS.index('goto')
    try:
        e.run()
    except RuntimeError:
        assert not e.charges['output_ends']
        assert not e.charges['acknowledgements'] and not e.charges['returns']
        assert e.charges['native_refusals'] == 1
        return {'mutation': 'executing End replaced by goto',
                'detected': True, 'output_before_refusal': e.output,
                'End': e.charges['output_ends'],
                'Ack': e.charges['acknowledgements'], 'Return': e.charges['returns'],
                'native_steps': e.native_steps, 'field_actions': e.ticks,
                'scope': 'actual installed primitive phase guard component, not completed evaluator'}
    raise AssertionError('Return phase advanced without End/Ack')


def widening_falsifier():
    arena = Arena(registers=5)
    expected = 0
    widths = []
    original = arena.alloc
    def counted(kind):
        nonlocal expected
        before_width, slots = arena.width, arena.slots()
        cell = original(kind)
        if before_width != arena.width:
            expected += slots*(before_width+arena.width)
            widths.append((before_width, arena.width, slots))
        return cell
    arena.alloc = counted
    buf = Buffer(arena)
    buf.extend('0101'*40)
    buf.clear()
    assert arena.charges['reference_copy_bits'] == expected > 0
    mutant = Counter(arena.charges)
    mutant['reference_copy_bits'] -= 1
    rejected(lambda: equal_copy(mutant['reference_copy_bits'], expected))
    return {'width_changes': widths, 'actual_reference_copy_bits': expected,
            'omitted_copy_bit_rejected': True, 'cleanup_fields': arena.charges['cleanup_fields'],
            'free_list_writes': arena.charges['free_list_writes'], 'retained': arena.census(),
            'scope': 'actual Arena width/cleanup component; not complete evaluator session'}


def equal_copy(actual, expected):
    assert actual == expected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--theory-root', type=Path, required=True)
    ap.add_argument('--results', nargs='+', type=Path, required=True)
    ap.add_argument('--replay', nargs=2, type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    source_manifest(args.theory_root, HERE/'theory-candidate.json')
    ds = [json.loads(p.read_text()) for p in args.results]
    validation = [validate_attempt(d) for d in ds]
    replay = None
    if args.replay:
        a, b = [json.loads(p.read_text()) for p in args.replay]
        for d in (a,b):
            d.pop('seconds')
            d.pop('host_peak_rss_bytes')
            d['runtime'].pop('scope_descriptors_inherited')
        assert a == b
        replay = {'exact_after_removing_only_host_duration_RSS_and_scope_fd_count': True,
                  'results': [file_identity(p) for p in args.replay],
                  'same_primitive_counts_material_cuts_receipts_and_pending_state': True}
    kernel = BitKernel()
    try:
        mutants = native_mutants(kernel)
        decoder = decoder_charge_falsifier(kernel)
        phase = premature_phase_falsifier(kernel)
    finally:
        kernel.close()
    from integration_checks import source_binding_checks
    compatibility = source_binding_checks()
    reused = {}
    for name in ('service-results.json','service-domain-results.json','integration-results.json'):
        d = json.loads((HERE/name).read_text())
        for filename, identity in d['source_hashes'].items():
            assert hashlib.sha256((HERE/filename).read_bytes()).hexdigest() == identity
        reused[name] = {'result_identity': file_identity(HERE/name), 'matching_source_bindings': len(d['source_hashes']),
                        'scope': d.get('scope', 'supplied component or B author result; not fresh C acceptance')}
    regression = json.loads((HERE/'constructor/regression-results.json').read_text())
    assert [x['label_count'] for x in regression['DP_old_new_alias_falsifier']] == [0,2]
    interface = json.loads((HERE/'constructor/INTERFACE.json').read_text())
    for name, identity in interface['fixed_code'].items():
        assert file_identity(HERE/name) == identity
    result = {'schema':1, 'C_complete':False, 'D_complete':False, 'independent_acceptance':False,
        'whole_goal_complete':False, 'strict_missing_objects':STRICT,
        'actual_attempt_validation':validation, 'result_omission_mutants':result_mutants(ds[0]),
        'native_state_mutants':mutants, 'decoder_charge_falsifier':decoder,
        'premature_phase_falsifier':phase,
        'widening_charge_falsifier':widening_falsifier(), 'deterministic_replay':replay,
        'source_bound_reused_evidence':reused, 'fresh_native_compatibility':compatibility,
        'reused_DP_falsifier':regression['DP_old_new_alias_falsifier'],
        'complete_fixed_A_artifact_identities_verified':len(interface['fixed_code']),
        'source_identity':file_identity(Path(__file__)),
        'inputs':[{'path':str(p),**file_identity(p)} for p in args.results]}
    atomic(args.output, result)
    print(json.dumps({'author_prefix_checks':len(validation), 'C_complete':False, 'output':str(args.output)}))


if __name__ == '__main__':
    main()
