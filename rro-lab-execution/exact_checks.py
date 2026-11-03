"""Reproducible exact finite checks of sealed ML definitions; Python >=3.12.

Exit 0: implemented checks pass. --require-all exits 2 for explicit unexecuted
acceptance obligations, or 1 for an assertion/source-correspondence failure.
No successful host arithmetic check is presented as a compiled LAB evaluator.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

from framing import (BadFrame, NeedMore, Reader, Receive, Print, REPLIES, arithmetic_input,
                     arithmetic_output, parse_arithmetic, parse_arithmetic_output,
                     parse_reply, parse_command, rational, reply, u, v, words)
from lab import Arena, Buffer, Lab, Machine, reference_run
from copier_checks import contracts as copier_contracts, kernel_session
from laws import Evaluator, capacity_witness, conjoin, labels, repeated_full_direct
from nwmg import compile_non_erasing, encode
from pins import (BASE, EXPERIMENT_PIN, POSITIVE_PIN, ROOT, SOURCES, SUPPLIER_PIN,
                  THEORY_PATH, THEORY_PIN, source_manifest, suppliers)
from session import OriginalPort, Session, classify, prepare_source

assert sys.version_info >= (3, 12)
assert __debug__, 'Assertions must be enabled.'
HERE = Path(__file__).resolve().parent


def refused(call):
    try:
        call()
    except (AssertionError, BadFrame):
        return True
    raise AssertionError('Forbidden transition unexpectedly succeeded')


def compare(program, bits, bound):
    m = Machine(program, bits)
    ref, locations = reference_run(program, bits, bound)
    m.run(bound)
    assert [(t['head'], t['pc'], t['marked']) for t in m.trace] == ref
    assert m.commits == len(locations)
    return m


def structural_views(lab):
    """AA17.3 on actual physical image: support <-> boundary <-> logical tree."""
    support = set()
    def visit(node, path):
        support.add(path)
        for field, side in (('left', 'L'), ('right', 'R')):
            child = node.fields[field]
            if child is not None:
                visit(child, path+side)
    visit(lab.root, '')
    boundary = {p for p in support if p+'L' not in support and p+'R' not in support}
    closed = {word[:j] for word in boundary for j in range(len(word)+1)}
    assert closed == support
    def parse(path):
        if path+'L' not in support:
            assert path+'R' not in support
            return 'alpha'
        if path+'R' not in support:
            assert path+'LL' not in support and path+'LR' not in support
            return 'beta'
        return parse(path+'L'), parse(path+'R')
    assert parse('') == lab.logical()
    return {'support_vertices': len(support), 'terminal_boundary': len(boundary)}


def execution_checks():
    specs = {
        'MLX-cross': (('L','M','R','M','L','L','M','R','R','M'), '', 10),
        'MLX-branch': ((('J',3),'M',('J',4),'R','M'), '', 10),
        'MLX-frontier': (('R','R','L','L','L','L','R'), '', 7),
        'MLX-loop': (('M',('J',1)), '', 65),
        'MLX-bare': (('R','M'), '00', 2),
        'MLX-input': ((), '0101', 0),
    }
    output = {}
    for name, (program, bits, bound) in specs.items():
        m = compare(program, bits, bound)
        last = m.trace[-1]
        output[name] = {'commits': m.commits, 'halted': m.pc == len(program),
                        'final': last, 'views': structural_views(m.lab),
                        'routing': m.lab.routing(), 'charges': dict(sorted(m.lab.a.charges.items())),
                        'retained': m.lab.a.census(), 'peak_bits': m.lab.a.peak_bits,
                        'internal_cuts': len(m.lab.cuts)}
        if name == 'MLX-cross':
            assert [last[k] for k in ('P','S','V','H','D')] == [5,4,14,4,12]
            assert last['objects'] == 30 and last['reference_slots'] == 76
            assert last['head'] == 0 and last['marked'] == [-2,-1,0]
        elif name == 'MLX-branch':
            assert reference_run(program, bits, bound)[1] == [0,1,2,4]
            assert m.commits == 4 and last['P'] == 3
        elif name == 'MLX-frontier':
            assert [last[k] for k in ('P','S','V','H','D')] == [2,6,15,6,10]
        elif name == 'MLX-loop':
            assert m.pc == 1 and not output[name]['halted']
            assert len({(t['P'],t['S'],t['D']) for t in m.trace[1:]}) == 1
        elif name == 'MLX-input':
            assert m.output_bits(-1,4) == '001010'
            assert (last['P'], last['S']) == (4,5)
            # Actual halted tape observation is framed and physically printed.
            payload = words([u(m.pc), '0'+u(0), m.output_bits(-1,4)])
            flow = Print(m.lab.a, payload)
            assert flow.finish() == payload
            output[name]['framed_output'] = payload
    bare = Machine(('R','M'), '00')
    first = (bare.lab.logical(), bare.lab.routing(), structural_views(bare.lab))
    bare.step()
    assert first == (bare.lab.logical(), bare.lab.routing(), structural_views(bare.lab))
    assert bare.trace[0]['head'] == 0 and bare.trace[1]['head'] == 1
    assert bare.program[0] == 'R' and bare.program[bare.pc] == 'M'
    # Refusals preserve the actual tree and event word, including wrong guards.
    guard = compare(('M',), '', 1)
    previous = guard.lab.logical(), tuple(guard.lab.events), guard.lab.persistent()
    assert not guard.lab.grow('P', guard.head.fields['right'])
    assert not guard.lab.grow('S', guard.head)
    assert previous == (guard.lab.logical(), tuple(guard.lab.events), guard.lab.persistent())
    refused(lambda: guard.output_bits(0, -1))
    # Actual edit-before-maintenance cut retains routing obligation and live material.
    cuts = [c for c in guard.lab.cuts if c['phase'] == 'EditComplete']
    assert cuts and all(c['pending'] for c in cuts)
    later = [c for c in guard.lab.cuts if c['phase'] == 'AppendInitialized']
    assert later and all(c['pending'] for c in later)
    return output


def exhaustive_b_checks():
    counts = Counter()
    inputs = [''.join(x) for n in range(3) for x in product('01', repeat=n)]
    for q in range(4):
        instructions = ('L','R','M')+tuple(('J',j) for j in range(q))
        for program in product(instructions, repeat=q):
            for bits in inputs:
                m = compare(program, bits, 8)
                structural_views(m.lab)
                counts['instances'] += 1
                counts['instruction_commits'] += m.commits
                counts['logical_events'] += len(m.lab.events)
                counts['primitive_pending_cuts'] += len(m.lab.cuts)
                counts['halted' if m.pc == len(program) else 'nonhalting_boundaries'] += 1
    assert counts['instances'] == 1722
    return {'program_lengths': [0,1,2,3], 'input_lengths': [0,1,2], 'bound': 8,
            'all_jump_targets': True, **dict(counts)}


def compiler_checks():
    count, commits, variants = 0, 0, []
    # Exercise all six rule shapes, both initial symbols and existing/new neighbors.
    for write0, dir0, dir1 in product((0,1), ('L','R'), ('L','R')):
        transitions = (((write0, dir0, 1), (1, dir1, 1)),)
        w = compile_non_erasing(transitions)
        assert len(w) == 14
        for n in (1,2,3):
            for raw in product((0,1), repeat=n):
                bits = ''.join(map(str,raw))
                m = compare(w, encode(raw), 9)
                assert m.pc == len(w)
                values, head, _, _ = m.decode()
                read = raw[0]
                write, direction, _ = transitions[0][read]
                target = -1 if direction == 'L' else 1
                expected = {i: b for i,b in enumerate(raw)}
                expected[0] = write
                assert head == 2*target
                for j in range(-1, max(2,n)):
                    assert values.get(2*j+1,0) == expected.get(j,0)
                assert values.get(head,0) == 1
                assert 7 <= m.commits <= 9  # <=8 per rule, then terminal M.
                count += 1
                commits += m.commits
        variants.append([write0,dir0,dir1])
    # A multistate finite run crosses both tape sides before its terminal M.
    t = (((1,'R',1),(1,'R',1)), ((0,'L',2),(1,'L',2)),
         ((0,'L',3),(1,'L',3)), ((1,'R',4),(1,'R',4)))
    m = compare(compile_non_erasing(t), encode([0]), 40)
    assert m.pc == len(m.program) and m.trace[-1]['head'] == 0
    return {'source': 'NWMG arXiv:1304.0053v2 pp.6-9 equations (1)-(3)',
            'pdf_sha256': '5487d3f1bf8400cabc203153e3dc60b2ae133b3112ad7299c11fc3fc37e0112e',
            'single_transition_instances': count, 'rule_table_variants': variants,
            'single_transition_B_commits': commits,
            'multistate_B_commits': m.commits,
            'limit': 'Non-erasing binary transition compiler only; ML-T copier/evaluator not executed.'}


def framing_checks():
    count = 0
    for length in range(9):
        for bits in product('01', repeat=length):
            word = ''.join(bits)
            framed = v(word)
            r = Reader(framed)
            assert r.v() == word
            r.end()
            for end in range(len(framed)):
                try:
                    Reader(framed[:end]).v()
                except NeedMore:
                    pass
                else:
                    raise AssertionError('proper V prefix completed')
            count += 1
    for n in range(257):
        r = Reader(u(n))
        assert r.u() == n
        r.end()
    for value in (F(0),F(1),F(-1),F(1,3),F(-123,1024),F(2**80+1,3**40)):
        r = Reader(rational(value))
        assert r.rational() == value
        r.end()
    mask_cases = 0
    for mask in range(1 << 15):
        for parity in (0,1):
            data = (tuple((h,(h+parity) % 2) for h in range(15) if mask & (1 << h)),(),(),())
            frame = arithmetic_input(4,F(1,3),data,[])
            assert parse_arithmetic(frame) == (4,F(1,3),data,[])
            mask_cases += 1
    # A real stream prefix waits; there is no response from unreceived suffix bits.
    a = Arena()
    receive = Receive(a)
    command = '01'+v('0101')
    for bit in command[:-1]:
        receive.feed(bit)
        assert receive.seal() is None
    receive.feed(command[-1])
    assert receive.seal() == ('01','0101')
    refused(lambda: parse_command('11'+'0'))
    refused(lambda: Reader(v('00')).u())
    refused(lambda: Reader('0'+u(1)+u(0)).rational())
    for kind, payload in ((0,''),(1,words([rational(F(1,3))])),(2,'11'),(3,'')):
        wire = reply(kind,payload)
        assert parse_reply(wire) == (kind,payload)
        flow = Print(a,wire)
        assert flow.bit()
        assert flow.phase == 'Print' and flow.cursor is not None
        refused(flow.ack)
        refused(flow.end)
        assert flow.finish() == wire
    return {'V_words': count, 'U_values': 257, 'rational_values': 6,
            'within_slot_mask_frames':mask_cases,'within_slot_masks':1 << 15,
            'selected_value_patterns_per_mask':2,
            'proper_frame_prefixes_wait': True, 'reply_kinds': 4,
            'premature_End_Ack_refused': True, 'charges': dict(sorted(a.charges.items())),
            'retained': a.census(), 'peak_bits': a.peak_bits}


def arithmetic_checks(old, general, evaluator):
    empty = tuple(() for _ in range(4))
    one = old.conjoin(empty,0,0,0)
    left = (tuple(old.complete(0,0)[:7]),(),(),())
    same = (left[0],left[0],(),())
    distinct = (left[0],tuple(old.complete(1,0)[:7]),(),())
    full = tuple(old.complete(0,0) for _ in range(4))
    expected = {'empty': F(1), 'one_bit': F(1,2), 'one_left': F(1,128),
                'same_left_indices': F(1,72*128),
                'distinct_left_indices': F(71,72*128*127),
                'full_repeated_pair': F(1,1944**3*128*256)}
    cases = dict(zip(expected, (empty,one,left,same,distinct,full)))
    seen = {}
    for name, data in cases.items():
        z = evaluator.normalizer(data)
        assert z == (expected[name],)*2
        seen[name] = {'Z0':str(z[0]),'Z1':str(z[1])}
    direct, direct_counts = repeated_full_direct(old,full)
    assert direct == evaluator.normalizer(full) == general.normalizer(full)
    square = tuple(old.complete(a,b) for a in (0,1) for b in (0,1))
    square_z = evaluator.normalizer(square)
    literal_lifts = []
    for i in (0,1):
        rows = old.literal(i)
        lift = sum(len(rows[a] & rows[b])*(len(rows[a] & rows[b])-1)
                   for a in range(72) for b in range(72) if a != b)
        literal_lifts.append(lift)
        assert square_z[i] == F(lift,1944**4*128*127*256*255)
    assert square_z[0] != square_z[1] and min(square_z) > 0
    assert general.normalizer(square) == square_z
    prefix_square = square[:3]+((),)
    square_mixtures = {}
    for theta in (F(0),F(1,7),F(1,3),F(1)):
        values,ratios = evaluator.evaluate(4,theta,prefix_square,[square])
        assert ratios is not None and 0 < ratios[0] < 1
        assert values[5] == (1-theta)*square_z[0]+theta*square_z[1]
        square_mixtures[str(theta)] = {'joint':str(values[5]),'conditional':str(ratios[0])}
    # A genuine zero-normalizer input: both complete sides force K_(2,2), which
    # is positive, whereas 73 distinct complete left labels exceed capacity.
    data73 = tuple(tuple(old.complete(j,0)[:7]) for j in range(73))
    assert evaluator.normalizer(data73) == (F(0),F(0))
    capacity, capacity_data = capacity_witness(old)
    # Do not enumerate Bell(30)^2 or misreport that this shortcut ran ML-A.
    capacity_array = tuple(old.complete(a,b) for a in range(3) for b in range(10))
    for theta in (F(0),F(1,7),F(1,3),F(1)):
        assert (1-theta)*capacity[0]+theta*capacity[1] == theta*capacity[1]
    # Test the exact old-array update independently at small widths against
    # direct injective label assignments, then against 7/8-bit suppliers.
    dp_cases = 0
    for width in (1,2,3):
        masks = [(), ((0,0),), ((0,1),)]
        for r in range(1,4):
            for chosen in product(masks, repeat=r):
                direct_count = sum(len(set(values)) == r and all(
                    all((value >> (width-1-h)) & 1 == bit for h,bit in restrictions)
                    for value,restrictions in zip(values,chosen))
                    for values in product(range(1 << width),repeat=r))
                assert labels(chosen,width) == direct_count
                dp_cases += 1
    for width in (7,8):
        for chosen in (((),), ((),()), (((0,0),),((0,1),)), (None,()),
                       (tuple((h,0) for h in range(width)),)*2):
            assert labels(chosen,width) == old.label_dp(chosen,width)
            dp_cases += 1
    # Mixed masks at all 15 original coordinates; same-index vs distinct-index
    # and shared injection laws are compared by the supplier on its exact domain.
    priors = (F(0),F(1,7),F(1,3),F(1))
    framed = 0
    for n in (4,5,7):
        for h in range(15):
            data = tuple(((h,h % 2),) if j == 0 else () for j in range(n))
            same_c = tuple(((h,h % 2),) if j == 0 else () for j in range(n))
            contrary = tuple(((h,1-h % 2),) if j == 0 else () for j in range(n))
            other = tuple(((h,h % 2),) if j == 1 else () for j in range(n))
            for theta in priors:
                word = arithmetic_input(n,theta,data,[same_c,contrary,other,None])
                args = parse_arithmetic(word)
                assert args == (n,theta,data,[same_c,contrary,other,None])
                values, ratios = evaluator.evaluate(*args)
                assert ratios[:2] == [F(1),F(0)] and ratios[3] == 0
                population, alphabet = (72,128) if h < 7 else (243,256)
                assert ratios[2] == F(1,population)+(1-F(1,population))*F(alphabet//2-1,alphabet-1)
                out = arithmetic_output(values,ratios)
                flow = Print(Arena(),reply(1,out))
                returned = flow.finish()
                kind,payload = parse_reply(returned)
                assert kind == 1 and parse_arithmetic_output(payload) == (values,ratios)
                framed += 1
    values, ratios = evaluator.evaluate(73,F(1,3),data73,[data73])
    assert ratios is None and values == [F(0)]*6
    assert parse_arithmetic_output(arithmetic_output(values,ratios)) == (values,ratios)
    # Full-input retained denominator test: same data at N=4 and N=5 differs.
    five = tuple(old.complete(0,0) for _ in range(5))
    assert evaluator.normalizer(five)[0] == expected['full_repeated_pair']/1944
    return {'closed_form_cases':seen, 'ML_A_direct_N4':direct_counts,
            'ML_A_direct_result':list(map(str,direct)), 'label_DP_checks':dp_cases,
            'framed_host_arithmetic_roundtrips':framed+1,
            'priors':list(map(str,priors)), 'partial_mask_N':[4,5,7],
            'capacity':capacity_data, 'capacity_array_sha256':hashlib.sha256(
                arithmetic_input(30,F(1,3),capacity_array,[]).encode()).hexdigest(),
            'zero_capacity_N':73, 'normalizer_checks':dict(evaluator.counts),
            'positive_source_distinction_N4':{'literal_C4_injective_lifts':literal_lifts,
                                               'Z0':str(square_z[0]),'Z1':str(square_z[1]),
                                               'prior_mixtures':square_mixtures},
            'limit':'Host exact arithmetic and output flow, not ML-T/NWMG/ML-B evaluator execution.'}


def port_menu_checks(old, positive):
    counts = Counter()
    arrays = [tuple([0]*60),tuple([1]*60)]
    for i in (0,1):
        edges = [(a,t) for a,row in enumerate(old.literal(i)) for t in sorted(row)]
        samples = [edges[0],edges[0],edges[1],edges[-1]]
        arrays.append(tuple(bit for a,t in samples for _,bit in old.complete(a,t)))
    menu = set(''.join(p) for n in range(8) for p in product('LR',repeat=n))
    for b in range(60):
        for length in range(5):
            menu.update('R'*b+'L'+''.join(p) for p in product('LR',repeat=length))
    menu.update('R'*60+''.join(p) for length in range(6) for p in product('LR',repeat=length))
    menu.update(('R'*130,'R'*60+'L'*20,'LL'+'R'*100))
    observations = []
    for bits in arrays:
        source = prepare_source(bits)
        original = source.logical(),tuple(source.events),source.routing()
        oracle = positive.tree(bits)
        port = OriginalPort(source)
        replies = Counter()
        for word in sorted(menu):
            a = Arena()
            before = source.a.charges['field_reads']
            latch = port.call(word,a)
            actual = REPLIES[int(latch.bits(),2)]
            expect = positive.read(oracle,tuple(int(x == 'R') for x in word))
            assert actual == expect
            bit,pair = classify(4,word)
            assert actual == pair[0 if bit is None else bits[bit]]
            e = port.calls[-1][2]
            assert source.a.charges['field_reads']-before == 2*e
            counts['calls'] += 1
            counts['operand_bits'] += len(word)
            counts['provider_vertices'] += e
            replies[actual] += 1
        assert original == (source.logical(),tuple(source.events),source.routing())
        assert set(replies) == set(REPLIES)
        observations.append({'z_sha256':hashlib.sha256(bytes(bits)).hexdigest(),
                             'replies':dict(sorted(replies.items())),
                             'source_P':source.a.charges['P'], 'source_S':source.a.charges['S'],
                             'source_persistent':source.persistent(),
                             'preparation_events':len(source.events)})
    return {'N':4,'fixtures':observations,'distinct_operands_per_fixture':len(menu),
            'all_slot_suffixes_through_length':4,'marker_suffixes_through_length':5,
            'all_root_words_through_length':7,'maximum_operand_length':max(map(len,menu)),
            **dict(counts)}


def session_checks(old, evaluator):
    observations = []
    for bit in (0,1):
        bits = [bit]*60
        source = prepare_source(bits)
        preparation_charges = Counter(source.a.charges)
        port = OriginalPort(source)
        session = Session(4,F(1,3),port)
        operands = ['','LL','LRR','LLR','LL','LR','LLL','LRL']+[
            'R'*(b-1)+'LL' for b in range(2,16)]+['LL','R'*15+'LL','R'*59+'LL','R'*60+'LLL']
        predictions, acquired = [], set()
        for word in operands:
            before = session.cache_tuple()
            target, pair = classify(4,word)
            cylinders = []
            for response in REPLIES:
                if response not in pair:
                    c = None
                elif target is None:
                    c = tuple(() for _ in range(4))
                else:
                    c = [() for _ in range(4)]
                    c[target//15] = ((target % 15,pair.index(response)),)
                    c = tuple(c)
                cylinders.append(c)
            calls_before = len(port.calls)
            _, probs = evaluator.evaluate(4,F(1,3),before,cylinders)
            assert len(port.calls) == calls_before and probs is not None and sum(probs) == 1
            # The actual return comes only after the unchanged real port call.
            session.start(word)
            old_receipts = len(session.records)
            old_cache = session.cache_tuple()
            refused(session.finish)
            session.capture()
            assert len(session.records) == old_receipts and session.cache_tuple() == old_cache
            assert session.pending['handle'] is port
            refused(lambda:session.start('LL'))
            session.append()
            session.update()
            framed = session.finish()
            kind,payload = parse_reply(framed)
            actual = REPLIES[int(payload,2)]
            assert kind == 2 and actual == pair[0 if target is None else bit]
            assert probs[REPLIES.index(actual)] > 0
            if target in acquired:
                assert probs[REPLIES.index(actual)] == 1
            if target is not None:
                acquired.add(target)
            predictions.append([str(x) for x in probs])
        assert len(dict(session.cache_tuple()[0])) == 15
        assert len(session.cache_tuple()[1]) == 1 and len(session.cache_tuple()[3]) == 1
        assert len(session.cache_tuple()[2]) == 0
        assert sum(w == 'LL' for w,_,_ in session.records) == 3
        assert len({id(b.header) for _,_,b in session.records}) == len(operands)
        unlinked = [x for x in session.cuts if x['phase'] == 'ReceiptUnlinked']
        linked = [x for x in session.cuts if x['phase'] == 'ReceiptLinked']
        captured = [x for x in session.cuts if x['phase'] == 'Captured']
        assert len(unlinked) == len(linked) == len(captured) == len(operands)
        assert all(x['pending'] for x in unlinked+linked+captured)
        observations.append({'constant_source_bit':bit,'original_calls':len(port.calls),
                             'original_operands':operands,
                             'actual_replies':[val for _,val,_,_ in port.calls],
                             'LL_occurrences':3,'receipt_census':session.check_receipts(),
                             'final_cache':[dict(row) for row in session.cache_tuple()],
                             'mathematical_predictions':predictions,
                             'charges':dict(sorted(session.a.charges.items())),
                             'retained':session.a.census(), 'peak_bits':session.a.peak_bits,
                             'pending_cuts':len(session.cuts),
                             'inherited_source_persistent':source.persistent(),
                             'inherited_preparation_charges':dict(sorted(preparation_charges.items())),
                             'provider_primitive_charges':dict(sorted((source.a.charges-preparation_charges).items())),
                             'provider_observations':dict(sorted(port.observations.items()))})
    return {'N':4,'theta':'1/3','fixtures':observations,
            'limit':'Predict values corroborated by exact host evaluator; no actual compiled evaluator job claimed.'}


def null_law_session(old):
    # The promise C_N includes arrays outside both probability supports.
    # Nine distinct left labels with one common right label force K_(9,1),
    # while every literal right endpoint has degree eight in both sources.
    n = 9
    data = tuple(old.complete(a,0) for a in range(n))
    for i in (0,1):
        rows = old.literal(i)
        columns = [sum(t in row for row in rows) for t in range(243)]
        assert set(columns) == {8}
        assert sum(old.falling(degree,9) for degree in columns) == 0
    bits = [value for row in data for _,value in row]
    source = prepare_source(bits)
    preparation_charges = Counter(source.a.charges)
    port = OriginalPort(source)
    session = Session(n,F(0),port)
    for b in range(15*n):
        session.read('R'*b+'LL')
    assert session.cache_tuple() == data
    # Null law is not permission to fabricate or suppress the actual old reply.
    before = len(port.calls)
    output = session.read('LL')
    assert len(port.calls) == before+1 and parse_reply(output) == (2,'01')
    null_output = arithmetic_output([F(0)]*3,None)
    assert parse_arithmetic_output(null_output) == ([F(0)]*3,None)
    return {'N':n,'literal_forced_relation':'K_(9,1)', 'prior':'0',
            'Z0':'0','Z1':'0','actual_calls':len(port.calls),
            'original_indices_preserved':n,'cache_selected_bits':15*n,
            'same_original_port':session.port is port,
            'receipt_census':session.check_receipts(),
            'retained':session.a.census(),'peak_bits':session.a.peak_bits,
            'source_persistent':source.persistent(),
            'inherited_preparation_charges':dict(sorted(preparation_charges.items())),
            'provider_primitive_charges':dict(sorted((source.a.charges-preparation_charges).items())),
            'session_and_call_charges':dict(sorted(session.a.charges.items())),
            'source_events':len(source.events),
            'limit':'All-promise real-port witness outside both law supports; null output is host-built, not a compiled job.'}


def build_result(theory_root, candidate_manifest):
    source_manifest(theory_root, candidate_manifest)
    old, general, positive = suppliers()
    evaluator = Evaluator(old,general)
    result = {'schema':2, 'sealed_theory_pin':THEORY_PIN,
              'theory_sha256':SOURCES['candidate:'+THEORY_PATH]['sha256'],
              'experiment_base':EXPERIMENT_PIN, 'python_minimum':'3.12',
              'execution':execution_checks(), 'exhaustive_B':exhaustive_b_checks(),
              'non_erasing_compiler':compiler_checks(),
              'explicit_copier':copier_contracts(),
              'acquired_arithmetic_kernel_session':kernel_session(),
              'framing':framing_checks(),
              'arithmetic':arithmetic_checks(old,general,evaluator),
              'original_port_menu':port_menu_checks(old,positive),
              'real_sessions':session_checks(old,evaluator),
              'null_law_real_session':null_law_session(old)}
    result['source_manifest'] = dict(sorted(SOURCES.items()))
    result['mismatches'] = []
    result['unresolved_obligations'] = [
        {'id':'ML-N-evaluator-execution','interface':'Definitions 4.4/4.6; Theorem 4.5; ML-X arithmetic',
         'status':'not-executed',
         'finding':'The copier now generates a complete non-erasing binary table for every supplied erasing binary table, lowers it through the primary thirteen-instruction compiler, and decodes real LAB output. The acquired successor kernel executes on this path. ML-A still has no instantiated complete ML-T evaluator table, hence no complete evaluator W_U or decoded full-law run. Exact host law values and an arithmetic kernel do not discharge that obligation.',
         'required':'Run a public-input ML-A instance through the actual complete ML-T/NWMG/ML-L/ML-B path and decode/print its exact result; preserve every N/mask/prior domain.'},
        {'id':'ML-Cost-full-evaluator-session','interface':'Theorems 6.4/6.6; ML.16',
         'status':'not-executed',
         'finding':'Per-root work/routing censuses now remain valid across sequential jobs in one arena. An acquired successor kernel, its real original-port continuation and a simultaneous partial storage peak execute together. The complete evaluator job, installed finite-driver/code/control, and every framing/conversion/decoder/counter backing are still not instantiated. These partial vectors cannot certify complete ML.16.',
         'required':'Census one actual full compiled arithmetic and same-port continuation run, including installed control/code, backing for every counter/cursor/phase, conversion, decoder, buffers, grants, cleanup and width changes.'},
    ]
    result['acceptance_complete'] = False
    result['limits'] = [
        'Finite exact checks corroborate ordinary proofs; no infinite-domain completeness or Lean verification.',
        'P/S and provider are simulations of explicit paid ML-I premises, not physical allocation or hidden-source access.',
        'No operational minimality, real serial-history recovery, C117 quantum frontend or source reset/fork is claimed.',
        'Host arithmetic, bit operations, abstract B commits and LAB material/navigation coordinates remain distinct.',
        'N30/N73 capacities are independent exact combinatorial checks; Bell(30)^2/Bell(73)^2 ML-A schedules were not run.',
    ]
    result['program_hashes'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(HERE.glob('*.py'))}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--theory-root',type=Path,required=True,
                        help='Actual candidate checkout; both authorized candidate files must match the content manifest.')
    parser.add_argument('--candidate-manifest',type=Path,default=HERE/'theory-candidate.json')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--check',type=Path)
    parser.add_argument('--require-all',action='store_true')
    args = parser.parse_args()
    result = build_result(args.theory_root, args.candidate_manifest)
    data = json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+'\n'
    if args.output:
        args.output.write_text(data)
    if args.check:
        assert args.check.read_text() == data, 'Saved result differs from exact rerun'
    print(json.dumps({'checks':'pass','acceptance_complete':result['acceptance_complete'],
                      'unresolved_obligations':[x['id'] for x in result['unresolved_obligations']],
                      'result_sha256':hashlib.sha256(data.encode()).hexdigest()},sort_keys=True))
    if args.require_all and not result['acceptance_complete']:
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
