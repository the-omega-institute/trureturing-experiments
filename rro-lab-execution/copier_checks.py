"""Executable contracts for binary-TM copying, B lowering and decoded LAB I/O."""
from __future__ import annotations

from itertools import product
from fractions import Fraction as F
from session import OriginalPort, Session, prepare_source
from copier import Copier, decode_records, initial_word
from nwmg import compile_non_erasing, encode
from lab import Arena, Machine, reference_run
from framing import Reader, Receive, Print, u, v


def tm_reference(rows, word, bound):
    tape = {i: int(b) for i, b in enumerate(word)}
    head = state = 0
    for _ in range(bound):
        if state == len(rows):
            break
        write, direction, state = rows[state][tape.get(head, 0)]
        tape[head] = write
        head += 1 if direction == 'R' else -1
    return {i for i, b in tape.items() if b}, head, state


def run_nm(copier, word, copies):
    tape = {i for i, b in enumerate(initial_word(word)) if b == '1'}
    head = state = steps = 0
    commits = []
    while len(commits) < copies and state != len(copier.rows):
        write, direction, state = copier.rows[state][int(head in tape)]
        assert write >= int(head in tape)
        if write:
            tape.add(head)
        head += 1 if direction == 'R' else -1
        steps += 1
        if state in copier.boundaries:
            decoded = decode_records(lambda p: int(p in tape), head)
            assert decoded['state'] is None
            decoded['state'] = copier.boundaries[state]
            commits.append((decoded, steps))
        assert steps < 1000000
    return tape, head, state, commits


def contracts():
    # Erasure, both moves, first/last head positions, halting and looping.
    machines = [
        (((0, 'L', 1), (0, 'L', 1)),),
        (((1, 'R', 1), (0, 'R', 1)),),
        (((1, 'L', 1), (0, 'R', 1)), ((0, 'R', 2), (0, 'L', 2))),
        (((1, 'R', 0), (0, 'L', 0)),),
    ]
    instances = steps = 0
    for rows in machines:
        c = Copier(rows)
        assert all(w >= b for row in c.rows for b, (w, _, _) in enumerate(row))
        for size in range(3):
            for bits in product('01', repeat=size):
                word = ''.join(bits)
                tape, h, q, commits = run_nm(c, word, 3)
                for index, (d, count) in enumerate(commits, 1):
                    expected = tm_reference(rows, word, index)
                    assert (set(d['marked']), d['head'], d['state']) == expected
                    assert d['copies'] == index
                    assert d['left'] == -index
                    assert d['right'] == max(1, len(word))-1+index
                    assert d['end']+8 == 8*(index+1)*(max(1, len(word))+index+2)
                assert commits
                steps += commits[-1][1]
                instances += 1
    stopped = Copier(())
    initial = initial_word('01')
    stopped_decoded = decode_records(lambda p: int(initial[p]), 0)
    assert stopped.rows == () and stopped_decoded['marked'] == [1]
    assert stopped_decoded['head'] == 0 and stopped_decoded['copies'] == 0
    c = Copier(machines[0])
    program = compile_non_erasing(c.rows)
    # Actually received finite input is copied/erased and decoded from LAB;
    # the model is a general erasing binary TM, not a law-evaluator fixture.
    arena = Arena()
    # Parsing and framing remain the separately identified paid-I/O model.
    receiver = Receive(arena)
    receiver.feed('00'+v('1'))
    assert receiver.seal() == ('00', '1')
    actual_input = Reader(receiver.buffer.bits()[2:]).v()
    m = Machine(program, encode(initial_word(actual_input)), arena, audit=False)
    ref_tape, ref_head, ref_state, commits = run_nm(c, actual_input, 1)
    limit = 8*commits[-1][1]+1
    while m.pc != len(program):
        assert m.commits < limit
        m.step()
    m.boundary()
    bits, bhead, _, _ = m.decode(charged=True)
    assert bhead == 2*ref_head
    d = decode_records(lambda pos: bits.get(2*pos+1, 0), bhead//2)
    assert d['head'] == -1 and d['marked'] == []
    # Explicit finite return interval, read from the actual last complete copy.
    out = ''.join(str(int(j in d['marked'])) for j in range(d['left'], d['right']+1))
    wire = Print(arena, v(out)).finish()
    decoded = Reader(wire)
    assert decoded.v() == '000'
    decoded.end()
    receiver.buffer.clear()
    # A second job in the same arena must have its own persistent census and
    # per-job counters; previous work and routing remain allocated.
    second = Machine(('M',), '', arena)
    second.run(1)
    m.boundary()
    assert not (m.lab.owned & second.lab.owned)
    return {'instances': instances, 'non_erasing_steps': steps,
            'states': len(c.rows), 'B_code_length': len(program),
            'LAB_erasing_input': actual_input, 'LAB_decoded_output': out,
            'LAB_B_commits': m.commits, 'LAB_NM_steps': commits[-1][1],
            'LAB_persistent': m.lab.persistent(), 'shared_arena_job_census': True,
            'partial_retained': arena.census(), 'partial_peak_bits': arena.peak_bits,
            'immediately_halted_source': True,
            'full_evaluator': False, 'full_resource_account': False}


def kernel_session():
    """An acquired integer successor and real original-port continuation.

    This is an ML-T arithmetic *kernel*, not the complete ML-A evaluator.
    Its counters are all on the copied TM tape; framing/decoder/driver backing
    remains outside this partial ledger and cannot certify full ML.16.
    """
    source = prepare_source([1]*60)
    inherited = source.a.charges.copy()
    port = OriginalPort(source)
    session = Session(4, F(1, 3), port)
    session.read('LL')
    acquired = str(dict(session.cache_tuple()[0])[0])
    assert acquired == '1'
    # Little-endian successor: erase carries, then mark the first zero.
    c = Copier((((1, 'L', 1), (0, 'R', 0)),))
    program = compile_non_erasing(c.rows)
    receiver = Receive(session.a)
    receiver.feed('00'+v(acquired))
    assert receiver.seal() == ('00', acquired)
    received = Reader(receiver.buffer.bits()[2:]).v()
    m = Machine(program, encode(initial_word(received)), session.a, audit=False)
    _, _, _, commits = run_nm(c, received, 2)
    while m.pc != len(program):
        assert m.commits < 8*commits[-1][1]+1
        m.step()
    m.boundary()
    values, head, _, _ = m.decode(charged=True)
    d = decode_records(lambda p: values.get(2*p+1, 0), head//2)
    digits = ''.join(str(int(j in d['marked'])) for j in range(2))
    returned = Print(session.a, v(digits)).finish()
    r = Reader(returned)
    assert r.v() == '01'  # Little-endian two, obtained from the actual copied tape.
    r.end()
    receiver.buffer.clear()
    frozen = dict(m.lab.persistent())
    session.read('LL')
    session.read('LR')
    session.read('LL')
    assert session.port is port and port.root is source.root
    assert [w for w, _, _, _ in port.calls] == ['LL', 'LL', 'LR', 'LL']
    assert dict(session.cache_tuple()[0]) == {0: 1}
    assert m.lab.persistent() == frozen
    m.boundary()
    return {'acquired_original_index_coordinate': [0, 0],
            'input_little_endian': acquired, 'decoded_little_endian': digits,
            'NM_steps': commits[-1][1], 'B_commits': m.commits,
            'copier': c.description, 'work': frozen,
            'same_original_port': True, 'receipts': session.check_receipts(),
            'observer_partial_charges': dict(sorted(session.a.charges.items())),
            'provider_primitive_charges': dict(sorted((source.a.charges-inherited).items())),
            'inherited_preparation_charges': dict(sorted(inherited.items())),
            'source_retained': source.a.census(),
            'observer_partial_retained': session.a.census(),
            'simultaneous_partial_peak_bits': source.a.census()['bits']+session.a.peak_bits,
            'scope': 'Acquired successor kernel only; full evaluator, code/control and I/O/decoder backing not certified.',
            'full_evaluator': False, 'full_resource_account': False}


if __name__ == '__main__':
    import sys
    failures = 0
    for name, run in [('copier-contracts', contracts), ('acquired-successor-session', kernel_session)]:
        try:
            result = run()
        except AssertionError as error:
            print('[FAIL]', name, repr(error))
            failures += 1
        else:
            print('[PASS]', name, result)
    raise SystemExit(int(bool(failures)))
