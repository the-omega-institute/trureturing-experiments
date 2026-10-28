#!/usr/bin/env python3
"""Exact independent checks for Finch arXiv:2610.34300v1, Remark 20.

Positive indices and the recurrence are taken literally from section 6.
This program supplies an integer certificate, not a Lean proof.
"""


def source_sequence(k, count):
    """Direct recurrence; x[0] is an unused sentinel."""
    x = [0] * (count + 1)
    if k <= count:
        x[k] = 1
    for n in range(k + 1, count + 1):
        x[n] = (x[n - 1] + x[n - k] ** 2) % 3
    return x


def pack(digits):
    s = 0
    for digit in digits:
        assert digit in (0, 1, 2)
        s = 3 * s + digit
    return s


def is_run(x, k, r, g):
    return r >= k + 1 and g >= 1 and all(
        (x[r + i] == 0) == (i % 4 == 0) for i in range(4 * g + 1)
    )


def suffix_has_zero_mask(s, g):
    for i in range(4 * g + 1):
        if ((s % 3) == 0) != (i % 4 == 0):
            return False
        s //= 3
    return True


def check_source_fidelity():
    # Example 18's entire displayed 60-term table, read from qrm3.tex.
    rows = [
        [0, 0, 0, 1, 1, 1, 1, 2, 0, 1],
        [2, 0, 0, 1, 2, 2, 2, 0, 1, 2],
        [0, 0, 1, 2, 2, 2, 0, 1, 2, 0],
        [0, 1, 2, 2, 2, 0, 1, 2, 0, 0],
        [1, 2, 2, 2, 0, 1, 2, 0, 0, 1],
        [2, 2, 2, 0, 1, 2, 0, 0, 1, 2],
    ]
    x = source_sequence(4, 64)
    assert x[1:61] == [v for row in rows for v in row]
    zeros = [n for n in range(5, 61) if x[n] == 0]
    assert zeros == [9, 12, 13, 18, 21, 22, 27, 30, 31, 36, 39, 40,
                     45, 48, 49, 54, 57, 58]
    assert all(not is_run(x, 4, r, 1) for r in range(5, 61))
    print('FIDELITY Example 18: all 60 terms and all 18 listed zeros match')


def check_counterexample():
    k = 61
    modulus = 3 ** (k - 1)
    state_bound = 3 ** k

    def transition(s):
        return 3 * (s % modulus) + (s % 3 + (s // modulus) ** 2) % 3

    # Discover the first repeated state rather than assuming the supplied collision.
    states, seen = [], {}
    s = 1
    while s not in seen:
        assert len(states) < 100000
        assert 0 <= s < state_bound
        seen[s] = len(states)
        states.append(s)
        s = transition(s)
    mu, stop = seen[s], len(states)
    assert (mu, stop, stop - mu) == (739, 9685, 8946)
    assert modulus == 42391158275216203514294433201
    assert s == states[739] == 92196811829875935758365998383
    assert states[2759] == 64673413418723788697046593874
    # The listed set is forward invariant, including its final transition.
    assert all(transition(states[j]) == states[j + 1] for j in range(stop - 1))
    assert transition(states[-1]) == states[mu]
    for gaps in (11, 15):
        count = sum(suffix_has_zero_mask(t, gaps) for t in states)
        assert count == 0
        print(f'CERTIFICATE k=61: {len(states)} states; {gaps}-gap masks={count}')

    # Independent implementation: direct scalar recurrence, without transition().
    period = stop - mu
    x = source_sequence(k, stop + k + 2 * period)
    for j, packed in enumerate(states):
        assert pack(x[j + 1:j + k + 1]) == packed
    assert pack(x[stop + 1:stop + k + 1]) == states[mu]
    assert all(x[n + period] == x[n] for n in range(mu + 1, len(x) - period))
    # Every r>=62 maps to suffix window S_(r-17); this checks representatives,
    # including starts near the state-cycle boundary and their continuations.
    representatives = range(k + 1, stop + 17)
    assert all(not is_run(x, k, r, 11) for r in representatives)
    assert all(not is_run(x, k, r, 15) for r in representatives)
    ten_gap_starts = [r for r in representatives if is_run(x, k, r, 10)]
    assert 2780 in ten_gap_starts
    assert x[2780:2821] == [0, 1, 2, 2] * 10 + [0]
    assert 3 ** 41 == 36472996377170786403
    assert states[2759] % (3 ** 41) == 7750511730148792110
    assert int('0122' * 10 + '0', 3) == 7750511730148792110
    assert (k - 1) // 4 == 15
    print(f'COLLISION S_739=S_9685={s}; first state-cycle length={period}')
    print('PERIOD FIDELITY: Finch (2002), k=61, period 8946 reproduced')
    print('DIRECT CHECK: all 9685 packed states agree with scalar recurrence')
    print(f'DIRECT CHECK: 10-gap representative starts={ten_gap_starts}')
    print('ATTAINMENT r=2780: (0,1,2,2)^10 followed by 0; suffix=7750511730148792110')
    print('RESULT k=61: exact maximal gap count 10 < source upper bound 15')
    print('INFINITE SCOPE: invariant state list + forbidden suffix covers the entire orbit')


def main():
    check_source_fidelity()
    check_counterexample()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
