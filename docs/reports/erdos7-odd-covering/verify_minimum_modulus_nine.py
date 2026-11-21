#!/usr/bin/env python3
"""Replay the finite prefix for the least-modulus-nine theorem."""

from fractions import Fraction as F
from math import isqrt

T = 6997
S = 10**24
BANDS = tuple(zip(
    (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 101, 211, 419, 839, 1679, 3359, T),
    (0, 0, 0, 12, 18, 20, 18, 24, 26, 24, 20, 28, 26, 26, 28, 28, 30, 30, 30, 30, 28, 26),
))
CASES = (
    {5: (1,), 7: (1,), 11: (1,)},
    {5: (1, 2), 7: (1,)},
)


def ceil(x):
    return -(-x.numerator // x.denominator)


def prime(p):
    return p > 1 and all(p % d for d in range(2, isqrt(p) + 1))


def evaluate(case):
    width = 1850
    probabilities = [[0] * (width + 1) for _ in range(3)]
    probabilities[0][1] = S
    mean = moment = S
    mass_debit = 0
    for p in range(3, T + 1, 2):
        if not prime(p):
            continue
        delta = F(next(n for end, n in BANDS if p <= end), 100)
        d = p - 1 - (F(1, p) if p in (3, 5, 7) else 1)
        cap = (p - 1) / (d * (1 - delta))
        if not 0 < cap <= p:
            raise ValueError(f"invalid comparison cap at {p}")
        thresholds = tuple(1 + delta*d + F(p - 1, p) * sum(a <= s for a in case.get(p, ())) for s in range(3))
        if max(thresholds) >= width + 1:
            raise ValueError(f"insufficient retained range at {p}")
        mean_threshold = 1 + delta*d + sum(
            F(p - 1, p) * F(6, 5*3**a) for a in case.get(p, ())
        )
        correction = sum(
            (t - j)*probability_list[j]
            for t, probability_list in zip(thresholds, probabilities)
            for j in range(1, int(t)+1)
            if probability_list[j]
        )
        loss = max(0, ceil((mean - S*mean_threshold + correction)/(d*(1-delta))))
        mass_debit += loss
        mean = ceil(mean*(1 + cap/(p-1)))
        moment = ceil(moment*(1 + cap*F(3*p-1,(p-1)**2)))
        factors = [0, ceil(S*(1-cap/p))]
        tail_probability = cap*(p-1)/p**2
        for _ in range(2, width+1):
            upper = ceil(S*tail_probability)
            factors.append(upper)
            if upper == 1:
                factors.extend([1]*(width+1-len(factors)))
                break
            tail_probability /= p
        next_state = [[0]*(width+1) for _ in range(3)]
        for s, states in enumerate(probabilities):
            for value, count in enumerate(states[1:], start=1):
                if count:
                    for j in range(1, width//value+1):
                        dest = min(j-1,2) if p == 3 else s
                        next_state[dest][value*j] += count*factors[j]
        probabilities = [[(a+S-1)//S for a in row] for row in next_state]
    return mass_debit, moment


def verify(tail_budget):
    expected_debits = (938773600338587870460242, 935988366744625136574530)
    expected_moment = 1908984849362909824839827392
    missing_labels = ((15, 21, 33), (15, 21, 45))
    for case, labels, expected_debit, margin in zip(
        CASES, missing_labels, expected_debits, (9, 12)
    ):
        debit, moment = evaluate(case)
        if (debit, moment) != (expected_debit, expected_moment):
            raise ValueError(f"incorrect prefix ledger for {labels}")
        residual = S * (S - debit) - tail_budget * moment
        if 1000 * residual <= margin * S**2:
            raise ValueError(f"insufficient survivor mass for {labels}")
        print(f"Missing moduli {labels} give live mass greater than {margin}/1000")
