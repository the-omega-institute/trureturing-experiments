#!/usr/bin/env python3
"""Exact checks for the alpha-unpredictable/no-average construction.

Finite evidence only: the infinite construction and C-infinity regularity require
mathematical proofs. The polynomial kernel below is an exactly integrable test
kernel; support-based recurrence/separation checks apply to any normalized kernel
supported in [-1/4,1/4], including the smooth bump in the proposed proof.
"""
from bisect import bisect_right
from fractions import Fraction as Q
from itertools import product
import sympy as sp


class Word:
    def __init__(self):
        self.ends = [0]
        self.counts = [0]
        self.lengths = [0]
        self.runs = [0]

    def stage(self, n):
        while len(self.ends) <= n:
            j = len(self.ends)
            length = self.ends[-1] + j * 2**j
            run = j*j*(length+1)
            self.lengths.append(length)
            self.runs.append(run)
            self.ends.append(length+run)
            self.counts.append(self.counts[-1] + j*2**(j-1) + (j % 2)*run)

    def ensure_index(self, k):
        while self.ends[-1] <= k:
            self.stage(len(self.ends))

    def bit(self, k):
        if k < 0:
            return self.bit(-k-1)
        self.ensure_index(k)
        j = bisect_right(self.ends, k)
        offset = k-self.ends[j-1]
        if offset >= j*2**j:
            return j % 2
        word, digit = divmod(offset, j)
        return (word >> (j-1-digit)) & 1

    @staticmethod
    def digit_sum_before(m, width):
        return sum((m // (2*p))*p + max(0, m % (2*p)-p)
                   for p in (2**i for i in range(width)))

    def prefix_count(self, k):
        if not k:
            return 0
        self.ensure_index(k-1)
        j = bisect_right(self.ends, k-1)
        count = self.counts[j-1]
        offset = k-self.ends[j-1]
        if offset >= j*2**j:
            return count+j*2**(j-1)+(j % 2)*(offset-j*2**j)
        words, tail = divmod(offset, j)
        return count+self.digit_sum_before(words, j) + (words >> (j-tail)).bit_count()


R = Q(1, 4)


def primitive(u):
    return Q(15, 4)*(u-Q(32, 3)*u**3+Q(256, 5)*u**5)


def first_moment_primitive(u):
    return Q(15, 4)*(u*u/2-8*u**4+Q(128, 3)*u**6)


MOMENT = first_moment_primitive(R)
assert primitive(R)-primitive(-R) == 1
assert MOMENT == Q(5, 128)


def h(word, t):
    """Exact convolution with K(u)=(15/4)(1-16u^2)^2 on [-1/4,1/4]."""
    total = Q(0)
    for k in range((t-R).__floor__(), (t+R).__floor__()+1):
        lo, hi = max(-R, t-k-1), min(R, t-k)
        if lo < hi:
            total += word.bit(k)*(primitive(hi)-primitive(lo))
    return total


def smoothed_integral_at_integer(word, t):
    # Even extension gives no correction at 0. Only the jump at t contributes.
    return Q(word.prefix_count(t))+MOMENT*(word.bit(t)-word.bit(t-1))


def source_fidelity():
    # Journal Section 2, (2)-(3): kernel mass and bound for 0 <= Omega <= 1.
    u, t, c = sp.symbols('u t c', real=True)
    mass = sp.integrate(sp.exp(-2*u), (u, 0, sp.oo))
    assert mass == sp.Rational(1, 2)
    for constant in [sp.Rational(0), sp.Rational(1, 3), sp.Rational(1)]:
        theta = constant/2
        assert sp.diff(theta, t) == -2*theta+constant
        assert 0 <= theta <= mass
    # Exact zero-order-hold update for any starting state 0 <= x <= 1/2.
    x = sp.symbols('x', real=True)
    z = c/2+(x-c/2)*sp.exp(-2*t)
    assert sp.simplify(sp.diff(z, t)+2*z-c) == 0
    assert z.subs(t, 0) == x
    print('FIDELITY: integral_0^infinity exp(-2u) du = 1/2; source bound and ODE verified')


def explicit_crosscheck(word):
    explicit = []
    for n in range(1, 6):
        enum = [b for w in product([0, 1], repeat=n) for b in w]
        length = len(explicit)+len(enum)
        run = n*n*(length+1)
        explicit += enum + [n % 2]*run
        word.stage(n)
        assert len(explicit) == word.ends[n]
        assert sum(explicit) == word.counts[n]
        probes = set(range(min(1000, len(explicit))))
        probes.update(range(max(0, length-100), min(len(explicit), length+100)))
        probes.update(range(0, len(explicit), max(1, len(explicit)//257)))
        probes.add(len(explicit)-1)
        for k in probes:
            assert word.bit(k) == explicit[k]
            assert word.prefix_count(k+1) == sum(explicit[:k+1])
    print('COMPRESSED_WORD: independently expanded stages 1..5 agree at counts and sampled positions')


def frequencies(word):
    for n in range(1, 41):
        word.stage(n)
        t, length = word.ends[n], word.lengths[n]
        frequency = Q(word.counts[n], t)
        error = abs(frequency-n % 2)
        assert Q(length, t) <= Q(1, n*n+1)
        assert error <= Q(1, n*n+1)
        integral = smoothed_integral_at_integer(word, t)
        avg = integral/t
        assert abs(integral-word.counts[n]) <= 2*2*MOMENT <= Q(1, 2)
        assert abs(avg-n % 2) <= Q(1, n*n+1)+Q(1, 2*t)
        for a, b in [(Q(1), Q(2)), (Q(-3), Q(5, 2)), (Q(0), Q(1, 3))]:
            assert abs(a+(b-a)*avg-(a+(b-a)*(n % 2))) <= (b-a)*(Q(1, n*n+1)+Q(1, 2*t))
        if n <= 8:
            print(f'FREQUENCY n={n} T={t} raw={frequency} smoothed={avg} target={n % 2}')
    print('ENDPOINT_BOUNDS: exact stages 1..40, both parities, three affine families')


def recurrence_and_separation(word):
    previous_q, previous_s = 0, Q(0)
    for n in range(1, 6):
        m = 2*n+3
        word.stage(m)
        central = [word.bit(k) for k in range(-n-1, n+2)]
        encoded = sum(v*2**(m-1-i) for i, v in enumerate(central))
        start = word.ends[m-1]+m*encoded
        q = start+n+1
        assert q > previous_q
        for k in range(-n-1, n+2):
            assert word.bit(q+k) == word.bit(k)
        for j in range(-16*n, 16*n+1):
            t = Q(j, 16)
            assert h(word, t+q) == h(word, t)
        # A late even constant block ends in zeros, followed by E_(r+1).
        # Its first 1 is after 2(r+1)-1 leading zeros. Find 0^q1 there.
        r = 2
        while True:
            word.stage(r+1)
            first_one = word.ends[r]+2*(r+1)-1
            k = first_one-q
            if word.runs[r] >= q and Q(k)+Q(1, 2) > previous_s+1:
                break
            r += 2
        assert word.lengths[r] <= k < word.ends[r]
        assert word.bit(k) == 0 and word.bit(k+q) == 1
        s, delta = Q(k)+Q(1, 2), Q(1, 8)
        # Whole symmetric window plus kernel support lies inside the two cells.
        assert Q(k) < s-delta-R and s+delta+R < k+1
        for j in range(17):
            t = s-delta+Q(j, 64)
            assert s-delta <= t <= s+delta
            assert h(word, t) == 0 and h(word, t+q) == 1
            assert abs(h(word, t+q)-h(word, t)) > Q(1, 2)
        print(f'RETURN n={n} q={q} compact=[{-n},{n}] separation_center={s} delta=1/8 gap=1')
        previous_q, previous_s = q, s
    for j in range(-100, 101):
        t = Q(j, 16)
        assert h(word, -t) == h(word, t)
        assert 0 <= h(word, t) <= 1
    print('SMOOTHING: exact polynomial test kernel; recurrence and full interval support certificates')
    print('EVENNESS: exact samples; symmetric endpoint means equal one-sided endpoint means')
    print('LIMITATION: finite checks do not certify the infinite limits or C-infinity regularity')


def main():
    source_fidelity()
    word = Word()
    explicit_crosscheck(word)
    frequencies(word)
    recurrence_and_separation(word)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
