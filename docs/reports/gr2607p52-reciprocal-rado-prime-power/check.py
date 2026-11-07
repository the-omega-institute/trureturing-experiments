"""Exact, independently written checks for Gaiser--Ramezanpour Problem 5.2."""
from bisect import bisect_left, bisect_right
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import lcm
from collections import Counter


def vp(n, p):
    assert n > 0 and p > 1
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def ceilq(x):
    return -(-x.numerator // x.denominator)


def monochromatic_solutions(k, colors, p=None):
    """Complete sorted-tuple search, including repetitions.

    Pruning uses exact necessary sum bounds and a residual denominator test.
    If exactly one available denominator has maximal positive p-valuation,
    split exhaustively by its multiplicity 0..j before searching the rest.
    This is a generic multiset enumeration optimization, not a distinctness rule.
    """
    stats = Counter()
    N = len(colors) - 1

    def dfs(j, ds, residual):
        stats['nodes'] += 1
        if j == 0:
            return () if residual == 0 else None
        if residual <= 0 or not ds:
            return None
        cap = residual - Q(j - 1, ds[-1])
        if cap <= 0:
            stats['size_cut'] += 1
            return None
        low = max(ds[0], ceilq(1 / cap))
        high = min(ds[-1], (j / residual).numerator // (j / residual).denominator)
        begin = bisect_left(ds, low)
        end = bisect_right(ds, high, begin)
        if begin == end:
            stats['size_cut'] += 1
            return None
        ds = ds[begin:]
        if p:
            heights = [vp(d,p) for d in ds]
            h = max(heights)
            if vp(residual.denominator,p) > h:
                stats['valuation_cut'] += 1
                return None
            peaks = [i for i,e in enumerate(heights) if e == h]
            if h > 0 and len(peaks) == 1:
                i = peaks[0]
                d = ds[i]
                other = ds[:i] + ds[i+1:]
                stats['multiplicity_splits'] += 1
                for t in range(j+1):
                    ans = dfs(j-t, other, residual-Q(t,d))
                    if ans is not None:
                        return tuple(sorted((d,)*t + ans))
                return None
        if j == 1:
            d = 1 / residual
            i = bisect_left(ds, d.numerator)
            if d.denominator == 1 and i < len(ds) and ds[i] == d.numerator:
                return (d.numerator,)
            return None
        for i,d in enumerate(ds):
            if d > high:
                break
            # Reuse index i, so equal consecutive denominators remain legal.
            ans = dfs(j-1, ds[i:], residual-Q(1,d))
            if ans is not None:
                return (d,) + ans
        return None

    for a in range(1,N+1):
        ds = tuple(d for d in range(a+1,N+1) if colors[d] == colors[a])
        if not ds:
            continue
        ans = dfs(k,ds,Q(1,a))
        if ans is not None:
            assert len(ans) == k and all(colors[b] == colors[a] for b in ans)
            assert sum((Q(1,b) for b in ans),Q(0)) == Q(1,a)
            yield a,ans,dict(stats)
    return dict(stats)


def avoiding(k, colors, p=None):
    g = monochromatic_solutions(k, colors, p)
    try:
        a, bs, st = next(g)
        raise AssertionError(('monochromatic', k, a, bs, st))
    except StopIteration as e:
        return e.value


def coloring(k, tail=2):
    M = 3 * k * k
    return [None] + [int(3 <= d < 3*k or d == M) for d in range(1, M + tail + 1)]


def main():
    # Cross-check pruning against direct enumeration, independently of the proof.
    cross = 0
    for N in range(3,13):
        for k in (2,3):
            for scheme in range(3):
                colors = [None]+[((d*d+scheme*d+scheme)%5)%2 for d in range(1,N+1)]
                brute = set()
                for a in range(1,N+1):
                    ds = [d for d in range(a+1,N+1) if colors[d] == colors[a]]
                    if any(sum((Q(1,b) for b in bs),Q(0)) == Q(1,a)
                           for bs in combinations_with_replacement(ds,k)):
                        brute.add(a)
                for prime in (None,3,5,7):
                    got = {a for a,_,_ in monochromatic_solutions(k,colors,prime)}
                    assert got == brute, (N,k,scheme,prime,brute,got)
                    cross += 1
    print('ENUMERATOR_CROSSCHECK direct brute-force comparisons:',cross)
    checked = 0
    for p in (3, 5, 7):
        for length in range(1, 5):
            for ds in combinations_with_replacement(range(1, 19), length):
                D = lcm(*ds)
                s = sum((Q(1, d) for d in ds), Q(0))
                h = max(vp(d, p) for d in ds)
                assert vp(D, p) == h
                assert D % s.denominator == 0
                assert vp(s.denominator, p) <= h
                checked += 1
    print('DENOMINATOR_LEMMA exact multisets checked:', checked)
    # Independent unpruned tuple enumeration checks source Remark 3.2 at k=3.
    R = set((1, 2)) | set(range(9, 27))
    cols = [None] + [int(d not in R) for d in range(1, 35)]
    tuples = 0
    for a in range(1, 35):
        ds = [d for d in range(a+1, 35) if cols[d] == cols[a]]
        for bs in combinations_with_replacement(ds, 3):
            tuples += 1
            assert sum((Q(1, b) for b in bs), Q(0)) != Q(1, a)
    print('FIDELITY Remark 3.2 k=3 N=34 avoids; tuples:', tuples)
    # Published Theorem 3.1 colouring, checked with p-adic cuts disabled at k=3,5.
    for p in (3, 5):
        st = avoiding(p, coloring(p, 0))
        print('FIDELITY Theorem 3.1 k=%d N=%d stats=%s' % (p, 3*p*p, st))
    # Negative control: the old unmodified block colouring fails at M.
    for k in (3, 5, 7):
        colors = [None] + [int(3 <= d < 3*k) for d in range(1, 3*k*k+1)]
        a, bs, _ = next(monochromatic_solutions(k, colors, k))
        assert a == 3*k and bs == (3*k*k,) * k
    print('NEGATIVE_CONTROL block colouring fails at M with repeated M, k=3,5,7')
    for p in (3, 5, 7):
        for m in (1, 2):
            k = p**m
            st = avoiding(k, coloring(k), p)
            print('CANDIDATE p=%d m=%d k=%d N=%d no monochromatic solution stats=%s' % (p,m,k,3*k*k+2,st))
    for k in range(3, 101):
        M = 3*k*k
        e1 = Q(1,M-2*k) + Q(k-1,M+2) - Q(k,M)
        e2 = Q(2,M-k) + Q(k-2,M+2) - Q(k,M)
        assert e1 == Q(10*k*k, M*(M-2*k)*(M+2)) > 0
        assert e2 == Q(14*k*k, M*(M-k)*(M+2)) > 0
        assert Q(k,M+2) > Q(1,3*k+1)
        assert [x for x in range(M-2*k+1,M+3) if x%k==0] == [M-k,M]
    print('TAIL_IDENTITIES_AND_MULTIPLICITIES exact k=3..100: OK')
    for p in (3,5,7,11,13):
        for m in range(1,5):
            k=p**m
            residual=Q(1,3*k)-Q(1,k*(3*k-1))
            assert residual == Q(3*k-4,3*k*(3*k-1))
            assert vp(residual.denominator,p)==m+vp(3,p)
            assert all(vp(d,p)<m for d in range(3*k*k-2*k+1,3*k*k+3) if d%k)
    print('RESIDUAL_VALUATIONS p=3,5,7,11,13 m=1..4: OK')
    # Sharpness for the nominated colouring: a red solution appears on adding M+3.
    # This is a boundary check of this colouring, not an upper bound on f_2(7).
    k=7; a,bs,_=next(monochromatic_solutions(k,coloring(k,3),7))
    print('BOUNDARY candidate colouring k=7 N=150 has solution:',a,bs)
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
