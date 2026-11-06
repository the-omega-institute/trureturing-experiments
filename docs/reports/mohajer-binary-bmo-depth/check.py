#!/usr/bin/env python3
"""Exact checks for the finite binary-tree BMO depth-three argument.

Source fidelity: Mohajer, arXiv:2610.14459v1, explicit depth-four
counterexample (Section 4.4). All arithmetic uses Fraction or integers.
The symbolic coefficient tests certify the six universal linear bounds;
the finite grids are supplementary tests, not a proof over the reals.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from random import Random


def nodes(m):
    return [tuple(range(j * (2 ** (m-d)), (j+1) * (2 ** (m-d))))
            for d in range(m+1) for j in range(2**d)]


def mean(x, I):
    return sum((x[i] for i in I), Q(0)) / len(I)


def deviation(x, I):
    a = mean(x, I)
    return sum((abs(x[i]-a) for i in I), Q(0)) / len(I)


def norm(x):
    m = len(x).bit_length()-1
    assert len(x) == 2**m
    return max(deviation(x, I) for I in nodes(m))


def linear_coeff(I, U, n=8):
    I, U = set(I), set(U)
    assert U <= I and I
    return tuple(Q(2, len(I)) * ((i in U) - Q(len(U), len(I)))
                 if i in I else Q(0) for i in range(n))


def expand(terms):
    return tuple(sum((w * linear_coeff(I, U)[i] for w, I, U in terms), Q(0))
                 for i in range(8))


def auts(m):
    if m == 0:
        return [(0,)]
    h = 2**(m-1)
    return [tuple((b+h if flip else b) for b in left) +
                  tuple((b if flip else b+h) for b in right)
            for left in auts(m-1) for right in auts(m-1) for flip in (False, True)]


def tail_check(x):
    n = len(x)
    if mean(x, range(n)) < 0:
        x = tuple(-v for v in x)
    mu = mean(x, range(n))
    y = tuple(abs(v) for v in x)
    t = mean(y, range(n))
    S = {i for i in range(n) if y[i] > t}
    Sc = set(range(n))-S
    p = Q(len(S), n)
    N = sum((max(-x[i], Q(0)) for i in S), Q(0))/n
    R = sum((max(-x[i], Q(0)) for i in Sc), Q(0))/n
    D = sum((max(mu-x[i], Q(0)) for i in Sc), Q(0))/n
    U = sum((max(x[i]-mu, Q(0)) for i in Sc), Q(0))/n
    r = Q(sum(x[i] < 0 for i in S), n)
    dx, dy = deviation(x, range(n)), deviation(y, range(n))
    assert dx/2 == N+r*mu+D
    assert dy/2 == 2*(1-p)*N-2*p*R+D-U
    assert dy <= max(Q(1), 2*(1-p))*dx
    if not S or all(x[i] >= 0 for i in S) or all(x[i] <= 0 for i in S):
        assert dy <= dx
    T = sum((y[i] for i in S), Q(0))
    C = sum((y[i] for i in Sc), Q(0))
    assert dy == Q(2, n*n)*((n-len(S))*T-len(S)*C)
    return len(S), len({1 if x[i] > 0 else -1 for i in S})


def main():
    a = tuple(map(Q, (0,-8,-12,-12,0,0,0,0,6,1,1,1,1,10,14,14)))
    expected = [[Q(5)], [Q(5),Q(5)], [Q(4),Q(0),Q(15,8),Q(35,8)],
                [Q(4),Q(0),Q(0),Q(0),Q(5,2),Q(0),Q(9,2),Q(0)], [Q(0)]*16]
    actual = [[deviation(a, I) for I in nodes(4)[2**d-1:2**(d+1)-1]] for d in range(5)]
    assert actual == expected
    assert mean(a, range(16)) == 1
    assert norm(a) == 5 and norm(tuple(map(abs, a))) == Q(41,8)
    for d, row in enumerate(actual):
        print('FIDELITY input oscillations level', d, ':', ','.join(map(str,row)))
    print('FIDELITY depth-four input=5 output=41/8 ratio=41/40')

    # Exact symbolic forms; strings encode digit sets, with R expanded explicitly.
    R = '01234567'
    certs = {
        (1,2): ([(12,'01','1'),(8,'23','2'),(4,'0123','23')], (-7,5,5,-3,0,0,0,0)),
        (1,4): ([(10,'01','1'),(4,'45','4'),(8,'4567','4')], (-5,5,0,0,5,-3,-1,-1)),
        (2,3): ([(4,'01','1'),(4,'23','3'),(8,'0123','2'),(16,'0123','23')], (-7,-3,5,5,0,0,0,0)),
        (2,4): ([(8,'01','1'),(4,'23','2'),(8,'45','4'),(8,'0123','2'),(4,'4567','45')], (-5,3,5,-3,5,-3,-1,-1)),
        (4,5): ([(8,'01','1'),(8,'4567','45'),(16,R,'45')], (-5,3,-1,-1,5,5,-3,-3)),
        (4,6): ([(8,'01','1'),(4,'45','4'),(4,'67','6'),(16,R,'46')], (-5,3,-1,-1,5,-3,5,-3)),
    }
    sign_cones = 0
    for pair, (encoded, expected_coeff) in certs.items():
        terms = [(w, tuple(map(int,I)),tuple(map(int,U))) for w,I,U in encoded]
        coeff = expand(terms)
        assert coeff == expected_coeff
        weight = sum(w for w,_,_ in terms)
        assert all(w >= 0 for w,_,_ in terms) and weight <= 32
        assert coeff[0] <= -5 and all(coeff[i] == 5 for i in pair)
        other = sorted(set(range(8))-{0,*pair})
        assert all(abs(coeff[i]) <= 3 for i in other)
        # Write x_0=-z_0, x_i=z_i at selected positives, and x_l=s_l*z_l
        # elsewhere. H-Q must have nonnegative coefficients in z_l >= 0.
        for signs in product((-1,1), repeat=len(other)):
            margin = [-coeff[0]-5] + [coeff[i]-5 for i in pair]
            margin += [coeff[i]*s+3 for i,s in zip(other,signs)]
            assert all(v >= 0 for v in margin)
            sign_cones += 1
        print('CERTIFICATE pair', pair, 'coefficients', ','.join(map(str,coeff)), 'weight', weight)
    print('SYMBOLIC domination sign cones:', sign_cones)

    G = auts(3)
    node_sets = {frozenset(I) for I in nodes(3)}
    assert len(G) == len(set(G)) == 128
    for g in G:
        assert sorted(g) == list(range(8))
        assert {frozenset(g[i] for i in I) for I in nodes(3)} == node_sets
    G0 = [g for g in G if g[0] == 0]
    assert len(G0) == 16
    orbits = [{tuple(sorted((g[i],g[j]))) for g in G0} for i,j in certs]
    assert set.union(*orbits) == set(combinations(range(1,8),2))
    assert sum(map(len,orbits)) == 21
    # Cover all 168 signed-tail placements, not just the fixed-leaf presentation.
    all_marked = {(p, tuple(pair)) for p in range(8)
                  for pair in combinations([i for i in range(8) if i != p],2)}
    covered = {(g[0],tuple(sorted((g[i],g[j])))) for g in G for i,j in certs}
    assert covered == all_marked
    print('ORBITS automorphisms=128 stabilizer=16 sizes='+','.join(map(str,map(len,orbits)))+' marked_tails=168')

    # Size-two mixed tails: opposite halves have a universal coefficient bound.
    for p in range(4):
        for n in range(4,8):
            coeff = expand([(8,tuple(range(4)),(p,)),
                            (8,tuple(range(4,8)),tuple(i for i in range(4,8) if i != n))])
            assert coeff[p] == 3 and coeff[n] == -3
            assert all(abs(coeff[i]) <= 1 for i in range(8) if i not in (p,n))
    print('SIZE_TWO opposite-half certificates=16; same-half estimate is the triangle inequality')

    for m in (1,2,3):
        x = (Q(0),)*(2**(m-1))+(Q(2),)*(2**(m-1))
        assert norm(x) == norm(tuple(map(abs,x))) == 1
    print('SHARPNESS depths=1,2,3 input=output=1')

    tested, hard, ties = 0, {2:0,3:0}, 0
    cases = [tuple(map(Q,x)) for x in product((-1,0,1), repeat=8)]
    rng = Random(261014459)
    cases += [tuple(Q(rng.randrange(-16,17),rng.randrange(1,6)) for _ in range(8))
              for _ in range(4096)]
    for x in cases:
        M = norm(x)
        assert norm(tuple(map(abs,x))) <= M
        k, signs = tail_check(x)
        if signs == 2 and k in hard:
            hard[k] += 1
        t = mean(tuple(map(abs,x)),range(8))
        ties += any(abs(v) == t for v in x)
        for I in nodes(3):
            assert deviation(tuple(map(abs,x)),I) <= M
            # Check subset-deviation bounds for each node, all subsets.
            center = mean(x,I)
            for U in (tuple(i for i in I if x[i] > center), I, ()):
                L = Q(2,len(I))*(sum((x[i] for i in U),Q(0))-Q(len(U),len(I))*sum((x[i] for i in I),Q(0)))
                assert L <= deviation(x,I) <= M
        tested += 1
    print('EXACT_SAMPLES ternary_grid=6561 rational_samples=4096 total='+str(tested))
    print('EXACT_SAMPLES mixed_tail_2='+str(hard[2])+' mixed_tail_3='+str(hard[3])+' root_ties='+str(ties))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
