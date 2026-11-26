"""Minimal uniform register and chord cutoff for the signed stationary three-time filter."""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations, product
from math import gcd

from sympy import Matrix, Rational

words = list(product(range(2), repeat=3))


def signed_law():
    p, q = Rational(4,13), Rational(9,13)
    gains = [(-1/q, 1/p), (-2/q, 1/p)]
    # Each row acts on the eight ordered word masses and the common free scale.
    equations = [[1]*8+[0], [w[0] for w in words]+[0]]
    rhs = [1,p]
    for pair in product(range(2), repeat=2):
        equations.append([int(w[:2]==pair)-int(w[1:]==pair) for w in words]+[0])
        rhs.append(0)
    for degree, coefficient in enumerate([18,-19,-3,4]):
        equations.append([gains[0][a]*gains[1][b] if a+b+c==degree else 0
                          for a,b,c in words]+[-Rational(coefficient,13)])
        rhs.append(0)
    matrix = Matrix(equations)
    independent = matrix.T.rref()[1]
    assert len(independent) == 9
    solution = matrix[list(independent),:].inv()*Matrix(rhs)[list(independent),:]
    expected = [p**sum(w)*q**(3-sum(w)) for w in words]+[1]
    assert list(solution)==expected and matrix*solution==Matrix(rhs)
    assert gcd(624,13**3)==13 and Q(4,13)**3==Q(64,2197)


def distance(p, q):
    return sum((a-b)**2 for a,b in zip(p,q,strict=True))


def groups():
    outer = sorted(set(permutations((1,5,29))))
    inner = sorted(set(permutations((7,17,23))))
    old = sorted(set(permutations((11,11,25))))
    lifted = sorted(set(permutations((13,13,23))))
    assert all(min(distance(p,q) for q in inner+old+lifted)==152 for p in outer)
    assert Q(216,624)>Q(9,13)**3
    # Each four-type site divides its internal channels equally between two six-type neighbors.
    result = []
    for four, six, count, cutoff in ((old,outer,16,152),(lifted,inner,48,52)):
        used = []
        for p in four:
            neighbors = [q for q in six if distance(p,q)==cutoff]
            assert len(neighbors)==2
            for i,q in enumerate(neighbors):
                used.append(q)
                four_channels = [(p,j) for j in range(i*count//2,(i+1)*count//2)]
                result.append(([(q,j) for j in range(36)],four_channels))
        assert sorted(used)==six
    assert len({ch for group in result for side in group for ch in side})==624
    return result


def memoryless_cutoff():
    four = sorted((p,k) for t,h in (((11,11,25),1),((13,13,23),3))
                  for p in set(permutations(t)) for k in range(h))
    six = sorted(p for t in ((1,5,29),(7,17,23)) for p in permutations(t))
    lifted = [i for i,(p,_) in enumerate(four) if sorted(p)==[13,13,23]]
    near = {j for i in lifted for j,q in enumerate(six) if distance(four[i][0],q)<244}
    assert len(lifted)==9 and len(near)==6
    match = {}
    def augment(i, seen):
        for j,q in enumerate(six):
            if j not in seen and distance(four[i][0],q)<=244:
                seen.add(j)
                if j not in match or augment(match[j], seen):
                    match[j] = i
                    return True
        return False
    assert all(augment(i,set()) for i in range(12))
    assert len(match)==12 and len(set(match.values()))==12


def reversible_witness():
    counts = Counter()
    for zeros, ones in groups():
        # Number per individual ordered word with 0, 1, 2, 3 ones, respectively.
        profile = [3380,1352,0,0] if len(ones)==8 else [2452,1240,1152,512]
        pools = [[(ch,seed) for ch in channels for seed in range(169)] for channels in (zeros,ones)]
        succ = {}
        color = {state:i for i,pool in enumerate(pools) for state in pool}
        for amount, word in zip(profile, [(0,), (0,0,1), (0,1,1), (1,)], strict=True):
            for _ in range(amount):
                cycle = [pools[bit].pop() for bit in word]
                succ.update(zip(cycle,cycle[1:]+cycle[:1],strict=True))
        assert not any(pools) and set(succ)==set(succ.values())
        for a,b in succ.items():
            assert distance(a[0][0],b[0][0])<=152
            counts[color[a],color[b],color[succ[b]]] += 1
    total = 624*169
    assert sum(counts.values())==total
    for w in words:
        assert Q(counts[w],total)==Q(4,13)**sum(w)*Q(9,13)**(3-sum(w))
    print('Eight stationary word counts:', [counts[w] for w in words])


if __name__ == '__main__':
    signed_law()
    memoryless_cutoff()
    reversible_witness()
