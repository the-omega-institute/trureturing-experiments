"""Exact cutoff-200 base spectra, Fourier simplex and optimal visible-to-hidden gaps."""

from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

import sympy as sp


def roots(matrix, variable):
    poly = matrix.charpoly(variable).as_poly()
    intervals = poly.intervals(eps=sp.Rational(1,10**15))
    assert sum(m for _,m in intervals)==matrix.rows
    print(sp.factor(poly.as_expr()))
    print("Distinct root multiplicities:", [m for _,m in intervals])
    return [interval for interval,_ in intervals]


def verify_joint_gap(adjacency, transformed):
    # In hidden adjacency coordinates, both operators use a common forbidden set.
    Q = sp.Rational
    forbidden = sorted(set(adjacency + transformed))
    assert all(left[1] < right[0] for left, right in zip(forbidden, forbidden[1:]))
    lower, upper = Q("-1.88"), Q("1.54")
    assert forbidden[0][1] < lower < upper < forbidden[-1][0]
    quarter = Q(1, 4)
    intervals = []
    endpoints = []
    for left, right in zip(forbidden, forbidden[1:]):
        lo = max(lower, left[0] + quarter)
        hi = min(upper, right[1] - quarter)
        if lo <= hi:
            intervals.append((lo, hi))
            endpoints.append((left, right))

    assert len(intervals) == 3
    ((minus_one, r), (zero, beta), (beta_again, _)) = endpoints
    assert minus_one == (-1, -1) and zero == (0, 0) and beta_again == beta
    assert r in adjacency and beta in transformed
    # Enclose the true feasible coordinates outward, so empty intersections
    # with the zero-sum plane exclude the corresponding gap assignments.
    possible = []
    for triple in combinations_with_replacement(range(3), 3):
        if sum(intervals[i][0] for i in triple) <= 0 <= sum(intervals[i][1] for i in triple):
            possible.append(triple)
    assert possible == [(0, 1, 1)]

    gap_lo = (r[0] + 2*beta[0])/3
    gap_hi = (r[1] + 2*beta[1])/3
    assert Q(".25029998053739") < gap_lo <= gap_hi < Q(".25029998053740")
    # At equality the negative coordinate is r-gap and both positive
    # coordinates are beta-gap. Test the opposite endpoints of their gaps.
    u_lo, u_hi = (beta[0]-r[1])/3, (beta[1]-r[0])/3
    assert u_lo > gap_hi and 1 - 2*u_hi > gap_hi
    assert Q("1.54") * u_hi < 1  # w_j = (1 - u*C_{negative row,j})/3
    print("Optimal joint distance in (0.25029998053739, 0.25029998053740).")


def verify():
    local = run_path(str(Path(__file__).with_name("ch39-local-product-contact.py")))
    sites, _ = local["lifted_sites"]()
    old = sorted(set(sites))
    distance = local["geometry"]["squared"]
    graph = sp.Matrix([[int(0 < distance(p,q) <= 200) for q in old] for p in old])
    sizes = Counter(sites)
    degree = list(graph * sp.ones(19,1))
    assert all(degree[i]==8 for i,p in enumerate(old) if sizes[p]==3)
    counts = Counter(sizes[p]*sizes[q] for i,p in enumerate(old)
                     for j,q in enumerate(old[:i]) if graph[i,j])
    assert counts == {1:30, 3:18, 9:3}
    last = Counter(sizes[p]*sizes[q] for i,p in enumerate(old)
                   for q in old[:i] if distance(p,q)==200)
    assert last == {1:3, 9:3}
    multiplicities_by_length = defaultdict(set)
    for i, p in enumerate(old):
        for q in old[:i]:
            length = distance(p, q)
            if 0 < length <= 200:
                multiplicities_by_length[length].add(sizes[p]*sizes[q])
    assert {length for length, kinds in multiplicities_by_length.items()
            if len(kinds) > 1} == {200}
    assert multiplicities_by_length[200] == {1, 9}

    t = sp.symbols("t")
    field = sp.QQ.cyclotomic_field(9)
    z = field.unit
    indices = (1,2,4)
    C = [[z**((k*j)%9)+z**((-k*j)%9) for j in indices] for k in indices]
    assert all(sum(row,field.zero)==field.zero for row in C)
    for i in range(3):
        for j in range(3):
            assert sum((C[k][i]*C[k][j] for k in range(3)),field.zero)==field.convert(9*int(i==j)-3)
    Q = sp.Rational
    cosine = sp.Poly(t**3-3*t+1,t).intervals(eps=Q(1,10**10))
    for ((lo,hi),mult),(a,b) in zip(cosine,[(-188,-187),(34,35),(153,154)]):
        assert mult==1 and Q(a,100)<lo<=hi<Q(b,100)

    adjacency = roots(graph,t)
    laplacian = roots(sp.diag(*degree)-graph,t)
    # A hidden coordinate lies in [-1.88,1.54]. Retain only gaps that
    # can leave it at least 0.359 from adjacent visible eigenvalues.
    retained = []
    for i,(left,right) in enumerate(zip(adjacency,adjacency[1:])):
        lo = max(Q(-188,100),left[0]+Q(359,1000))
        hi = min(Q(154,100),right[1]-Q(359,1000))
        if lo<=hi:
            retained.append((i,lo,hi))
    assert adjacency[0][1]<Q(-188,100) and adjacency[-1][0]>Q(154,100)
    assert [i for i,_,_ in retained]==[2,7]
    assert retained[0][2]<Q(-16,10) and retained[1][2]<Q(93,100)
    neg_lo,neg_hi=adjacency[7]
    assert Q("0.3595181495938") < -neg_hi <= -neg_lo < Q("0.3595181495939")
    assert min(abs(x) for i,r in enumerate(adjacency) if i!=7 for x in r)>-neg_lo
    assert 3-Q("1.5")*Q("1.6")+Q(".35")*Q(".93")-Q("1.87")*Q(".67")<0

    transformed=[(8-b,8-a) for a,b in reversed(laplacian)]
    near=[r for r in transformed if r[0]>=0 and r[0]<2]
    assert len(near)==4 and near[0]==(0,0) and near[-1][0]>Q("1.54")
    left,right = near[1]
    assert Q(".2948708841598") < left/2 <= right/2 < Q(".2948708841600")
    assert all(r[1]-l[0]<left for l,r in zip(near[1:],near[2:]))
    verify_joint_gap(adjacency, transformed)
    previous = max(b for a,b in transformed if b<0)
    assert previous < -3*right/2 and Q("1.54")*right/2<1
    print("Optimal adjacency distance in (0.3595181495938, 0.3595181495939).")
    print("Optimal Laplacian distance in (0.2948708841598, 0.2948708841600).")


if __name__ == "__main__":
    verify()
