#!/usr/bin/env python3
"""Exact checks for arXiv:2610.32974v1, Conjecture 1.

The finite domain is every decreasing n in {0,1/8,...,1}^d, d=2,3,
and every admissible partition. For each input, ALL points of its continuous
feasible polytope are checked by subdivision into ordering chambers.
This is evidence for the proposed argument, not an arbitrary-dimension proof.
Only integers and Fraction arithmetic are used, including linear systems.
"""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))


def solve(rows, rhs):
    n = len(rhs)
    a = [list(map(Q, row)) + [Q(b)] for row, b in zip(rows, rhs)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return None
        a[j], a[pivot] = a[pivot], a[j]
        c = a[j][j]
        a[j] = [x/c for x in a[j]]
        for i in range(n):
            if i != j:
                c = a[i][j]
                a[i] = [x-c*y for x, y in zip(a[i], a[j])]
    return tuple(row[-1] for row in a)


def chamber_vertices(inequalities, forms):
    # inequalities: a.x <= b; forms: b_j(x)=a_j.x+c_j.
    # The bounded polytope is cut by EVERY b_i=b_j hyperplane.
    m = len(forms[0][0])
    lines = list(inequalities)
    for (a, c), (b, e) in combinations(forms, 2):
        lines.append((tuple(x-y for x, y in zip(a, b)), e-c))
    points = set()
    for rows in combinations(lines, m):
        x = solve([a for a, _ in rows], [b for _, b in rows])
        if x is not None and all(dot(a, x) <= b for a, b in inequalities):
            points.add(tuple(dot(a, x)+c for a, c in forms))
    assert points
    return points


def polytope(n, lam):
    total = sum(n)
    if len(n) == 2:
        # Schur-Horn: n2 <= b1 <= n1; b2=total-b1.
        return [((Q(1),), n[0]), ((Q(-1),), -n[1])], [
            ((Q(1),), Q(0)), ((Q(-1),), total)]
    forms = [((Q(1), Q(0)), Q(0)), ((Q(0), Q(1)), Q(0)),
             ((Q(-1), Q(-1)), total)]
    if lam == (2, 1):
        # Source Appendix on d=3: the COMPLETE Cauchy rectangle.
        return [((Q(1), Q(0)), n[0]), ((Q(-1), Q(0)), -n[1]),
                ((Q(0), Q(1)), n[1]), ((Q(0), Q(-1)), -n[2])], forms
    assert lam == (1, 1, 1)
    # In dimension three Schur-Horn is equivalent to n3 <= b_i <= n1,
    # together with the common total. These six halfspaces are exact.
    inequalities = []
    for a, c in forms:
        inequalities += [(a, n[0]-c), (tuple(-x for x in a), c-n[2])]
    return inequalities, forms


def optimum(n, lam):
    if lam == (2, 1):
        a = (n[0]+n[2])/2
        return tuple(sorted((a, n[1]), reverse=True)) + (a,)
    return (sum(n)/len(n),)*len(n)


def set_partitions(items):
    if not items:
        yield ()
        return
    x, *rest = items
    for parts in set_partitions(rest):
        yield ((x,),) + parts
        for j in range(len(parts)):
            yield parts[:j] + ((x,)+parts[j],) + parts[j+1:]


def subset_average_signatures(n):
    answer = set()
    for parts in set_partitions(list(range(len(n)))):
        v = []
        for part in parts:
            v += [sum(n[j] for j in part)/len(part)]*len(part)
        answer.add(tuple(sorted(v, reverse=True)))
    return answer


def check_majorized(x, y):
    assert sum(x) == sum(y)
    a, b = sorted(x, reverse=True), sorted(y, reverse=True)
    for r in range(1, len(x)+1):
        assert sum(a[:r]) <= sum(b[:r]), (x, y, r)


def fidelity():
    # Source VI.1: n=(1,0), lambda=(1,1), bmax=(1/2,1/2).
    # The real Hadamard unitary has entries +/-1/sqrt(2); its conjugated
    # rank-one matrix is A below. Exact idempotence certifies eigenvalues 1,0.
    a = ((Q(1, 2), Q(1, 2)), (Q(1, 2), Q(1, 2)))
    assert tuple(tuple(sum(a[i][k]*a[k][j] for k in range(2))
                       for j in range(2)) for i in range(2)) == a
    assert a[0][0]+a[1][1] == 1
    assert a[0][0]*a[1][1]-a[0][1]*a[1][0] == 0
    assert optimum((Q(1), Q(0)), (1, 1)) == (Q(1, 2), Q(1, 2))
    # Source d=3, (2,1): matrix [[a,0,c],[0,n2,0],[c,0,a]].
    # Its characteristic polynomial is (z-n1)(z-n2)(z-n3), since
    # a+c=n1 and a-c=n3; this avoids any numerical diagonalization.
    n = (Q(1), Q(1, 4), Q(0))
    av, off = (n[0]+n[2])/2, (n[0]-n[2])/2
    assert av+off == n[0] and av-off == n[2]
    assert optimum(n, (2, 1)) == (Q(1, 2), Q(1, 4), Q(1, 2))
    print('FIDELITY: source VI.1 and Table III / Appendix d=3 reproduced exactly')


def main():
    fidelity()
    instances = vertices = inequalities = homotopy = 0
    endpoint_inputs = repeated_inputs = 0
    for d in (2, 3):
        partitions = [(1, 1)] if d == 2 else [(2, 1), (1, 1, 1)]
        vectors = [tuple(reversed(a)) for a in combinations_with_replacement(
            [Q(j, 8) for j in range(9)], d)]
        for n in vectors:
            endpoint_inputs += int(0 in n or 1 in n)
            repeated_inputs += int(len(set(n)) < d)
            for lam in partitions:
                p, forms = polytope(n, lam)
                b = optimum(n, lam)
                # Membership in the source's complete feasible polytope.
                assert all(dot(a, b[:-1]) <= c for a, c in p)
                assert tuple(sorted(b, reverse=True)) in subset_average_signatures(n)
                points = chamber_vertices(p, forms)
                for c in points:
                    check_majorized(b, c)
                    inequalities += d
                    # Supporting-gradient certificates for f0=-x^2 and
                    # f1=-x^2-x^4. A linear certificate at every chamber
                    # vertex covers the entire feasible polytope.
                    # Their convex homotopy is covered for EVERY t in [0,1].
                    for grad in [tuple(-2*x for x in b),
                                 tuple(-2*x-4*x**3 for x in b)]:
                        assert dot(grad, tuple(y-x for x, y in zip(b, c))) <= 0
                        homotopy += 1
                instances += 1
                vertices += len(points)
    print('DOMAIN: d=2,3; every admissible partition; n_j in {0,1/8,...,1}')
    print('OCCUPATIONS: 210 inputs; endpoint inputs=%d; repeated inputs=%d' %
          (endpoint_inputs, repeated_inputs))
    print('CONTINUOUS_POLYTOPES: %d; ordering-chamber vertices=%d' % (instances, vertices))
    print('MAJORIZATION: %d exact partial-sum comparisons including equal totals' % inequalities)
    print('HOMOTOPY: %d exact supporting-gradient checks for both strictly concave endpoints' % homotopy)
    print('COVERAGE: chamberwise affine partial sums certify every feasible b for these inputs')
    print('LIMIT: finite occupation grid and dimensions 2,3; no general-dimensional proof or entropy numerics')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
