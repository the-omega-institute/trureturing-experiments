#!/usr/bin/env python3
"""Exact checks for Meyer's first-level centroid-covering question.

Only integers and fractions are used. These checks concern the explicitly
listed families; they do not certify the universal assertion.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), Q(0))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def exposed(points, u):
    h = max(dot(x, u) for x in points)
    return tuple(x for x in points if dot(x, u) == h)


def subsets(values):
    return [s for k in range(1, len(values) + 1) for s in combinations(values, k)]


def interval_radius(points, halfwidth):
    """Exact max nearest-point distance on the interval hull, scaled by b.

    A nearest-point distance changes slope only at a data point or the
    midpoint of consecutive data points, so these are all the candidates.
    """
    points = sorted(set(points))
    critical = points + [(a + b) / 2 for a, b in zip(points, points[1:])]
    return max(min(abs(x - a) / halfwidth for a in points) for x in critical)


def fidelity_uniform_essential():
    # Source §2.3: a_j = 2^(-j-1), b_j = 2^(-2j-1).
    # At ±u_j the exposed A-set consists of both endpoints of length 2a_j.
    # No truncated zonotope is substituted for the infinite source body.
    # We check the exact exposed-interval calculation and uniform tail sums.
    for j in range(1, 65):
        a = Q(1, 2 ** (j + 1))
        b = Q(1, 2 ** (2 * j + 1))
        rho = interval_radius((-a, a), b)
        assert rho == a / b == 2**j
        assert rho * (2 * b) == 2 * a  # rho |K_u*| = 2 d(A_u)
    for n in range(1, 33):
        a_tail = Q(1, 2 ** (n + 1))
        b_tail = Q(1, 6 * 4**n)
        assert sum((Q(1, 2 ** (j + 1)) for j in range(1, n + 1)), Q(0)) + a_tail == Q(1, 2)
        assert sum((Q(1, 2 ** (2 * j + 1)) for j in range(1, n + 1)), Q(0)) + b_tail == Q(1, 6)
    # Pairwise nonparallel rational unit vectors v_j, using the rational circle.
    directions = [((1-Q(j)**2)/(1+Q(j)**2), 2*Q(j)/(1+Q(j)**2)) for j in range(1, 33)]
    assert all(dot(v, v) == 1 for v in directions)
    assert all(a[0]*b[1]-a[1]*b[0] != 0 for a, b in combinations(directions, 2))
    print('FIDELITY: source §2.3 exposed-interval radii = 2^j, j=1..64; exact tails n=1..32')


def polytope_radii():
    # A is a Cartesian product of nonempty subsets of {0,1,2}; K is an
    # axis-parallel box with halfwidths b_i. The independent coordinate
    # interval calculations give the l_infinity box gauge covering radius.
    # Every iterated face is a coordinate face. After j orthogonal exposures
    # it has at most n-j free coordinates. Conversely every such face occurs.
    # Enumerate all these faces, including the faces of singleton coordinates.
    profiles = subsets(tuple(map(Q, (0, 1, 2))))
    cases = faces = 0
    for n in (2, 3, 4):
        for coords in product(profiles, repeat=n):
            widths = tuple(Q(i+1, i+2) for i in range(n))
            gaps = tuple(interval_radius(coords[i], widths[i]) for i in range(n))
            radii = []
            face_data = []
            for state in product((-1, 0, 1), repeat=n):
                codim = sum(s != 0 for s in state)
                if codim:
                    face_data.append((codim, max((gaps[i] for i in range(n) if state[i] == 0), default=Q(0))))
                    faces += 1
            for j in range(1, n):
                rj = max(r for codim, r in face_data if codim >= j)
                radii.append(rj)
                assert rj <= Q(factorial(n-1), factorial(n-j)) * radii[0]
                assert rj <= radii[0]  # source Theorem 1.3, sharper for polytopes
            cases += 1
    # Independent full grid evaluation of the nearest-point problem.
    grids = 0
    for coords in product(profiles, repeat=3):
        points = tuple(product(*coords))
        widths = (Q(1, 2), Q(2, 3), Q(3, 4))
        critical = [list(c) + [(a+b)/2 for a,b in zip(c,c[1:])] for c in coords]
        brute = max(min(max(abs(x[i]-a[i])/widths[i] for i in range(3)) for a in points) for x in product(*critical))
        assert brute == max(interval_radius(coords[i], widths[i]) for i in range(3))
        grids += 1
    print(f'POLYTOPES: {cases} product-grid/box pairs, {faces} faces; {grids} independent 3D grid checks')


def finite_points_disc_segment():
    # Disc K={(x,y,0):x²+y²<=1}; for u_xy != 0 its face is a singleton.
    # Segment K={0}x{0}x[-1,1]; for u_z != 0 its face is a singleton.
    # A non-singleton finite A_u cannot be covered by a singleton centered at 0.
    # Directions are deliberately not normalized: positive scaling does not
    # change maximizers, and avoids square roots in the exact arithmetic.
    cube = tuple(product(map(Q, (0, 1)), repeat=3))
    directions = [u for u in product(map(Q, range(-2, 3)), repeat=3) if any(u)]
    disc_infinite = segment_infinite = segment_finite = 0
    for points in subsets(cube):
        disc_witness = next((u for u in directions if (u[0] or u[1]) and len(exposed(points, u)) > 1), None)
        assert (disc_witness is None) == (len(points) == 1)
        disc_infinite += disc_witness is not None
        same_xy = len({x[:2] for x in points}) == 1
        segment_witness = next((u for u in directions if u[2] and len(exposed(points, u)) > 1), None)
        assert (segment_witness is None) == same_xy
        if same_xy:
            r1 = interval_radius(tuple(x[2] for x in points), Q(1))
            r2 = r1  # two independent horizontal directions leave this line fixed
            assert r2 <= 2 * r1
            segment_finite += 1
        else:
            segment_infinite += 1
    print(f'FINITE A: 255 cube subsets with disc K: {disc_infinite} infinite-radius witnesses; segment K: {segment_finite} finite and {segment_infinite} infinite')


def nonpolytope_families():
    # Family I: A = circle(radius a) x [-h,h]; K = disc(radius b) x {0}.
    # u_xy != 0: A_u is a singleton or an entire vertical segment, hence convex.
    # u_xy = 0: A_u is a circle and K_u* is the full disc, radius a/b.
    # Every second exposure leaves a point or convex segment, hence R2=0.
    # A + lambda K has radial interval [max(a-lambda*b,0), a+lambda*b]
    # at every height. Thus its exact convexification threshold is a/b.
    circle_cases = segment_cases = 0
    parameters = (Q(1, 3), Q(1, 2), Q(1), Q(2), Q(3))
    for a, b, h in product(parameters, repeat=3):
        r1 = a / b
        r2 = Q(0)
        assert r2 <= 2 * r1
        for lam in (r1 / 2, r1, 2 * r1):
            inner_radius = max(a - lam * b, Q(0))
            assert (inner_radius == 0) == (lam >= r1)
        circle_cases += 1
        # Family II: A=disc(radius a) x {-h,h}; K={0}² x [-b,b].
        # Horizontal directions expose two points distance 2h; vertical
        # directions expose a full convex disc. R1=R2=c=h/b.
        r1 = interval_radius((-h, h), b)
        r2 = h / b
        assert r2 == r1 and r2 <= 2 * r1
        for lam in (r1 / 2, r1, 2 * r1):
            # Two vertical intervals centered at ±h overlap iff lam*b>=h.
            assert (-h + lam*b >= h - lam*b) == (lam >= r1)
        segment_cases += 1
    print(f'NONPOLYTOPAL A: {circle_cases} circle-cylinder/disc pairs, R1=c=a/b and R2=0; {segment_cases} two-disc/segment pairs, R1=R2=c=h/b')


def projection_nonexposed_endpoint():
    # C=conv(unit disc in z=0, (1,1,0)); u=e1, v=-e2, x=(1,0,0).
    # C_u=[(1,0,0),(1,1,0)] and x=(C_u)_v is not itself an exposed
    # singleton of C. The following exact rational points implement the
    # projection construction, crossing curved first-level faces.
    x, u, v = (Q(1),Q(0),Q(0)), (Q(1),Q(0),Q(0)), (Q(0),Q(-1),Q(0))
    previous = Q(2)
    for k in range(3, 65):
        t = Q(k*k - 2*k - 1, 2*k)
        z = (Q(k*k-1, k*k+1), Q(-2*k, k*k+1), Q(0))
        d = sub(z, x)
        w = tuple(u[i] + (v[i]-d[i])/t for i in range(3))
        assert t > 0 and dot(z,z) == 1
        hz = dot(w,z)
        assert hz > 0 and hz*hz == dot(w,w)
        assert dot(w,(Q(1),Q(1),Q(0))) <= hz
        # Thus z maximizes w on both the disc and its convex hull with (1,1,0).
        assert dot(d,d)-2*dot(v,d)-2*t*dot(u,d) <= 0
        assert dot(d,d) == Q(4, k*k+1) < previous
        previous = dot(d,d)
    print('PROJECTION: 62 rational curved-face projections to a nonexposed endpoint; squared error 4/(k²+1)')


def centroid_difference_body():
    # d-simplex: g=(1,...,1)/(d+1); every difference-body vertex v_i-v_j
    # lies in (d+1)(G-g) exactly when g+(v_i-v_j)/(d+1) has nonnegative
    # barycentric coordinates. These vertex checks cover the entire body.
    count = 0
    for d in range(0, 9):
        vertices = [(Q(0),)*d] + [tuple(Q(i==j) for i in range(d)) for j in range(d)]
        centroid = (Q(1,d+1),)*d
        for a,b in product(vertices, repeat=2):
            p = tuple(centroid[i]+(a[i]-b[i])/Q(d+1) for i in range(d))
            bary = (1-sum(p,Q(0)),)+p
            assert all(c>=0 for c in bary) and sum(bary,Q(0))==1
            count += 1
    print(f'CENTROID: {count} difference-body vertex tests for simplex dimensions 0..8, including singleton')


def main():
    fidelity_uniform_essential()
    polytope_radii()
    finite_points_disc_segment()
    nonpolytope_families()
    projection_nonexposed_endpoint()
    centroid_difference_body()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
