"""Signed Regge product families. Angles are enclosed in dyadic whole turns."""
from collections import Counter
from functools import cache
from itertools import combinations, permutations, product
from math import factorial

from flint import arb, ctx, fmpq, fmpq_mat

ctx.prec = 192
GRID = 1 << 42
SIDE = 4


@cache
def cosine_tick(k):
    """Dyadic rotations using the half-angle identity and outward square roots."""
    if k == 0:
        return arb(1)
    if k == GRID // 2:
        return arb(-1)
    if k == GRID // 4:
        return arb(0)
    if k > GRID // 4:
        return -cosine_tick(GRID // 2 - k)
    return ((1 + cosine_tick(2*k))/2).sqrt()


def angle_turns(cosine):
    lo, hi = 0, GRID // 2
    assert cosine > -1 and cosine < 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        candidate = cosine_tick(mid)
        if cosine < candidate:
            lo = mid
        elif cosine > candidate:
            hi = mid
        else:
            assert cosine < cosine_tick(mid-1) and cosine > cosine_tick(mid+1)
            lo, hi = mid-1, mid+1
            break
    return arb(fmpq(lo+hi, 2*GRID), fmpq(hi-lo, 2*GRID))


def squared_length(mode, parity, mask, t):
    d = (mask & 7).bit_count()
    base = fmpq(d*(4-d), 4)
    if mode == "parity":
        base *= (1 + fmpq(parity,20)) * (1 + fmpq(parity*(-1)**d,20))
    elif d == 3:
        base += fmpq(1,20)
    assert 0 <= base < 2
    return base/(t*t) + (t**6 if mask & 8 else 0)


def gram(distances, vertices):
    start, *ends = vertices
    return fmpq_mat([[(distances[start,i]+distances[start,j]-distances[i,j])/2
                      for j in ends] for i in ends])


@cache
def local_geometry(mode, parity, order, t):
    n = len(order)
    masks = tuple(sum(1 << i for i in order[:j]) for j in range(n+1))
    idx = tuple(range(n+1))
    distances = {(i,i): fmpq(0) for i in idx}
    for i,j in combinations(idx,2):
        sign = parity * (-1)**(masks[i] & 7).bit_count()
        q = squared_length(mode, sign, masks[j]-masks[i], t)
        distances[i,j] = distances[j,i] = q
    g = gram(distances, idx)
    if n == 3:
        # Consecutive coordinate steps, before the cumulative vertex basis change.
        diff = fmpq_mat([[int(i==j)-int(i==j+1) for j in range(n)] for i in range(n)])
        steps = diff*g*diff.transpose()
        for k in range(1,4):
            minor = fmpq_mat([[steps[i,j]-fmpq(int(i==j),8) for j in range(k)]
                              for i in range(k)])
            assert minor.det() > 0
    assert g.det() > 0
    d = fmpq_mat([[-1]*n] + [[int(i==j) for j in range(n)] for i in range(n)])
    normals = d*g.inv()*d.transpose()
    hinges = []
    for i,j in combinations(idx,2):
        face = [v for v in idx if v not in (i,j)]
        size = arb(gram(distances,face).det()).sqrt()/factorial(n-2)
        cos = -arb(normals[i,j])/arb(normals[i,i]*normals[j,j]).sqrt()
        steps = tuple(masks[b]-masks[a] for a,b in zip(face, face[1:]))
        hinges.append((masks[face[0]], steps, size, angle_turns(cos)))
    return g.det()/factorial(n)**2, hinges


def triangulation(mode, n, t):
    faces, counts = {}, Counter()
    volume = arb(0)
    reference_squared, _ = local_geometry(mode, 1, tuple(range(n)), t)
    for point in product(range(SIDE), repeat=n):
        parity = (-1)**sum(point[:3])
        for order in permutations(range(n)):
            squared, hinges = local_geometry(mode, parity, order, t)
            assert squared == reference_squared
            volume += arb(squared).sqrt()
            for start, steps, size, angle in hinges:
                origin = tuple((x+((start>>i)&1)) % SIDE for i,x in enumerate(point))
                key = (origin, steps)
                if key not in faces:
                    faces[key] = (size, arb(1))
                fixed_size, deficit = faces[key]
                assert fixed_size.overlaps(size)
                faces[key] = (fixed_size, deficit-angle)
                counts[key] += 1
    assert len(faces) == SIDE**n * (7 if n == 3 else 50)
    assert set(counts.values()) == {4,6}
    action = sum(area*deficit for area,deficit in faces.values())
    return volume, action, faces


def verify():
    tol = fmpq(1, 10**8)
    for mode, bounds, vols in [
        ("parity", (152,153), (31,32)), ("body", (-13,-12), (32,33))
    ]:
        v3, s3, base = triangulation(mode, 3, 1)
        assert vols[0] < v3 and v3 < vols[1]
        assert fmpq(bounds[0],1000) < s3 and s3 < fmpq(bounds[1],1000)
        print(mode, "V3", v3, "S3 in turns", s3)
        for t in (1,2):
            v4, s4, faces = triangulation(mode, 4, t)
            assert abs(v4-4*v3) < tol and abs(s4-4*t*t*s3) < fmpq(1,10**6)
            for (origin, steps), (area, deficit) in faces.items():
                if 8 in steps:
                    horizontal = next(s for s in steps if s != 8)
                    edge_length, expected = base[origin[:3], (horizontal,)]
                    assert abs(area-edge_length*t*t/2) < tol
                else:
                    expected = 0
                assert abs(deficit-expected) < tol
            print("  t", t, "V4", v4, "S4 in turns", s4)
    eps = 192*fmpq(1,2**40)
    incidence = SIDE**4 * factorial(4) * 10
    # Nine eps t^2 per triangle bounds |area'|; 1 + degree/2 bounds |deficit|.
    assert 9*eps*(50*SIDE**4 + incidence//2) < fmpq(1,1000)
    assert 4*fmpq(12,1000)-fmpq(1,1000) > fmpq(1,25)


if __name__ == "__main__":
    verify()
