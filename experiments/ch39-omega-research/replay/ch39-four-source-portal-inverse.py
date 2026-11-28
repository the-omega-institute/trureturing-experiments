"""Finite radial sampling and sparse four-layer phase reconstruction checks."""

from itertools import combinations, permutations, product

from flint import arb, arb_mat, ctx, fmpq


def radial_samples():
    ctx.prec = 256
    rows = []
    for k in range(5):
        r = arb(fmpq(100 + k, 100)).sqrt()
        row = []
        for mass in (1, 2, 3, 4):
            row.append(mass * (mass*r).bessel_k(1) / (4*arb.pi()**2*r))
        rows.append(row)
    triangle = arb_mat(rows[:4])
    quartet = arb_mat([[b-a for a, b in zip(left, right)]
                       for left, right in zip(rows, rows[1:])])
    for name, matrix in (("triangle", triangle), ("quartet differences", quartet)):
        determinant = matrix.det()
        assert not determinant.contains(0)
        print(name, "determinant", determinant)


def closed_cycles():
    cycles = []
    for size in (3, 4):
        for subset in combinations(range(4), size):
            for tail in permutations(subset[1:]):
                if tail[0] < tail[-1]:
                    path = (subset[0],) + tail
                    cycles.append(list(zip(path, path[1:] + path[:1])))
    return cycles


def support_check():
    edges = list(combinations(range(4), 2))
    cycles = closed_cycles()
    cosines = (1, 0, -1, 0)
    classes = 0
    for included in product((False, True), repeat=len(edges)):
        support = {e for e, present in zip(edges, included) if present}
        components = [{i} for i in range(4)]
        chords = []
        for e in sorted(support):
            left, right = [next(c for c in components if v in c) for v in e]
            if left == right:
                chords.append(e)
            else:
                left.update(right)
                components.remove(right)
        words = [word for word in cycles
                 if all(tuple(sorted(e)) in support for e in word)]
        seen = {}
        for phase in product(range(4), repeat=len(chords)):
            angles = dict(zip(chords, phase))
            key = []
            for word in words:
                winding = sum((1 if i < j else -1)*angles.get(tuple(sorted((i,j))), 0)
                              for i, j in word)
                key.append(cosines[winding % 4])
            orbit = min(phase, tuple(-x % 4 for x in phase))
            assert seen.setdefault(tuple(key), orbit) == orbit
        classes += len(seen)
    print("support graphs", 2**len(edges), "quarter-turn phase orbits", classes)


if __name__ == "__main__":
    radial_samples()
    support_check()
