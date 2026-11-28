"""Regge action, volume and fatness on independent edge boxes of the 4^4 torus."""

from collections import Counter
from itertools import combinations
from math import comb, factorial

from flint import arb, arb_mat, ctx, fmpq


SIDE = 4
FACES = []
for i, j in combinations(range(5), 2):
    a, b, c = sorted(set(range(5)) - {i, j})
    p, q = b-a, c-b
    valence = factorial(p)*factorial(q)*factorial(5-p-q)
    FACES.append((i, j, a, b, c, valence))


# The first tuple specifies squared lengths by step count. The second
# coefficient multiplies the sum of the endpoint checkerboard signs.
PROFILES = {
    "equilateral": ((1, 1, 1, 1), 0),
    "checkerboard": ((fmpq(4, 5), fmpq(6, 5), fmpq(6, 5), fmpq(4, 5)), fmpq(1, 25)),
    "negative": ((1, 1, fmpq(3242533, 10**7), 1), 0),
    "positive": ((1, fmpq(13, 10), fmpq(59, 40), fmpq(6, 5)), fmpq(11, 100)),
}


def chain(parity, width, profile):
    lengths, amplitude = PROFILES[profile]
    edges = [[arb(0) for _ in range(5)] for _ in range(5)]
    ceiling = arb(0)
    for i, j in combinations(range(5), 2):
        step = j-i
        base = lengths[step - 1] * (1 + amplitude * parity * ((-1)**i + (-1)**j))
        edges[i][j] = edges[j][i] = arb(base, width)
        ceiling = max(ceiling, edges[i][j].upper())

    gram = arb_mat([[(edges[0][i]+edges[0][j]-edges[i][j])/2 for j in range(1, 5)]
                    for i in range(1, 5)])
    for k in range(1, 5):
        assert arb_mat([[gram[i, j] for j in range(k)] for i in range(k)]).det() > 0
    determinant = gram.det()
    fat = determinant/(ceiling**4)
    assert fat > fmpq(1, 16)
    change = arb_mat([[-1]*4]+[[int(i == j) for j in range(4)] for i in range(4)])
    normals = change*gram.inv()*change.transpose()
    action = arb(0)
    for i, j, a, b, c, valence in FACES:
        ab, ac, bc = edges[a][b], edges[a][c], edges[b][c]
        area = (ab*ac-((ab+ac-bc)/2)**2).sqrt()/2
        cosine = -normals[i, j]/(normals[i, i]*normals[j, j]).sqrt()
        assert -1 < cosine < 1
        action += area*(2*arb.pi()/valence-cosine.acos())
    return action, determinant.sqrt()/24, fat


def verify_box(profile, width):
    a, v, fat_a = chain(1, width, profile)
    b, w, fat_b = chain(-1, width, profile)
    block = SIDE**4 * factorial(4)//2
    action = block*(a+b)
    volume = block*(v+w)
    shape = action/volume.sqrt()
    assert 0 < volume < 256
    print(profile, "half width", width,
          "T endpoints", shape.lower().str(12), shape.upper().str(12), "V", volume.str(10),
          "fatness lower", min(fat_a.lower(), fat_b.lower()).str(9))
    return shape


def main():
    ctx.prec = 192
    # Each triangle's 2pi budget is divided among its incident simplices.
    incidence = Counter(n for *_, n in FACES)
    assert incidence == {4: 5, 6: 5}
    assert sum(fmpq(1, n) for *_, n in FACES)*24 == 50

    for profile, bounds in (("equilateral", ("-20.281", "-20.280")),
                            ("checkerboard", ("0.62804", "0.62805")),
                            ("negative", ("-153.826", "-153.825")),
                            ("positive", ("2.23", "2.231"))):
        value = verify_box(profile, 0)
        assert arb(bounds[0]) < value < arb(bounds[1])
    cases = (("equilateral", 5, -19, 20456), ("equilateral", 4, -7, 16616),
             ("checkerboard", 8, fmpq(3, 5), 31976), ("checkerboard", 6, fmpq(1, 3), 24296),
             ("negative", 10, fmpq(-153825, 1000), 39656), ("positive", 8, 2, 31976))
    for profile, digits, bound, exponent in cases:
        width = fmpq(1, 10**digits)
        value = verify_box(profile, width)
        assert value < bound if bound < 0 else value > bound

        # Scale the whole box by 1/16 in squared lengths to ensure V<=1.
        # Fatness bounds the enclosing weighted ellipsoid in its positive orthant.
        edges = 15*SIDE**4
        half_dim = edges//2
        ellipsoid_volume = (arb.pi()*960/4)**half_dim/factorial(half_dim)
        for steps in range(1, 5):
            valence = factorial(steps)*factorial(5-steps)
            count = SIDE**4*comb(4, steps)
            ellipsoid_volume /= arb(valence)**(count//2)
        mass = (arb(width)/8)**edges / ellipsoid_volume
        assert mass > arb(10)**(-exponent)
        print("probability exceeds 10^", -exponent)


if __name__ == "__main__":
    main()
