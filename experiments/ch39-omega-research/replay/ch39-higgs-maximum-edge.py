"""Rank requirements for a uniform edge tolerance on the rough five-periodic background."""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, arb_poly, ctx

base = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))
momenta = ((4, 1, 2, 1), (4, 0, 3, 2), (4, 3, 4, 0))
coordinates = (
    ((636412708781, 0), (48625353908, -166016867087), (138635324155, -112192968678),
     (-521546001003, 50904768426), (-222664805177, -457235914494)),
    ((209062754470, 0), (129847868393, 525080035751), (-306746095537, -534348288516),
     (-237298552124, -367253649170), (301595307294, -44231393886)),
    ((571254099659, 0), (388435911453, -58115750138), (109152019914, 25239549353),
     (366127783063, 307630377781), (277024106436, 448808755695)),
)
single = ((-61395848, 0), (-13052123, 20195291), (-15258426, 14055249),
          (56223357, 2784466), (20974794, 40164803))
adjoint = base["adjoint"]


def norm(vector):
    return (adjoint(vector)*vector)[0, 0].real.sqrt()


def normalized(entries):
    v = acb_mat([[acb(a, b)] for a, b in entries])
    return v / norm(v)


def verify():
    ctx.prec = 180
    lower, upper = arb("0.50898378605"), arb("0.50898378607")
    rank_cut = arb("0.50898379")
    gap_cut = arb("0.54")
    weights = [arb(n)/10**12 for n in (190228603135, 393680674736, 194536323233, 221554398896)]
    mixture = [arb(n)/10**12 for n in (701252351895, 57318735330, 241428912775)]
    c = (arb(5).sqrt()-1)/4
    zeta = acb(c, (1-c*c).sqrt())
    parts = []
    for link in base["field"]()[:4]:
        phases = [link*zeta**j for j in range(5)]
        parts.append([2*base["identity"](5)-u-adjoint(u) for u in phases])

    vectors = [normalized(row) for row in coordinates]
    energies = {}
    for k in product(range(5), repeat=4):
        actions = [parts[a][k[a]] for a in range(4)]
        H = sum((actions[a]*weights[a] for a in range(1, 4)), actions[0]*weights[0])
        poly = arb_poly([z.real for z in H.charpoly().coeffs()])
        assert base["positive_roots"](poly, lower) == 5
        assert base["positive_roots"](poly, gap_cut) == 5-int(k in momenta)
        if k in momenta:
            v = vectors[momenta.index(k)]
            energies[k] = [(adjoint(v)*op*v)[0, 0].real for op in actions]
            # All but the target line have eigenvalue > 0.54 in this block.
            assert norm((H-base["identity"](5)*lower)*v) < arb("0.0000000001")
            if k == momenta[0]:
                u = normalized(single)
                assert all((adjoint(u)*op*u)[0, 0].real < arb("0.53652") for op in actions)

    for a in range(4):
        assert sum(mixture[i]*energies[k][a] for i, k in enumerate(momenta)) < upper
    assert 2*arb("0.0000000001")/(gap_cut-lower) < arb("0.0000001")
    for i, j in combinations(range(3), 2):
        assert abs((adjoint(vectors[i])*vectors[j])[0, 0]) > arb("0.065")
    differences = {tuple((x-y)%5 for x, y in zip(k, l)) for k in momenta for l in momenta if k != l}
    assert len(differences) == 6

    tail = arb("0.00036")
    assert (rank_cut-lower)/(gap_cut-lower) < tail**2
    assert (2*tail+tail**2)/arb("0.064") < arb("0.012")
    assert rank_cut-(rank_cut-lower)/min(weights) > arb("0.50898")
    # Each directional energy differs from its low-space contribution by at most 8*tail.
    target = arb("0.509", "0.003")
    assert target.contains(arb("0.50898")-8*tail) and target.contains(rank_cut+8*tail)
    matrix = arb_mat([[1]*3] + [
        [energies[k][a] + arb(0, "0.000001") for k in momenta] for a in (0, 2)
    ])
    rhs = arb_mat([[arb(1, "0.000001")], [target], [target]])
    masses = matrix.solve(rhs)
    assert all(masses[i, 0] > arb("0.04") for i in range(3))
    print("Squared maximum edge cost for three copies in (0.50898378605, 0.50898378607).")
    print("At most two copies exceed 0.50898379; one copy attains a cost below 0.53652.")


if __name__ == "__main__":
    verify()
