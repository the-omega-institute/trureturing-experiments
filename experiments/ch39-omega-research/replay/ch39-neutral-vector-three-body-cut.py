"""十九点壳三个轮换不变中性矢量的两圈超荷产生阈值。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat, fmpz_poly

shell = run_path(str(Path(__file__).with_name("ch39-unit-seventy-five-crossing.py")))


def invariant_masses():
    soft = shell["soft"]
    h = acb.zeta_zero(11431).imag
    base, neutral, _ = soft["vector_pencil"](h)
    rows = [[int(t.real.unique_fmpz()) for t in row] for row in (base + neutral).tolist()]
    rotation = soft["star"]["charged"]["rotation"]
    assert all(rotation[rotation[rotation[i]]] == i for i in range(19))
    orbits = sorted({tuple(sorted({i, rotation[i], rotation[rotation[i]]})) for i in range(19)})
    assert len(orbits) == 7
    assert all(rows[i][j] == rows[rotation[i]][rotation[j]] for i, j in product(range(19), repeat=2))
    spokes = set(soft["star"]["orbit"])
    assert len(spokes) == 3 and {rotation[i] for i in spokes} == spokes
    assert rotation[soft["pole"]] == soft["pole"]
    quotient = fmpz_mat([[sum(rows[group[0]][j] for j in target) for target in orbits]
                         for group in orbits])
    x = fmpz_poly([0, 1])
    assert quotient.charpoly() == x*(x-11)*(x-17)*(x-19)**2*(x*x-28*x+191)
    # On the invariant subspace, the same three-spoke star perturbation is
    # negative semidefinite, with norm at most four times the edge decrement.
    decrement = arb(23)/25 * arb(2)**-20 / (h*h + soft["point"]["study"]["offset"])
    assert 0 < 4*decrement < arb("1e-12")
    r = soft["geo"]["theta"](804) / soft["geo"]["theta"](12)
    loss = arb(0).union(4*decrement)
    return [(mass_sq-loss).sqrt()/r for mass_sq in (arb(11), 14-arb(5).sqrt(), 14+arb(5).sqrt())]


def verify():
    ctx.prec = 224
    a, b, c = invariant_masses()
    assert 0 < a < b < c
    threshold = a+b+c
    assert arb("1.26253555") < threshold < arb("1.26253557")
    # Each F.G at rest equals -2 E_F E_G times the spatial polarization dot
    # product. Keep the incident source along x and sum three final spins.
    spin_sum = sum(int(i == 0 and j == k)**2 for i, j, k in product(range(3), repeat=3))
    assert spin_sum == 3
    amplitude = 4*threshold*a*b*c
    phase_space = ((a*b*c/threshold**3).sqrt())/(64*arb.pi()**2)
    weight = spin_sum*amplitude**2*phase_space
    assert arb(".00012599292") < weight < arb(".00012599295")
    # Maximize z*(1-z**3) for z=M_threshold**2/Lambda**2 in [0,1].
    optimum = (arb(1)/4)**(arb(1)/3)
    cut_bound = optimum*(1-optimum**3)/(32*arb.pi()**3*threshold**2)
    assert cut_bound < arb(".000299") < arb(".0628")
    print("Three neutral masses / matching scale", *(mass.str(13) for mass in (a,b,c)))
    print("Three-body threshold / matching scale", threshold.str(13))
    print("Leading phase-space integral coefficient in matching units", weight.str(13))
    print("Upper envelope of finite-cut inverse hypercharge shift per c squared", cut_bound.str(13))


if __name__ == "__main__":
    verify()
