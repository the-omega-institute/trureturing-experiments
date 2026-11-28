"""把夸克十二次有理作用不再校准地施于原壳中性读出。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    model = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    ctx.prec = 256
    roots = [arb(r) for r in model["ROOTS"]]
    _, block, _ = pump["rotation_bases"](shell)
    adj = pump["adjoint"]
    base, links = pump["clock_matrices"](shell,0)
    edges = links[0]+links[1]
    flat = base+edges+adj(edges)
    soft = shell["soft"]
    pole = soft["pole"]
    radial = acb_mat(19,19)
    for site in soft["star"]["orbit"]:
        radial[pole,pole] += 1
        radial[site,site] += 1
        radial[pole,site] -= 1
        radial[site,pole] -= 1
    maximum = arb(23)/25*arb(2)**-20/(acb.zeta_zero(1).imag**2+soft["point"]["study"]["offset"])
    assert 0 < maximum < arb("4e-9")
    hamiltonian = flat-arb(0).union(arb("4e-9"))*radial
    vector, gap = pump["first_line"](adj(block)*hamiltonian*block)
    vector = block*acb_mat([[v] for v in vector])
    energy = (adj(vector)*hamiltonian*vector)[0,0].real
    overlap = (adj(vector)*radial*vector)[0,0].real
    derivative = 2*prod((energy-r)**2 for r in roots)*sum(1/(energy-r) for r in roots)
    assert gap > arb("2.75") and overlap > arb(".45617")
    assert arb("-1.758") < derivative < arb("-1.757")
    print("中性非平凡字符能量、径向期望、十二次势导数",energy.str(11),overlap.str(11),derivative.str(11))


if __name__ == "__main__":
    verify()
