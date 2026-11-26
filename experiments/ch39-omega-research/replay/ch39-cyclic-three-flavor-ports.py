"""三个等价味端口在同一十九点壳上相干传播时的全相位反例。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

point = run_path(str(Path(__file__).with_name("ch39-port-coupling-competition.py")))
study, native = point["study"], point["native"]
geo = study["geo"]
rotation = [geo["sites"].index((y, z, x)) for x, y, z in geo["sites"]]
outer = (0, rotation[0], rotation[rotation[0]])
inside = [i for i in range(19) if i not in outer]


def responses(index, shift):
    height = acb.zeta_zero(index).imag
    shell = geo["matrix"](0, geo["phases"](height))
    block = acb_mat([[shell[i, j]+(shift if i == j else 0) for j in inside] for i in inside])
    omega = acb(-1, arb(3).sqrt())/2
    matrices = []
    for character in range(3):
        source = acb_mat(len(inside), 2)
        near, far = native["portals"]
        for r, origin in enumerate(outer):
            for arm, site in enumerate((near, far)):
                factor = omega**(character*r)*shell[site, origin]/arb(3).sqrt()
                for k in range(2):
                    source[inside.index(site), k] += factor*native["columns"][k, arm]
            near, far = rotation[near], rotation[far]
        matrices.append(point["adjoint"](source)*block.solve(source))
    return matrices, 1/(height**2+study["offset"])


def check_three_columns():
    height = acb.zeta_zero(11431).imag
    shell = geo["matrix"](0, geo["phases"](height))
    block = acb_mat([[shell[i, j] for j in inside] for i in inside])
    port = acb_mat([[arb(3)/5], [acb(0, 4)/5]])
    near, far = native["portals"]
    first = acb_mat(len(inside), 1)
    for k, site in enumerate((near, far)):
        first[inside.index(site), 0] = (native["columns"].transpose()*port)[k, 0]*shell[site, 0]
    source = acb_mat(len(inside), 3)
    for i, start in enumerate(inside):
        site = start
        for r in range(3):
            source[inside.index(site), r] = first[i, 0]
            site = rotation[site]
    gram = point["adjoint"](source)*block.solve(source)
    matrices, _ = responses(11431, 0)
    omega = acb(-1, arb(3).sqrt())/2
    fourier = acb_mat([[omega**(i*j)/arb(3).sqrt() for j in range(3)] for i in range(3)])
    diagonalized = point["adjoint"](fourier)*gram*fourier
    for i in range(3):
        for j in range(3):
            expected = (point["adjoint"](port)*matrices[i]*port)[0, 0] if i == j else 0
            assert abs(diagonalized[i, j]-expected) < arb("1e-40")


def verify():
    ctx.prec = 192
    check_three_columns()
    shift = arb(2)**-20
    centers, inv = responses(11431, 0)
    assert all(m.trace().real+128*inv+arb("2e-6") < 1 for m in centers)
    others = {n: responses(n, shift)[0] for n in (330381, 657986, 2019542)}
    weights = ((330381, 657986, 608), (657986, 2019542, 188), (330381, 657986, 609))
    lower = arb("1.074e-6")
    for m, (first, second, numerator) in enumerate(weights):
        weight = arb(numerator)/1000
        difference = weight*others[first][m]+(1-weight)*others[second][m]-centers[m]
        difference[0, 0] -= lower
        difference[1, 1] -= lower
        assert difference[0, 0].real > 0 and difference.det().real > 0
    assert lower-2*inv > arb("1.05e-6")
    assert lower-128*inv > arb("3.3e-8")
    print("Three positive definite 2x2 certificates cover every normalized complex port.")
    print("Effective advantage at gauge coefficients 2 and 128",
          (lower-2*inv).str(12), (lower-128*inv).str(12))


if __name__ == "__main__":
    verify()
