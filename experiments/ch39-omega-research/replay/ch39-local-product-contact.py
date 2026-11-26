"""Two lifted 25-site shells with the shortest connected chord cutoff and unit edges."""
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpz_mat, fmpz_poly

contact = run_path(str(Path(__file__).with_name("ch39-624-clock-contact.py")))
geometry = run_path(str(Path(__file__).with_name("ch39-intrinsic-shell-comparisons.py")))


def lifted_sites():
    shell = [p for p in product(range(30),repeat=3) if sum(v*v for v in p)==867]
    pole = (17,17,17)
    distance = geometry["squared"]
    sites = [p for p in shell for _ in range(3 if distance(p,pole)==68 else 1)]
    assert len(sites)==25
    return sites, pole


def lifted_laplacian():
    sites, pole = lifted_sites()
    distance = geometry["squared"]
    def laplacian(cutoff):
        rows = [[-int(0<distance(p,q)<=cutoff) for q in sites] for p in sites]
        for i,row in enumerate(rows):
            row[i] = -sum(row)
        return fmpz_mat(rows)
    assert laplacian(151).rank()==21
    matrix = laplacian(152)
    assert matrix.rank()==24
    return matrix, sites.index(pole)


def graph_spectrum():
    matrix,port = lifted_laplacian()
    interior = [i for i in range(25) if i!=port]
    whole = matrix.charpoly()
    minor = fmpz_mat([[matrix[i,j] for j in interior] for i in interior]).charpoly()
    common = whole.gcd(minor)
    denominator,numerator = whole//common,minor//common
    assert denominator == fmpz_poly([0,-175,199,-29,1])
    assert numerator == fmpz_poly([-7,19,-11,1])
    for root,_ in whole.complex_roots():
        assert root.is_zero() or root.real > arb(".1")
    poles = []
    for root,mult in denominator.complex_roots():
        assert mult==1 and root.imag.is_zero()
        energy = root.real
        weight = numerator(energy)/denominator.derivative()(energy)
        assert weight>0
        poles.append((energy,weight))
    modes = []
    for i,j in combinations_with_replacement(range(len(poles)),2):
        x,p = poles[i]
        y,q = poles[j]
        if x.is_zero() and y.is_zero():
            assert arb(1)/625 in p*q
        else:
            modes.append((x+y,(1 if i==j else 2)*625*p*q))
    assert len(modes)==9
    assert all(0 not in modes[i][0]-modes[j][0] for i in range(9) for j in range(i))
    moments = [sum(w*e**k for e,w in modes)/625 for k in (1,2,3)]
    print("local product pole-port channel moments",*(m.str(10) for m in moments))
    return modes


def verify():
    ctx.prec=192
    modes = graph_spectrum()
    shift,effective,excited,probability,residual = contact["ground_readouts"](contact["clock_window"](),modes)
    assert arb(".00155526915") < effective < arb(".00155526916")
    assert shift>0 and excited<arb("1.446e-10")
    assert probability>arb(".999866") and residual<arb("2.97e-7")


if __name__=="__main__":
    verify()
