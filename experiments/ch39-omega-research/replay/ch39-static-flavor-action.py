"""Least imaginary antisymmetric action in the radial-preserving three-flavor model."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

shell = run_path(str(Path(__file__).with_name("ch39-shell-clock-ckm-blocks.py")))
p, cis = shell["radial"], shell["cis"]
v = shell["target_box"](arb(0), arb(0))
offsets, coeffs = shell["rephasing_equations"](v)
c2 = shell["axis_cos_sq"]
comm_norm = (c2*(1-c2)).sqrt()
seeds = (("3.778043", "1.714406"), ("5.799087", "1.646460"),
         (".595113", "4.613742"), ("2.386175", "4.555629"))


def equations(x, y, parameters=None):
    constants, rows = (offsets, coeffs) if parameters is None else parameters
    harmonics = (cis(x-y), cis(x), cis(y))
    values, slopes = [], []
    for a, row in zip(constants, rows):
        z, h, k = [coef*phase for coef, phase in zip(row, harmonics)]
        values.append([a+(z+h+k).real])
        slopes.append([-(z+h).imag, (z-k).imag])
    return arb_mat(values), arb_mat(slopes)


def around(point, r):
    return [u + arb(0, r.upper()) for u in point]


def certify_root(point, r):
    f, jac = equations(*point)
    frozen = jac.inv()
    _, varying = equations(*around(point, r))
    remainder = arb_mat([[1,0],[0,1]]) - frozen*varying
    resid = frozen*f
    lip = max(sum(abs(remainder[i,j]) for j in range(2)).upper() for i in range(2))
    assert lip < 1 and max(abs(resid[i,0]).upper() for i in range(2)) + lip*r < r
    return around(point, r)


def enclose_roots():
    coarse, fine = [], []
    for seed in seeds:
        point = list(map(arb, seed))
        coarse.append(certify_root(point, arb(".001")))
        for _ in range(4):
            f, jac = equations(*point)
            correction = jac.solve(f)
            point = [(u-correction[i,0]).mid() for i,u in enumerate(point)]
        tiny = certify_root(point, arb("1e-32"))
        assert all(a.contains_interior(b) for a,b in zip(coarse[-1],tiny))
        fine.append(tiny)
    pending = [(0,0,0)]
    boxes = 0
    while pending:
        level, i, j = pending.pop()
        boxes += 1
        unit = arb.pi()/(1 << level)
        x, y = around([(2*i+1)*unit, (2*j+1)*unit], unit)
        if any(a.contains_interior(x) and b.contains_interior(y) for a,b in coarse):
            continue
        f, _ = equations(x,y)
        if 0 in f[0,0] and 0 in f[1,0]:
            assert level < 18
            pending.extend((level+1,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("central rephasing roots", len(fine), "cover boxes", boxes)
    return fine


def rephased_matrix(root):
    right = (cis(root[0]), cis(root[1]), 1)
    rows = [[z*h for z,h in zip(row,right)] for row in v]
    left = [(sum(z*r for z,r in zip(row,p))/p[i]).conjugate() for i,row in enumerate(rows)]
    return [[l*z for z in row] for l,row in zip(left,rows)]


def oriented_coefficient(root):
    w = rephased_matrix(root)
    plane_trace = sum(w[i][i] for i in range(3)) - 1
    half_gap = (abs(plane_trace)/2).acos()
    skew = abs(w[0][2]-w[2][0])*arb(sum(z*z for z in p)).sqrt()/p[1]
    return half_gap/half_gap.sin()*skew/(2*comm_norm)


if __name__ == "__main__":
    full = shell["target_box"]()
    gap = abs(full[2][0]) - abs(full[0][2])
    robust = arb(sum(z*z for z in p)).sqrt()*gap/(2*p[1]*comm_norm)
    assert robust > arb(".03010")
    costs = [oriented_coefficient(root) for root in enclose_roots()]
    assert costs[3] > arb(".03898381516") and costs[3] < arb(".03898381517")
    assert all(cost > costs[3] for cost in costs[:3])
    print("full box lower bound", robust.lower().str(12))
    print("central branch coefficients", *(cost.str(13) for cost in costs))
