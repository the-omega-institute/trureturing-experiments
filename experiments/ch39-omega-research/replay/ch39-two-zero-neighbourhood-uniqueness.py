"""两部门各取中心前后六十四个实际零点，排除其余全部可调换基组合。"""
from pathlib import Path
from itertools import product
from runpy import run_path
from flint import acb, acb_mat, arb, ctx

def load(name):
    return run_path(str(Path(__file__).with_name(f"ch39-{name}.py")))

m = load("two-zero-quark-window")
load("moving-preparation-quark-shape")["target_gaps"]()
ctx.prec = 224
adj, eig, graph = (m[k] for k in ("adj", "eig", "shell"))

def prepare(q, n):
    midpoint = acb.zeta_zero(n).imag.mid()
    e = arb("1e-12")
    lo, hi = midpoint-e, midpoint+e
    assert [x.zeta_nzeros().unique_fmpz() for x in (lo, hi)] == [n-1, n]
    assert m["hardy"](lo)*m["hardy"](hi) < 0
    h = midpoint + arb(0, e.upper())
    assert h < 1311000
    L = graph["laplacian"](q, h)
    rows = []
    for Q in m["frames"]:
        levels, basis = eig(adj(Q)*L*Q)
        assert levels[0] < arb("8.685") and levels[1] > arb("9.6")
        v = acb_mat([[basis[i, 0]] for i in range(basis.nrows())])
        rows.append((adj(Q*v)*m["sources"]/(1+(levels[0]/9)**512)).entries())
    A = acb_mat(rows)
    inverse = A.inv()
    assert sum(abs(z)**2 for z in inverse.entries()) < arb(4430)**2
    return h, A, inverse

theta = graph["geo"]["theta"]
speeds = [theta(d)-theta(244) for d in (724, 804)]
tau = 2*arb.pi()
offsets = range(-64, 65)
up, down = [[prepare(q, n+i) for i in offsets]
            for q, n in ((-4, 2345666), (2, 964890))]
eps = (arb(9)/7)*4430/(1+(arb("9.6")/9)**512)
assert 0 < eps < arb("2.55e-11")
assert (arb(9)/7)**2*((1+eps)/(1-eps))**2 < arb("1.654")

count = 0
for i, j in product(offsets, repeat=2):
    if (i, j) == (0, 0):
        continue
    hu, A, _ = up[i+64]
    hd, _, inverse = down[j+64]
    angles = [w*(hu-hd) for w in speeds]
    turns = [(x/tau+arb(1)/2).floor().unique_fmpz() for x in angles]
    assert all(k is not None for k in turns)
    if any(abs(x-k*tau) < arb(".005") for x, k in zip(angles, turns)):
        count += 1
        continue
    a, b, c = m["singular"](A*inverse)
    assert b/a < arb("77.12")*arb(".107") or c/b < arb("61.32")*arb(".107")
assert count == 49
print(count, len(offsets)**2-1-count)
