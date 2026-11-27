"""在两枚有限编号实际零点上检验六质量和四项混合参数。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx

def load(name):
    return run_path(str(Path(__file__).resolve().parent / f"ch39-{name}.py"))

shell = load("three-ray-chiral-boundary")
pump = load("native-light-projector-pump")
model = load("moving-preparation-quark-shape")
refs = load("shell-flavor-matching-windows")
ctx.prec = 224
adj = pump["adjoint"]
sources = model["orthonormal_sources"](shell, adj)
frames = list(pump["rotation_bases"](shell))
radius = arb("1e-11")
pad = arb(0, radius.upper())
tau = 2 * arb.pi()

def theta(d):
    return (1 - arb(d)/1734).acos()

def hardy(t):
    phase = acb(arb(1)/4, t/2).lgamma().imag - t*arb.pi().log()/2
    return (acb(0, phase).exp()*acb(arb(1)/2, t).zeta()).real

def locate(q, x, y):
    a, b = [abs(q)*(theta(d)-theta(244)) for d in (724, 804)]
    if q < 0:
        x, y = -x, -y
    assert radius < x < tau-radius
    stop = int(((a*1400000-x+radius)/tau).upper().ceil().unique_fmpz()) + 1
    hits = []
    for k in range(stop):
        heights = (x+k*tau+pad)/a
        wraps = (b*heights-y+pad)/tau
        assert wraps.upper()-wraps.lower() < 1
        j = wraps.unique_fmpz()
        if j is not None:
            hits.append((k, int(j)))
    [(k, j)] = hits
    midpoint = ((x+k*tau)/a + (y+j*tau)/b)/2
    m = 10**13
    t = arb((m*midpoint+arb(1)/2).floor().unique_fmpz())/m
    e = arb("1e-12")
    heights = t + arb(0, e.upper())
    assert 500000 < heights < 1400000
    assert abs(a*heights-x-k*tau) < radius
    assert abs(b*heights-y-j*tau) < radius
    lo, hi = t-e, t+e
    n, N = [v.zeta_nzeros().unique_fmpz() for v in (lo, hi)]
    assert N == n+1 and hardy(lo)*hardy(hi) < 0
    assert int(N) == {-4: 2345666, 2: 964890}[q]
    print(q, k, j, n, N, t.str(29))

def eig(H):
    values, vectors = H.eig(right=True, algorithm="rump")
    order = sorted(range(len(values)), key=lambda i: float(values[i].real.mid()))
    roots = [values[i].real for i in order]
    assert all(values[i].imag.contains(0) and values[i].real > 0 for i in order)
    assert all(x < y for x, y in zip(roots, roots[1:]))
    norms = [sum(abs(vectors[i,j])**2 for i in order).sqrt() for j in order]
    basis = acb_mat([[vectors[i,j]/n for j,n in zip(order,norms)]
                     for i in range(len(order))])
    return roots, basis

radial = acb_mat(19, 19)
pole = shell["soft"]["pole"]
for i in shell["soft"]["star"]["orbit"]:
    radial[pole,pole] += 1
    radial[i,i] += 1
    radial[pole,i] -= 1
    radial[i,pole] -= 1
d = arb("2e-18")
deficit = d + arb(0, d.upper())
assert arb(23)/25 * arb(2)**-20 / 500000**2 < arb("4e-18")

def prepare(q, x, y):
    L, edges = pump["clock_matrices"](shell, q)
    L -= deficit*radial
    for edge, angle in zip(edges, (x+pad, y+pad), strict=True):
        term = edge*acb(0, angle).exp()
        L += term + adj(term)
    rows = []
    for Q in frames:
        levels, basis = eig(adj(Q)*L*Q)
        assert levels[0] < arb("8.685") and levels[1] > arb("9.6")
        v = acb_mat([[basis[i,0]] for i in range(basis.nrows())])
        rows.append((adj(Q*v)*sources/(1+(levels[0]/9)**512)).entries())
    return acb_mat(rows)

entries = [
    "3.81536789374417292e-5,0.00292058312666154686 "
    "0.0118502329710027280,0.00727882767491516284 "
    "-1.50338639415880113,1.87408266165223205",
    "-0.000824770774018664215,-0.00167027614241661673 "
    "-0.00937626282071660403,-0.00155790541331349012 "
    "1.32720223364068534,-0.509342463568852674",
    "0.000763827459989976886,-0.00147705916306281873 "
    "-0.00479444779678322068,-0.00759858211198344168 "
    "0.212250573161628303,-1.50767757434267292",
]
K = acb_mat([[acb(*map(arb, z.split(","))) for z in row.split()]
             for row in entries])
assert abs(K.det()) > 0
singular = load("first-zero-relative-quark-spectrum")["singular_values"]
source_shape = singular(K)
assert arb("14551.17") < source_shape[2]/source_shape[0] < arb("14551.18")
outputs = []
for q, x, y, scale in [
    (-4, "-0.680061861483031", "0.873223965353125", "7.06539940386823251"),
    ( 2, "0.565247837594779", "-0.630343671648397", "1"),
]:
    x, y = arb(x), arb(y)
    locate(q, x, y)
    prepared = prepare(q, x, y)
    if q == -4:
        shape = singular(prepared)
        required = arb(".8573")/arb(".00000629") * shape[0]/shape[2]
        assert required > arb("76.31")
    Y = prepared*K*arb(scale)
    outputs.append(eig(adj(Y)*Y))

(u, U), (d, D) = outputs

def mixing_angles(U, D):
    V = adj(U)*D
    s13 = abs(V[0,2])
    c13 = (1-s13**2).sqrt()
    s12, s23 = abs(V[0,1])/c13, abs(V[1,2])/c13
    quartet = V[0,1]*V[1,1].conjugate()*V[1,2]*V[0,2].conjugate()
    phase = (quartet/(s12*s23*s13*c13**2)+s12*s23*s13).arg()
    return [s12.asin(), s23.asin(), s13.asin(), phase]


values = [z.sqrt() for z in u+d] + mixing_angles(U,D)
ref = refs["up"] + refs["down"] + list(map(arb, (".22704", ".04275", ".00377", "1.139")))
tol = list(map(arb, ("1.4e-7", "5e-5", ".0043", "2e-7", "3e-6", ".00009",
                    ".00082", ".00042", ".00008", ".023")))
assert all(abs(x-y) < w/8 for x,y,w in zip(values, ref, tol, strict=True))
print(*(z.str(12) for z in values))
