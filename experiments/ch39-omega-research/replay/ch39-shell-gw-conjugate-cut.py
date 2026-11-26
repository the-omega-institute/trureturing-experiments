"""First spectral cut supporting a non-singlet profile together with its covariant adjoint."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

folder = Path(__file__).resolve().parent
a = run_path(str(folder/"ch39-shell-gw-asymmetric-bands.py"))
m = a["m"]
hs, adj, rows, count = m["hs"],m["adj"],m["rows"],m["count"]
sing, family, phases, gen = a["five"]["m"]["transported_spectra"]()
frame = a["band_frame"]


def matrix(cols):
    return acb_mat([[c[i,0] for c in cols] for i in range(cols[0].nrows())])


def singlet_phases():
    v = frame(sing,1,arb(".02"))[0]
    return [v[i,0]/abs(v[i,0]) for i in range(6)]


def conjugate(v,p):
    r = v.nrows()//6
    return acb_mat([[p[i//r]**2*v[i,0].conjugate()] for i in range(6*r)])


def lower_angle(v,p,bound):
    x = matrix(v)
    y = matrix([conjugate(q,p) for q in v])
    gram = adj(x)*x
    residual = y-x*gram.solve(adj(x)*y)
    form = adj(residual)*residual-bound**2*(adj(y)*y)
    assert count(rows(form),arb(0)) == 0


def rotation():
    # Cycle the spatial axes and conjugate the family fiber by the same cycle.
    P = acb_mat([[int(i==(j+1)%3) for j in range(3)] for i in range(3)])
    d = acb_mat([[hs(x,P*y*adj(P)) for y in gen] for x in gen])
    pts = m["m"]["pts"]
    orbit = [next(k for k in pts if pts[k]==(p[2],p[0],p[1])) for p in (pts[i] for i in range(6))]
    return acb_mat([[d[i%8,j%8] if i//8==orbit[j//8] else 0
                     for j in range(48)] for i in range(48)])


def multiplicities(v):
    V = matrix(v)
    T = rotation()
    g = adj(V)*V
    theta = acb(0,2*arb.pi()/3).exp()
    trace = m["m"]["trace"]
    out = []
    for eig in (acb(1),theta,theta.conjugate()):
        sector = (V+eig.conjugate()*(T*V)+eig*(T*(T*V)))/3
        count_char = trace(g.solve(adj(V)*sector))
        n = count_char.real.unique_fmpz()
        assert n is not None and count_char.imag.contains(0)
        out.append(int(n))
    return out


if __name__ == "__main__":
    p = singlet_phases()
    lo,hi = arb("3.9978200948"),arb("3.9978200951")
    assert [count(rows(sing),t) for t in (lo,hi)] == [2,2]
    assert [count(rows(family),t) for t in (lo,hi)] == [23,24]
    # In the two-dimensional singlet band the ground line is fixed by J;
    # rank three of these columns excludes a second common direction.
    v = frame(sing,2,lo)
    joined = matrix(v+[conjugate(v[1],p)])
    assert count(rows(adj(joined)*joined),arb(".1")) == 0
    lower_angle(frame(family,8,arb("1.105")),p,arb(".00096"))
    lower_angle(frame(family,23,lo),p,arb(".0000117"))
    # The singlet ground character is trivial: projection on the constant is nonzero.
    t = frame(sing,1,arb(".02"))[0]
    assert abs(sum(t[i,0] for i in range(6)))>1
    chars = multiplicities(frame(family,24,hi))
    assert chars == [7,8,9]
    print("first common adjoint at",lo.str(14),hi.str(14), "cycle dimensions",chars,
          "two common adjoint directions, plus the singlet ground line")
