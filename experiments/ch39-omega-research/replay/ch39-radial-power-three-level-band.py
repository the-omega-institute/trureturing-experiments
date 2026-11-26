"""Calibrate a radial power at two old sites and certify the top-three slope on a band."""
from collections import Counter
from fractions import Fraction as F
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx
from sympy import Poly, Rational, symbols

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))
ball = shell["ball"]


def rotation_bases(nodes):
    orbits = []
    for p in nodes:
        orbit = [p,p[1:]+p[:1],p[2:]+p[:2]]
        if p == min(orbit):
            orbits.append([nodes.index(q) for q in dict.fromkeys(orbit)])
    omega = acb(-1,arb(3).sqrt())/2
    for k in range(3):
        groups = [g for g in orbits if k == 0 or len(g) == 3]
        u = acb_mat(len(nodes),len(groups))
        for j,g in enumerate(groups):
            for i,v in enumerate(g):
                u[v,j] = omega**(k*i)/arb(len(g)).sqrt()
        yield u.conjugate().transpose(), u


def calibrate(rows):
    outer, inner = rows
    def readings(row,p):
        return [(n*(-2*p*l).exp(),l) for l,n in row]
    def log_mean(row, p):
        vals = readings(row,p)
        return sum(w*l for w,l in vals)/sum(w for w,l in vals)
    for j in range(128):
        lo, hi = ball(F(j,64)), ball(F(j+1,64))
        assert 2*(log_mean(outer,hi)-log_mean(inner,lo)) > arb(".35")
    assert 2*(min(l for l,n in outer)-log_mean(inner,arb(2))) > arb(".35")
    def ratio(p):
        a,b = [sum(w for w,_ in readings(row,p)) for row in rows]
        return b/a
    reference = ratio(arb(3))
    error = arb(".0025")
    stretch = ((1+error)/(1-error))**2
    assert ratio(arb("2.99")) < reference/stretch and reference*stretch < ratio(arb("3.01"))
    print("Cubic two-site ratio",reference.str(18))


def verify():
    ctx.prec = 192
    nodes, edges = shell["geometry"]()
    distances = [shell["LOCAL"]["squared"](nodes[i],nodes[j]) for i,j,_,_ in edges]
    logs = {h: (arb(h)/108).log() for h in set(distances)}
    rows = []
    for point in ((3,15,87),(39,39,69)):
        i = nodes.index(point)
        counts = Counter(h for h,(u,v,_,_) in zip(distances,edges) if i in (u,v))
        rows.append([(logs[h],m) for h,m in counts.items()])
    calibrate(rows)
    # Perturbation from M_3(0) on the two complex parameter domains.
    for radius, phase_factor in [(arb(".06"),arb(4)/3),(arb(".01"),arb(2))]:
        row_sums = [arb(0) for _ in nodes]
        for h,(i,j,w,k) in zip(distances,edges):
            stretch = (radius*logs[h]).exp()
            cost = ball(w)*(stretch-1 + abs(k)*phase_factor*stretch)
            row_sums[i] += cost
            row_sums[j] += cost
        assert max(row_sums) < arb(".0045")
    transformations = list(rotation_bases(nodes))
    t_nodes = [F(i,112) for i in range(5)]
    s_nodes = [t*t for t in t_nodes]
    variable = symbols("s")
    polynomials = []
    for j, sj in enumerate(s_nodes):
        basis = prod((variable-Rational(v))/Rational(sj-v) for i,v in enumerate(s_nodes) if i!=j)
        polynomials.append(Poly(basis.diff(variable),variable))
    T = s_nodes[-1]
    sbox = arb(0).union(ball(T))
    # On |t|<=1/2 the cluster sum is bounded by 3(.0796+.0045)<.253.
    radius = F(1,4)
    dist = radius-T
    s_remainder = ball(radius*F(".253")/dist**5*(5*T**4/dist+T**5/dist**2))
    p_nodes = [F(j,4) for j in range(-4,5)]
    slopes = []
    for v in p_nodes:
        p = 3+ball(v)/100
        weights = {h:(-p*l).exp()/19 for h,l in logs.items()}
        readings = []
        for t in t_nodes:
            it = acb(0,ball(t))
            phase = (1+it)/(1-it)
            mat = acb_mat(len(nodes),len(nodes))
            for h,(i,j,_,k) in zip(distances,edges):
                entry = weights[h]*phase**k
                mat[i,j],mat[j,i] = entry,entry.conjugate()
            total = arb(0)
            for adj,u in transformations:
                eig = (adj*mat*u).eig(multiple=True)
                high = [z.real for z in eig if z.real > arb(".075")]
                assert len(high) == 1 and all(z.real < arb(".075") or z.real > arb(".075") for z in eig)
                if v == 0 and t == 0:
                    assert arb(".0795") < high[0] < arb(".0796")
                    assert all(z.real < arb(".070") for z in eig if z.real < arb(".075"))
                total += high[0]
            readings.append(total)
        coeffs = [sum(y*ball(F(poly.nth(k))) for y,poly in zip(readings,polynomials)) for k in range(4)]
        derivative = arb(0)
        for vcoeff in reversed(coeffs):
            derivative = derivative*sbox+vcoeff
        slopes.append(derivative+arb(0,s_remainder))
    # Interpolate the s-derivative in the normalized exponent v=100(p-3).
    # On |v|=6 its modulus is <=.253/(.16-T) by the Cauchy derivative bound.
    bound = ball(F(".253")/(F(".16")-T))
    denominator = 6*prod(36-v*v for v in p_nodes if v>0)
    middle = arb("-7.4e-7")
    enclosed = []
    for j in range(128):
        x = ball(F(j,64)-1).union(ball(F(j+1,64)-1))
        omega = prod(x-ball(v) for v in p_nodes)
        remainder = 6*bound*abs(omega)/(5*ball(denominator))
        val = middle
        for k,v in enumerate(p_nodes):
            lagrange = prod((x-ball(q))/ball(v-q) for i,q in enumerate(p_nodes) if i!=k)
            val += lagrange*(slopes[k]-middle)
        enclosed.append(val+arb(0,remainder))
    assert all(v < arb("-5e-7") for v in enclosed)
    print("Uniform derivative upper bound",max(v.upper() for v in enclosed).str(12))


if __name__ == "__main__":
    verify()
