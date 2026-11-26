"""Sharp support endpoint for a positive spectrum imitating finite logarithmic data."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx, fmpq, fmpq_poly


def ball(value):
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def rational_log(value):
    # For the grids below, the odd atanh series has y <= 1/2.
    y = (value-1)/(value+1)
    assert 0 < y <= Fraction(1,2)
    count = 192
    term = 2*y
    partial = Fraction(0)
    for j in range(count):
        partial += term/(2*j+1)
        term *= y*y
    tail = term/((2*count+1)*(1-y*y))
    return ball(partial).union(ball(partial+tail))


def divided_moments(xs):
    count = len(xs)
    factors = [fmpq_poly([1,fmpq(a.numerator,a.denominator)]) for a in (x-1 for x in xs)]
    denominator = fmpq_poly([1])
    for factor in factors:
        denominator *= factor
    basis = [denominator//factor for factor in factors]
    moments = arb_mat([[row[j] for j in range(count)] for row in basis]).solve(
        arb_mat([[rational_log(x)/ball(x-1)] for x in xs]))
    return moments, arb_poly(denominator)


def gaussian_nodes(moments, order):
    hankel = arb_mat([[moments[i+j] for j in range(order)] for i in range(order)])
    coefficients = hankel.solve(arb_mat([[-moments[i+order]] for i in range(order)]))
    polynomial = arb_poly([*coefficients.entries(),arb(1)])
    nodes = []
    for z in polynomial.complex_roots(tol=1e-30):
        assert z.imag.contains(0) and 0 < z.real < 1
        nodes.append(z.real)
    nodes.sort(key=float)
    assert len(nodes)==order and all(a<b for a,b in zip(nodes,nodes[1:]))
    return nodes


def gaussian_spectrum(count):
    assert count >= 2 and count%2 == 0
    order = count//2
    xs = [1+Fraction(2*i,count) for i in range(1,count+1)]
    moments, den = divided_moments(xs)
    nodes = gaussian_nodes(moments.entries(),order)
    vandermonde = arb_mat([[s**i for s in nodes] for i in range(order)])
    weights = vandermonde.solve(arb_mat([[moments[i,0]] for i in range(order)])).entries()
    assert min(weights)>0
    for i in range(count):
        assert abs(sum(w*s**i for w,s in zip(weights,nodes))-moments[i,0]) < arb("1e-22")
    weighted = [w*den(s) for w,s in zip(weights,nodes,strict=True)]
    for x in xs:
        reconstructed = sum(w*ball(x-1)/(1+ball(x-1)*s) for s,w in zip(nodes,weighted))
        assert reconstructed.overlaps(rational_log(x)) and abs(reconstructed-rational_log(x))<arb("1e-20")
    return xs, nodes, weighted


def verify():
    shape = run_path(str(Path(__file__).with_name("ch39-heavy-family-shape.py")))
    mat, poles, q = (shape[k] for k in ("mat","poles","q"))
    time = poles["study"]["nominal"][2]
    certificates = []
    with ctx.workprec(400):
        for count in (2,4,8):
            xs,nodes,weights = gaussian_spectrum(count)
            critical_mass = (1/nodes[-1]-1).sqrt()
            print(count,"upper support limit M/M_N",critical_mass.str(18),
                  "quadrature atoms",len(nodes),flush=True)
            lo,hi = {2:("1.34864","1.34865"),4:(".69923",".69924"),8:(".37034",".37036")}[count]
            assert arb(lo)<critical_mass<arb(hi)
            certificates.append((xs,nodes,weights,critical_mass))
    xs,nodes,weights,_ = certificates[-1]
    trace = mat["charge_trace"]()
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        n1 = poles["reference_readings"](logs)[3]
        e,b,d = n1**2/scale,c/rho,eta/rho
        assert (3*e).sqrt() < masses[0] and max(row[3] for row in certificates)*e.sqrt() < 2*masses[0]
        k0 = 1/rho
        k = 1/(poles["offset"]+poles["study"]["comp"]["model"]["target"](11430,time)/(15*scale))
        h = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
        assert 0 < k < k0
        error = arb(0)
        for x in xs:
            P,V,_ = shape["probe"](rational_log(x),e,z,b,d,masses[0])
            bare = P + mat["hypercharge"]*(h/k0)
            # Use the physical point masses t_j and positive matrix residues.
            shifts = [e*(1/s-1) for s in nodes]
            projected_log = sum(ball(x-1)*e*(e*w/s**2)/((t+e)*(t+ball(x)*e))
                                for t,s,w in zip(shifts,nodes,weights,strict=True))
            triplet = V[0,0]-trace[0,0]*rational_log(x)/(16*arb.pi()**2)
            assert triplet > 0
            added = trace*(projected_log/(16*arb.pi()**2))
            added[0,0] += triplet
            added *= k0/k-1
            primary = bare.inv()-V*k0
            alternate = bare.inv()-(V+added)*k
            assert mat["positive"](primary) and mat["positive"](alternate)
            diff = primary.inv()-alternate.inv()
            error = max(error,max(abs(v).upper() for v in diff.entries()))
        assert error < q("1e-10")
        print("case",case,"eight full matrix residuals bounded by",error.str(12),flush=True)


if __name__ == "__main__":
    verify()
