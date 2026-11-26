"""Identify the arithmetic target from two visible charged poles and one quark partial width."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_poly

root = Path(__file__).parent
cross = run_path(str(root/"ch39-retuned-data-top-identification.py"))
retune, window, study, q = (cross[k] for k in ("retune","window","study","q"))
four = cross["fit"]["four"]
neutral = study["retuned"]["current"]["neutral"]
roots = study["enclose_roots"]
unit, common, colour, right = neutral["graph_constants"]()
offset = common/(48*arb.pi()**2)
bands = ((".025",".032"),(".35",".385"))


def charged_poles(z):
    return tuple(roots(lambda e: right(e,unit)+4*z, a,b) for a,b in bands)


def verify_ratio_monotonicity():
    _,_,_,(p,v) = neutral["graph_constants"](port_builder=neutral["pole_polynomials"])
    spectrum = []
    for value,multiplicity in v.complex_roots():
        assert value.imag == 0 and multiplicity == 1
        r = value.real
        weight = p(r)/v.derivative()(r)
        assert weight > 0
        spectrum.append((r*unit,weight))
    spectrum.sort(key=lambda pair: float(pair[0].mid()))
    assert spectrum[0][0] == 0 and q(".032")<spectrum[1][0]<q(".35")
    assert q(".385") < spectrum[2][0]
    bounds = []
    for lo,hi in zip(charged_poles(four["z_limit"]),charged_poles(arb(0)),strict=True):
        e = lo.union(hi)
        derivative = sum(w/(r-e)**2 for r,w in spectrum)
        bounds.append(1/(e*derivative))
    assert bounds[0] > bounds[1]
    return (4*(bounds[0]-bounds[1])).lower()


def reference_readings(logs):
    time = study["nominal"][2]
    g,_,w,_,_,_ = neutral["parameters"](logs)
    L = (2*(q(logs[2])-time)).exp()
    z = 1/(2*w*w*g*g)
    # Diagonalize the rank-one perturbation's reduced characteristic
    # polynomial independently of the secular bisection used to invert it.
    _,_,_,(p,v) = neutral["graph_constants"](port_builder=neutral["pole_polynomials"])
    polynomial = arb_poly(v) - arb_poly(p)*(1/(4*z*unit))
    eigenvalues = sorted((x.real*unit for x in polynomial.complex_roots(tol=1e-30)
                          if x.imag.contains(0)),key=lambda x: float(x.mid()))
    assert len(eigenvalues) == v.degree() and min(eigenvalues) > 0
    first,second = eigenvalues[:2]
    component = 16*z**2/right(first,unit,True)
    width_r = 3*g*g*component/(16*arb.pi())
    en,_,couplings,_ = study["retuned"]["current"]["light_current"](logs)
    width_n = cross["fit"]["widths"]["reference_width"](couplings)
    return [(L*e).sqrt() for e in (first,second)]+[width_r,(L*en).sqrt(),width_n]


def charged_geometry(mass1,mass2):
    observed = (mass2/mass1)**2
    def gap(z):
        e1,e2 = charged_poles(z)
        return e2/e1-observed
    z = roots(gap,0,four["z_limit"])
    e,_ = charged_poles(z)
    L = mass1**2/e
    return L,z,e


def recover_charged(mass1,mass2,width):
    L,z,e = charged_geometry(mass1,mass2)
    weight = 16*z*z/right(e,unit,True)
    rho = 3*weight/(16*arb.pi()*width)
    return 15*L*(rho-offset),L,z,rho


def recover_neutral(mass,width,L,z,rho):
    e = mass**2/L
    a = 3*colour(e,unit)/8+right(e,unit)/4
    derivative = 3*colour(e,unit,True)/8+right(e,unit,True)/4
    delta = -4*rho/(z+a)
    rD = colour(e,unit)*delta/(16*rho)
    rR = -right(e,unit)*delta/(8*rho)
    eta = (12+27*rD**2+3*(3*rD+rR)**2)/(24*arb.pi()*width)-derivative*delta**2/(4*rho)
    c = delta+eta*e
    assert min(eta,c,delta)>0 and e<unit/8
    return eta,c,four["physical_masses"](c,rho,z)


def verify():
    print("log mass-squared ratio derivative lower bound",verify_ratio_monotonicity().str(10))
    time = study["nominal"][2]
    targets = [study["comp"]["model"]["target"](n,time) for n in (11430,11431,11432)]
    graph_gap = window["graph_gap"]()
    for i,logs in enumerate(retune["log_sets"]):
        readings = reference_readings(logs)
        expected,_,_ = retune["predictions"](logs,time)
        value,L,z,rho = recover_charged(*readings[:3])
        recovered = (L,z,rho,*recover_neutral(*readings[3:],L,z,rho))
        for a,b in zip(recovered[:5]+recovered[5],expected[:5]+expected[5],strict=True):
            assert a.overlaps(b) and abs(a-b) < q("1e-11")
        assert value.overlaps(targets[1])
        data = [v*(1+arb(0,q(error).upper())) for v,error in
                zip(readings,("1e-8","1e-8","1e-6","1e-8","1e-6"),strict=True)]
        target,L,z,rho = recover_charged(*data[:3])
        assert targets[0]>target>targets[2]
        rho = rho.intersection(targets[1]/(15*L)+offset)
        eta,c,masses = recover_neutral(*data[3:],L,z,rho)
        print("case",i,"reference readings",*(v.str(13) for v in readings))
        print("target interval",target.lower().str(14),target.upper().str(14))
        print("absolute mass log errors",
              *(abs((m*L.sqrt()/(r*expected[0].sqrt())).log()).upper().str(9)
                for m,r in zip(masses,expected[-1],strict=True)))
        windows = window["window_data"]((L,z,rho,eta,c,masses),graph_gap)
        cross["exclude_other_orbits"](i,windows,q(logs[2]))


if __name__=="__main__":
    verify()
