"""First neutral mass window, with all light SM currents in the two-loop gauge logarithm."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
neutral = run_path(str(root/"ch39-shared-quartic-neutral.py"))

# Charges (D, 2T_R, family sign) for left-handed Weyl fields.
weyl_charges = {"Q":(1,0,1), "uc":(-1,-1,-1), "dc":(-1,1,-1),
                "L":(-3,0,1), "ec":(3,1,-1)}


def light_current(logs):
    g,u,w,m3,h,active,unit,gc,gr = neutral["parameters"](logs,include_graphs=True)
    start = neutral["bisect"](lambda e:1/(2*w*w)+active(e)+4*g*g/(18*g*g*u*u-e),
                              ".0001",".003")
    c=18*u*u
    norm = (active(start,True)+4*g*g/(c*g*g-start)**2).sqrt()
    alpha=-g*g*gc(start,unit)/(8*norm)
    beta=g*g*gr(start,unit)/(4*norm)
    family=-2*g*g/((c*g*g-start)*norm)

    # The charged Higgsed towers also have their first eigenvalue above m3.
    b=w*g/arb(28).sqrt()
    q = neutral["graph"]["squared"]
    points=[p for p in product(range(30),repeat=3) if sum(y*y for y in p)==867]
    left=[p for p in points if tuple(sorted(p)) in ((1,5,29),(7,17,23),(11,11,25),(17,17,17))]
    def left_edge(p,r):
        return q(p,r)!=136 or (11,11,25) in (tuple(sorted(p)),tuple(sorted(r)))
    neutral["pole_resolvent"](left,(32,152,56,136),27,edge_filter=left_edge)
    assert c*g*g>m3*m3 and unit/64>m3*m3
    spectrum=run_path(str(root/"ch39-published-split-adjoint.py"))
    data=spectrum["data"]
    exact=lambda v:data["q"](v)
    d_lower=((exact(logs[1])-exact(logs[0])).exp()**2-1).sqrt()*m3
    assert 1+(d_lower*d_lower/(h*h)+14*b*b)*gc(m3*m3,unit)<0
    assert 1+14*b*b*gr(m3*m3,unit)<0

    shifts = [v*(m3/start.sqrt()).log() for v in gauge_rates(alpha,beta,family)]
    scaled=[value/r for value,r in zip(shifts,(arb(24)/25,arb(16)/15,arb(1)))]
    return start,scaled,(alpha,beta,family),shifts


def gauge_rates(alpha,beta,family):
    # Y uses the unified normalization. Fields are Q, uc, dc, L, ec.
    indices = [
        (F(1,10),F(4,5),F(1,5),F(3,10),F(3,5)),
        (F(3,2),0,0,F(1,2),0),
        (1,F(1,2),F(1,2),0,0),
    ]
    shifts=[]
    factor = 2/(16*arb.pi()**2)**2
    for row,higgs in zip(indices,(F(3,10),F(1,2),0)):
        total=arb(0)
        for x in (1,1,-2):
            for t,(dc,r,sign) in zip(row,weyl_charges.values(),strict=True):
                total+=2*arb(t.numerator)/t.denominator*(dc*alpha+r*beta+sign*x*family)**2
        total+=4*arb(higgs.numerator)/higgs.denominator*beta**2
        shifts.append(total*factor)
    return shifts


def verify():
    _,scaled,(alpha,beta,family),shifts=light_current(neutral["primary_logs"])
    assert abs(abs(family)-arb(".2027318524"))<arb("1e-10")
    spectrum=run_path(str(root/"ch39-published-split-adjoint.py"))
    data=spectrum["data"]
    exact=data["q"]
    difference=(max(scaled)-min(scaled))/2
    assert arb("2.62164e-6") < difference < arb("2.62165e-6")
    # Add just this finite-window logarithm to the original trajectory.
    # Its magnitude ramps up monotonically between the two masses.
    comp=spectrum["comp"]
    loop=comp["loop"]
    upper=data["corner"](1,True)
    step=arb(33)/256
    for _ in range(256):
        upper,tube=comp["advance"](upper,step,True,need_tube=True)
        assert all(sum(tube[i:i+3])<3 for i in (3,6,9))
    rate_g=arb(0).union(arb(2))
    rate_y=arb(0).union(arb(3))
    for three,eight in ((False,False),(True,False),(True,True)):
        flow=spectrum["flow"](three,eight)
        speeds=[-2/(loop*r)*(b+sum(z*rate_g-t*rate_y for z,t in zip(row,d))/loop)
                for r,b,row,d in zip(spectrum["ratios"],flow["BETA"],flow["B"],flow["D"])]
        assert speeds[1]-speeds[0]>arb(".08") and speeds[2]-speeds[1]>arb(".031")
    raw_tol=arb("9.08e-8")
    trial=exact("26.81725683017028")
    lower_mass=exact("24.3145558253697")
    early=(trial-lower_mass)*arb(".031")
    left=scaled[2]-scaled[1]
    right=scaled[0]-scaled[1]
    radius=(min(left,right,early)-2*raw_tol)/2
    assert radius>arb("8.07e-7")>arb(2)**-21
    print("isolated logarithm: every common target and scale has residual above",radius.lower().str(9))
    print("site neutral couplings D, 2T3R, X",*(v.str(12) for v in (alpha,beta,family)))
    print("inverse gauge shifts (magnitudes)",*(v.str(12) for v in shifts))
    print("scaled shifts",*(v.str(12) for v in scaled),"minimax after common shift",difference.str(12))


if __name__=="__main__":
    verify()
