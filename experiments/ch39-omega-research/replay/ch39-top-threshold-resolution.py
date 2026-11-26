"""Certify the mass accuracy needed to distinguish the three discrete top inputs."""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import arb, arb_series

root = Path(__file__).parent
tables = run_path(str(root/"ch39-top-table-selection.py"))
retuned, empirical, q = tables["retuned"], tables["empirical"], tables["q"]
adjoint, comp = empirical["adjoint"], empirical["comp"]
nominal = tuple(map(q,retuned["logs"]))
reference = comp["model"]["target"](11431,nominal[2])
ports = retuned["current"]["neutral"]["parameters"](retuned["logs"],include_graphs=True)
site_coupling = ports[0]
unit, colour, right = ports[-3:]


def enclose_roots(function,lo,hi):
    lo,hi = arb(lo),arb(hi)
    assert function(lo)<0<function(hi)
    endpoints = []
    for lower_endpoint in (True,False):
        a,b = lo,hi
        for _ in range(55):
            mid = ((a+b)/2).mid()
            sign = function(mid)
            if sign<0 or (not lower_endpoint and not sign>0):
                a = mid
            else:
                b = mid
        endpoints.append(a if lower_endpoint else b)
    return endpoints[0].union(endpoints[1])


def jet(value,slope=0):
    return arb_series([value,slope],prec=2)


def shift_derivatives(h3,h8,t):
    g2 = 1/(1/site_coupling**2 + reference*((-2*(t-nominal[2])).exp()-1)/15)
    m3,m8 = (h3-t).exp(),(h8-t).exp()
    d = (m8*m8-m3*m3).sqrt()
    b = ((1-5*d.log())/21).exp()

    def a_squared(h):
        return ((1+4*d.log()+h.log()+arb(15).log()/2)/10).exp()

    h = enclose_roots(lambda z:15*a_squared(z)-14*b*b-(d/z)**2,".5","1.0")
    w2,c = 28*b*b,18*(m3/h)**2
    def secular(e):
        return 1/(2*w2)+3*colour(e,unit)/8+right(e,unit)/4+4/(c-e)
    assert c>arb(".003") and unit/64>m3*m3
    energy = enclose_roots(secular,".0001",".003")
    assert (m3/energy.sqrt()).log()>q(".1")
    assert 1+(d*d/(h*h)+w2/2)*colour(m3*m3,unit)<0
    assert 1+w2/2*right(m3*m3,unit)<0
    gc,gr = [tuple(f(energy,unit,i) for i in range(3)) for f in (colour,right)]
    norm2 = 3*gc[1]/8+gr[1]/4+4/(c-energy)**2
    assert norm2>0
    slopes = []
    for mode in range(3):
        dt = int(mode==2)
        d3,d8 = int(mode==0)-dt,int(mode==1)-dt
        dlog = (m8*m8*d8-m3*m3*d3)/(d*d)
        hlog = -dlog*(6*a_squared(h)+arb(20)/3*b*b-2*(d/h)**2)/(3*a_squared(h)/2+2*(d/h)**2)
        cdot = 2*c*(d3-hlog)
        wlog = -10*dlog/21
        edot = (wlog/(2*w2)+4*cdot/(c-energy)**2)/norm2
        diff = jet(c-energy,cdot-edot)
        norm = (3*jet(gc[1],gc[2]*edot)/8 + jet(gr[1],gr[2]*edot)/4 + 4/diff**2).sqrt()
        alpha = -jet(gc[0],gc[1]*edot)/(8*norm)
        beta = jet(gr[0],gr[1]*edot)/(4*norm)
        family = -2/(diff*norm)
        gt = target_at(t)
        coupling = jet(g2,2*g2*g2*gt*dt/15)
        factor = coupling*(jet(m3,m3*d3)**2/jet(energy,edot)).log()/(16*arb.pi()**2)**2
        weak = (24*family**2+36*alpha**2+2*beta**2)*factor*q("15/16")
        one = (q("2.5")*family**2+q("23.75")*alpha**2+30*alpha*beta+q("9.375")*beta**2)*factor
        three = (q("1.5")*family**2-q("21.75")*alpha**2+q("4.125")*beta**2)*factor
        slopes.append([weak[1],one[1],three[1]])
    return slopes


def target_at(time):
    return reference*(-2*(time-nominal[2])).exp()


@lru_cache(None)
def square_inputs(index,upper):
    bands = tables["table_bands"]()
    word = empirical["endpoint"]
    state = [q(Fraction(v,word[0])**2) for v in word[1:]]
    for j,row in enumerate(bands):
        state.append(q(row[index if j==2 else (-1 if upper else 0)][0]))
    flow = adjoint["flow"](False,False)
    for _ in range(160):
        state = comp["advance"](state,q(24)/160,upper,flow=flow)
    return state


def at_scale(index,upper,h3,h8,time):
    state = square_inputs(index,upper)
    for a,b,n,mode in ((q(24),h3,4,(False,False)),(h3,h8,28,(True,False)),
                       (h8,time,8,(True,True))):
        flow = adjoint["flow"](*mode)
        for _ in range(n):
            state = comp["advance"](state,(b-a)/n,upper,flow=flow)
    return [1/(r*s) for r,s in zip(adjoint["ratios"],state[:3],strict=True)]


def correction_bounds(radius,clock_radius):
    width=(radius,radius,clock_radius)
    def quarters(center,r):
        return [(center-r+r*i/2).union(center-r+r*(i+1)/2) for i in range(4)]
    first,second,third=[quarters(a,r) for a,r in zip(nominal,width,strict=True)]
    bounds = [[arb(0)]*3 for _ in range(3)]
    for a in first:
        for b in second:
            for c in third:
                slopes=shift_derivatives(a,b,c)
                for k,row in enumerate(slopes):
                    for j,v in enumerate(row):
                        bounds[k][j]=bounds[k][j].union(v)
    _,d,_,_=retuned["current"]["light_current"](retuned["logs"])
    centres=(d[1],d[0]-d[1],d[2]-d[1])
    return [v+arb(0,sum(abs(bounds[k][j])*width[k] for k in range(3)).upper())
            for j,v in enumerate(centres)]


def verify():
    a,b,t=nominal
    epsilon=arb(2)**-21
    cases=(
        (".00105",".00000266",("24.31743852737601","26.71339683172463","26.81726789380154")),
        (".00082",".00000212",("24.31716884321350","26.71327119000509","26.81726734934438")),
    )
    bands=tables["table_bands"]()
    for index,(radius,dt,witness) in enumerate(cases):
        width=q(radius)
        weak,d12,d32=correction_bounds(width,q(".001"))
        test=t+q(dt)
        lo=at_scale(index,False,a-width,b-width,test)
        hi=at_scale(index,True,a+width,b+width,test)
        weak_gap=target_at(test)-(hi[1]-weak)
        pair_gap=lo[1]-hi[0]+d12
        assert weak_gap>epsilon and pair_gap>2*epsilon
        assert max(weak,weak+d12,weak+d32)<q(".0001")
        print("log radius",radius,"excluded by weak residual and pair gap",
              weak_gap.lower().str(8),pair_gap.lower().str(8))
        h3,h8,pivot=map(q,witness)
        state=[]
        word=empirical["endpoint"]
        for upper in (False,True):
            gauge=[q(abs(Fraction(v,word[0]))) for v in word[1:]]
            state.append(gauge+[q(row[index if j==2 else (-1 if upper else 0)][0]).sqrt()
                                for j,row in enumerate(bands)])
        energy,offsets,_,_=retuned["current"]["light_current"](witness)
        onset=pivot+energy.log()/2
        segments=((arb(0),h3,160,adjoint["flow"](False,False)),
                  (h3,h8,28,adjoint["flow"](True,False)),
                  (h8,pivot,12,adjoint["flow"](True,True)),
                  (pivot,arb(33),64,adjoint["flow"](True,True)))
        empirical["triplet"]["verify_ground"](state[0],segments,2,pivot,q("4.3e-7"),q(".136"),
            upper_values=state[1],logarithm=(onset,h3,offsets),stationary_readout=True)
        assert abs(pivot-t)<q(".001")
        assert max(abs(h3-a),abs(h8-b))<q((".00123420",".00096452")[index])
        print("all 80 recovered at mass log radius",max(abs(h3-a),abs(h8-b)).str(12))


if __name__=="__main__":
    verify()
