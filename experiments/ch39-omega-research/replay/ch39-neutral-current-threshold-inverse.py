"""Recover both coloured masses from three Euclidean family-current readings."""
from pathlib import Path
from runpy import run_path

from flint import arb

study = run_path(str(Path(__file__).with_name("ch39-top-threshold-resolution.py")))
g,u,w,m3,_,active,_,_,_ = study["ports"]
q,roots = study["q"],study["enclose_roots"]
s,t = q("1/512"),q("1/160")
left,right = active(-s),active(-t)
beta_limit = g*g/(56*(q("2/21")).exp())


def samples():
    c,eta,beta = 18*u*u,1/g**2,1/(2*w*w)
    return [1/c,1/(c+s*eta+4/(beta+left)),1/(c+t*eta+4/(beta+right))]


def recover(readings):
    c,first,second = [1/v for v in readings]
    diff = (first-c)/s-(second-c)/t
    def gap(beta):
        return 4/(s*(beta+left))-4/(t*(beta+right))-diff
    assert t.sqrt()*(beta_limit+right) < s.sqrt()*(beta_limit+left)
    beta=roots(gap,0,beta_limit)
    eta=(second-c-4/(beta+right))/t
    k=1/eta.sqrt()
    d=((1-21*(g*g/(56*beta)).log()/2)/5).exp()
    b_squared=g*g/(56*beta)
    def equation(h):
        a2=((1+4*d.log()+h.log()+arb(15).log()/2)/10).exp()
        return 15*a2-14*b_squared-(d/h)**2
    h=roots(equation,".5","1.0")
    triplet=h*g*(c/18).sqrt()
    return (triplet,(triplet*triplet+d*d).sqrt(),k),(c,beta,eta)


def verify():
    nominal=study["nominal"]
    truth=[(v-nominal[2]).exp() for v in nominal[:2]]+[g]
    data=samples()
    exact,_=recover(data)
    for value,expected in zip(exact,truth,strict=True):
        assert value.overlaps(expected) and abs(value-expected)<q("1e-12")
    error=arb(0,q("2e-8").upper())
    result,(c,beta,eta)=recover([v*(1+error) for v in data])
    logs=[abs((v/v0).log()).upper() for v,v0 in zip(result,truth,strict=True)]
    assert all(x<q(b) for x,b in zip(logs,(".000718",".000193",".000009"),strict=True))
    triplet,octet,k=result
    exclude_lower_tops(triplet,octet,c,beta,eta,g*g,nominal[2])
    print("three current readings",*(v.str(12) for v in data))
    print("20 ppb current windows: log errors in M3, M8, g_X",*(v.str(10) for v in logs))


def first_window(triplet,octet,c,beta,eta,g2):
    # Insert the inferred coupling into the physical mass equation, allowing
    # it to vary independently of the site coupling.
    unit,colour,right_green=study["ports"][-3:]
    def response(e,derivative=False):
        return g2*(3*colour(e,unit,derivative)/8+right_green(e,unit,derivative)/4)
    energy=roots(lambda e:beta+response(e)+4/(c-eta*e),".0001",".003")
    assert c/eta>triplet**2 and unit/64>triplet**2 and energy<triplet**2
    norm=(response(energy,True)+4*eta/(c-eta*energy)**2).sqrt()
    a=-g2*colour(energy,unit)/(8*norm)
    b=g2*right_green(energy,unit)/(4*norm)
    f=-2/((c-eta*energy)*norm)
    prefactor=(triplet**2/energy).log()/(16*arb.pi()**2)**2
    weak=q("15/16")*prefactor*(24*f*f+36*a*a+2*b*b)
    pair=prefactor*(q("2.5")*f*f+q("23.75")*a*a+30*a*b+q("9.375")*b*b)
    colour_shift=prefactor*(24*f*f+12*a*a+6*b*b)
    assert max(weak,weak+pair,colour_shift)/(triplet/energy.sqrt()).log()<q(".001")
    d2=octet**2-triplet**2
    h2=triplet**2/(g2*c/18)
    ww=g2/(2*beta)
    assert 1+(d2/h2+ww/2)*colour(triplet**2,unit)<0
    assert 1+ww/2*right_green(triplet**2,unit)<0
    return energy,(a,b,f),weak,pair


def exclude_lower_tops(triplet,octet,c,beta,eta,g2,log_scale):
    _,_,weak,pair = first_window(triplet,octet,c,beta,eta,g2)
    for index,offset in enumerate((".00000266",".00000212")):
        time=study["nominal"][2]+q(offset)
        h3=log_scale+triplet.log()
        h8=log_scale+octet.log()
        assert h3 < time
        lo=study["at_scale"](index,False,h3.lower(),h8.lower(),time)
        hi=study["at_scale"](index,True,h3.upper(),h8.upper(),time)
        weak_gap=study["target_at"](time)-hi[1]+weak
        pair_gap=lo[1]-hi[0]+pair
        assert weak_gap>arb(2)**-21 and pair_gap>arb(2)**-20
        print("excluded group",index,"weak and pair gaps",weak_gap.lower().str(8),pair_gap.lower().str(8))


if __name__=="__main__":
    verify()
