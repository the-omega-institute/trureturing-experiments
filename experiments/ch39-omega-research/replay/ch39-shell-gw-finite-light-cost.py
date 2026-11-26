"""An interval upper bound on the light-mass cost at the measured heavy-family angle."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

tangent = run_path(str(Path(__file__).with_name("ch39-shell-gw-infinitesimal-mixing.py")))
model, rank, seven = tangent["model"], tangent["rank"], tangent["seven"]
adj, hs = model["adj"], model["hs"]
array = model["array"]


def directions(modes,u,v):
    B=np.array(list(map(array,modes)))
    x,y=array(u),array(v)
    P,Q=np.eye(3)-x@x.conj().T,np.eye(3)-y@y.conj().T
    gram_rows=np.vstack([(P@B@Q).reshape(5,9).T,(x.conj().T@B@y).reshape(1,5)])
    W=(P@B@y).reshape(5,3).T
    dual=np.linalg.solve(gram_rows.conj().T@gram_rows,W.conj().T)
    eigs,frame=np.linalg.eigh(W@dual)
    # Fixing any binary trial coefficients produces an exact direction in the band.
    base=dual@frame[:,-1]/np.sqrt(eigs[-1])
    for phase in np.exp(2j*np.pi*np.arange(16)/16):
        yield phase*base


def find_trial(index):
    modes,_=seven["spectral_band"](index)
    left,right=rank["mix"]["orbit_frame"](modes)
    modes=[adj(left)*b*right for b in modes[:5]]
    outside=model["complement"](modes)
    target=float((1-(arb(".04275").cos()*arb(".00377").cos())**2).sqrt())
    best=None
    for seed in rank["numeric_seeds"](outside,0):
        x,v=rank["certified_root"](outside,0,seed)
        u=x.conjugate()
        y=u*adj(v)
        for coords in directions(modes,u,v):
            H=sum((acb(float(z.real),float(z.imag))*b for z,b in zip(coords,modes)),acb_mat(3,3))
            A,base,axis=array(H),array(y),array(u)
            def reading(t):
                U,s,_=np.linalg.svd(base+t*A)
                return np.linalg.norm((np.eye(3)-axis@axis.conj().T)@U[:,0]),np.linalg.norm(s[1:])/s[0]
            lo,hi=0,.025
            for _ in range(32):
                mid=(lo+hi)/2
                if reading(mid)[0]>target:
                    hi=mid
                else:
                    lo=mid
            cost=reading(hi)[1]
            if best is None or cost<best[0]:
                best=(cost,hi,u,y,H)
    return best[1:]


def bounds(Y,u):
    T=Y*adj(Y)
    lam=float(np.linalg.eigvalsh(array(T))[-1])
    I=model["eye"](3)
    lo,hi=arb(lam)-arb("1e-8"),arb(lam)+arb("1e-8")
    assert seven["count"](seven["rows"](T),lo)==2 and seven["count"](seven["rows"](T),hi)==3
    largest=lo.union(hi)
    trace=sum(T[i,i].real for i in range(3))
    small_sum=trace-largest
    small_prod=T.det().real/largest
    top=(T*T-small_sum*T+small_prod*I)/(largest**2-small_sum*largest+small_prod)
    distance=(1-(adj(u)*top*u)[0,0].real).sqrt()
    return distance,(small_sum/largest).sqrt()


if __name__=="__main__":
    target=(1-(arb(".04275").cos()*arb(".00377").cos())**2).sqrt()
    for n in (43,195):
        t,u,y,H=find_trial(n)
        # Endpoint distances cross the exact target. Use singular-value Lipschitz
        # bounds from the center to control the whole bracket's cost.
        width=arb(".000003")
        left,_=bounds(y+(arb(t)-width)*H,u)
        right,_=bounds(y+(arb(t)+width)*H,u)
        assert left < target < right
        center=y+arb(t)*H
        _,tau=bounds(center,u)
        normH=hs(H,H).real.sqrt()
        floor=abs((adj(u)*center*adj(y)*u)[0,0])
        eta=width*normH/floor
        upper=(tau+eta)/(1-eta)
        assert upper < arb(".0123" if n==43 else ".0107")
        print("zero",n,"attained exact target",target.str(10),"light cost below",upper.upper().str(9))
