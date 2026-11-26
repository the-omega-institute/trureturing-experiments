"""Exclude the final eight-mode allocation by a CKM row relation and the charm line."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb_mat, arb

five = run_path(str(Path(__file__).with_name("ch39-shell-gw-five-up-three-down.py")))
m = five["m"]
four = m["four"]
adj, hs = m["adj"], m["hs"]
a = arb


def texture_bounds():
    raw = m["images"](5)
    mats = m["orb"]["orthonormalize"](raw)
    y0,y1 = map(four["float_matrix"],raw[:2])
    L = four["unitary_trial"](np.linalg.eig(y1 @ np.linalg.inv(y0))[1])
    R = four["unitary_trial"](y0.conj().T @ L)
    L,R = four["exact_frame"](L),four["exact_frame"](R)
    mats = [adj(L)*b*R for b in mats]
    def row(i,j,basis):
        return acb_mat([[b[i,j] for b in basis]])
    def stack(r,s):
        return acb_mat([[v[0,k] for k in range(v.ncols())] for v in (r,s)])
    D = acb_mat([[b[i,i] for b in mats[:3]] for i in range(3)])
    Dh = acb_mat([[b[i,i] for b in mats[3:]] for i in range(3)])
    def low(i,j):
        return row(i,j,mats[:3])*D.inv()
    def high(i,j):
        return row(i,j,mats[3:])-low(i,j)*Dh
    for p in range(3):
        i,j = [k for k in range(3) if k!=p]
        for s,t in ((i,j),(j,i)):
            l = low(s,t)
            assert abs(l[0,p]) < a(".000060") and max(abs(z) for z in l.entries()) < a(".00029")
        def transfer(i,j):
            col = stack(high(i,p),high(j,p))
            other = stack(high(i,j),high(j,i))
            T = col*other.inv()
            ell = stack(low(i,p),low(j,p))-T*stack(low(i,j),low(j,i))
            t = T[0,1]/T[1,1]
            b = T[0,0]-t*T[1,0]
            f = acb_mat([[ell[0,k]-t*ell[1,k] for k in range(3)]])
            return abs(t),abs(b),hs(f,f).real.sqrt()
        t,b,f = transfer(i,j)
        if t > 1:
            t,b,f = transfer(j,i)
        assert a(".183") < t < a(".185") and b < a(".185") and f < a(".00036")


def exclude():
    # Normalize each sector by its own top singular value.
    tol = a("1e-6")
    mult = (1+tol)/(1-tol)
    ru, rd, light, first = a(".00362"),a(".01914"),a(".0000072"),a(".00097")
    assert light > a(".00000615")/a(".8616")*mult
    assert a(".00311")/a(".8616")/mult > a(".00360")
    assert a(".0000135")/a(".01401")*mult < first
    d = a(".00042")
    top = a(2).sqrt()*d/(1-ru-4*d)
    assert top < a(".00061")
    # Off-diagonal entries of the residual light 2-by-2 block after deflation.
    minor_entry = a(".000060")+2*a(".00029")*(ru+2*d)+top**2
    cross = (minor_entry+light)/a(".00360")
    # The charm left vector has an entry >=sqrt((1-top²)/2) off the top axis.
    # The small entries in the opposite row and column then constrain both sides.
    large = ((1-top**2)/2).sqrt()
    right_minor = cross/large
    left_minor = cross/(1-top**2-right_minor**2).sqrt()
    charm = a(".0203")
    assert left_minor < charm
    first_on_charm = (charm**2+top**2).sqrt()
    top_mix = a(".044")
    theta23,theta13,theta12 = (a(t)+a(0,tol) for t in (".04275",".00377",".22704"))
    vcb,vub = theta23.sin()*theta13.cos(),theta13.sin()
    mix = (1-(theta23.cos()*theta13.cos())**2).sqrt()
    assert mix+top < top_mix
    small = vub+charm*vcb+top
    big = (1-charm**2-top**2).sqrt()*vcb-first_on_charm*vub-top
    strange_light = theta12.sin()*theta13.cos()+charm+top
    assert small < a(".00526") and big > a(".042") and strange_light < a(".247")
    small,big = a(".00526"),a(".042")
    K = five["kappa"]
    fixed = K*((rd**2+first**2).sqrt()+a(2).sqrt()*d)+a(2).sqrt()*d
    beta = (fixed+rd*top_mix)/((1-top_mix**2).sqrt()-K*top_mix)
    assert beta < a(".074")
    beta = a(".074")
    correction = a(".00036")*(1+rd**2+first**2).sqrt()
    for light_row in (True,False):
        r_i = ((rd*a(".247"))**2+first**2).sqrt() if light_row else rd
        bi = small if light_row else top_mix
        lhs = (a(".183")*big-small if light_row else big-a(".185")*small)*(1-beta**2).sqrt()
        rhs = a(".185")*(bi*beta+r_i)+(r_i+a(".185")*rd)*beta+correction
        assert lhs > rhs
        print("light row" if light_row else "charm row",
              "required mismatch >",lhs.lower().str(8),"available <",rhs.upper().str(8))


if __name__ == "__main__":
    five["certify_shapes"]()
    texture_bounds()
    exclude()
