"""Exclude nonnegative integer powers of the original spin overlap as two contact strengths."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

single = run_path(str(Path(__file__).with_name("ch39-single-contact-ckm-box.py")))
arc = single["arc"]
projectors = arc["arc_projectors"]()
q32,q152 = projectors[32],projectors[152]
base = arb(829)/859


def positive_check(index,contact):
    native = arc["native"]
    basis = arc["basis"]
    t = native["plane_response"](acb.zeta_zero(index).imag)
    tx,ty,tz = arc["components"](t)
    bx,_,bz = arc["components"](contact)
    def admissible(w):
        u = basis.transpose()*w*basis
        kx = (u[0,1]+u[1,0])/2
        # Radial fixation represents the skew matrix by cross product with p.
        # The two small CKM entries keep this enclosure tight at finite width.
        ky = (w[0,2]-w[2,0])*arb(867).sqrt()/(10*acb(0,1))
        kz = (u[0,0]-u[1,1])/2
        sign = ((ty*kx-tx*ky)*bx+(ty*kz-tz*ky)*bz)*ky.conjugate()
        return not sign.real < 0
    return admissible


def verify():
    # h152/h32=(kappa152/kappa32)^m for kappa_d=1-d/3468.
    assert base.overlaps((1-arb(152)/3468)/(1-arb(32)/3468))
    native = arc["native"]
    for index,radius in ((11431,".003"),(11439,".10"),(11737,".20")):
        t = native["plane_response"](acb.zeta_zero(index).imag)
        basis = arc["basis"]
        n,m = [basis*(t*q-q*t)*basis.transpose() for q in (q32,q152)]
        count = single["exclude_normals"]
        # Exclude a symmetric interval around unit strength ratio. Partition
        # the interval to retain correlations when forming the commutator.
        lo,hi = 1-arb(radius),1+arb(radius)
        near = [lo+(hi-lo)*j/128 for j in range(129)]
        count([n+u.union(v)*m for u,v in zip(near,near[1:])])
        for m_exp in range(46):
            b = q32+(base**m_exp)*q152
            sign = positive_check(index,b) if (index,m_exp)==(11737,32) else None
            count([n+(base**m_exp)*m],admissible=sign)
        # One continuous family also covers arbitrarily large integer m.
        high = base**46
        ends = [high*j/64 for j in range(65)]
        count([n+u.union(v)*m for u,v in zip(ends,ends[1:])])
        print(index,"equal-strength radius",radius,"and every integer overlap power excluded")
        if index in (11439,11737):
            for j in range(64):
                lo,hi = arb(j)/64,arb(j+1)/64
                b = q32+lo.union(hi)*q152
                ends = [lo+(hi-lo)*k/16 for k in range(17)]
                normals = [n+u.union(v)*m for u,v in zip(ends,ends[1:])]
                count(normals,admissible=positive_check(index,b))
            print(index,"all nonnegative, nondecreasing overlap laws excluded")


if __name__ == "__main__":
    verify()
