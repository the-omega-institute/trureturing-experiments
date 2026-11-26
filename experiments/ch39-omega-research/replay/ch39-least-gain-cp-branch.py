"""Prove the global least native gain at zero 11431 using only the three angle centers."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

model = run_path(str(Path(__file__).with_name("ch39-four-contact-channels.py")))
native, contact, identity, basis = (model[k] for k in ("native","contact","identity","basis"))
indices, targets = model["indices"], model["targets"]
radial = [model["radial"][i,j].real for i,j in indices]
eye3 = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])


def traceless_coordinates(matrix):
    reduced = matrix-matrix.trace().real/2*identity
    embedded = basis*reduced*basis.transpose()
    return model["fixed"]["arc"]["components"](reduced),[embedded[i,j] for i,j in indices]


tv,te = traceless_coordinates(native)
bv,be = traceless_coordinates(contact)


def phase_box(gain,height):
    v = [gain*t+height*b for t,b in zip(tv,bv)]
    gap = sum(z*z for z in v).nonnegative_part().sqrt()
    co,si = gap.cos(),gap.sinc()
    z = [-a*co+acb(0,1)*(gain*t+height*b)*si for a,t,b in zip(radial,te,be)]
    if any(not (abs(b)-a).union(abs(b)+a).overlaps(q.sqrt()) for a,b,q in zip(radial,z,targets)):
        return None
    rows = [(b.real,-b.imag,(q-a*a-abs(b)**2)/(2*a)) for a,b,q in zip(radial,z,targets)]
    cs = sn = arb(0,1)
    for i,j in ((0,1),(2,1),(0,2)):
        a,b,v = rows[i]
        k,l,w = rows[j]
        det = a*l-b*k
        if 0 in det:
            continue
        x,y = (v*l-b*w)/det,(a*w-k*v)/det
        if not cs.overlaps(x) or not sn.overlaps(y):
            return None
        cs,sn = cs.intersection(x),sn.intersection(y)
        if 1 not in cs**2+sn**2:
            return None
    if any(0 not in a*cs+b*sn-v for a,b,v in rows):
        return None
    return cs,sn


def neighborhood():
    root = model["center_fit"]((".7210227965417486","2.557578307275763","3.086954654848417"))
    mid = [x.mid() for x in root]
    outer = [v+arb(0,".00001") for v in mid]
    _,jac = model["observables"](mid,True)
    _,tube = model["observables"](outer,True)
    rem = eye3-jac.inv()*tube
    assert max(sum(abs(rem[i,j]) for j in range(3)).upper() for i in range(3)) < 1
    assert all(a.contains_interior(b) for a,b in zip(outer,root))
    return root,outer


def contract(box,outer):
    for _ in range(5):
        if all(a.contains_interior(b) for a,b in zip(outer,box)):
            return True
        mid = [v.mid() for v in box]
        f,jac = model["observables"](mid,True)
        a = jac.inv()
        _,tube = model["observables"](box,True)
        update = -a*f+(eye3-a*tube)*arb_mat([[v-c] for v,c in zip(box,mid)])
        new = [c+update[i,0] for i,c in enumerate(mid)]
        if any(not v.overlaps(n) for v,n in zip(box,new)):
            return True
        box = [v.intersection(n) for v,n in zip(box,new)]
    return False


def verify():
    root,outer = neighborhood()
    ceiling = arb(".722")
    assert root[0] < ceiling
    s,t,u = (model[k] for k in ("s","t","u"))
    modulus_gap = s*t-(1-s*s).sqrt()*(1-t*t).sqrt()*u-u
    assert modulus_gap > 0
    # |W13-W31| is bounded by 2*gain*|Im(T13)|/gap independently of phase.
    gap_max = 2*ceiling*abs(te[2].imag)/modulus_gap
    tnorm,bnorm = [(sum(x*x for x in v)).sqrt() for v in (tv,bv)]
    assert (gap_max+ceiling*tnorm)/bnorm < 16
    pending = [(0,0,0)]
    boxes = 0
    while pending:
        depth,i,j = pending.pop()
        boxes += 1
        gain = ceiling*arb(i).union(arb(i+1))/(1<<depth)
        height = 16*arb(j).union(arb(j+1))/(1<<depth)
        phase = phase_box(gain,height)
        if phase is None:
            continue
        cs,sn = phase
        if depth >= 10 and cs > 0 and sn < 0:
            alpha = acb(cs,sn).arg()-(gain*native.trace().real+height*contact.trace().real)/2+2*arb.pi()
            if contract([gain,height,alpha],outer):
                continue
        assert depth < 20
        pending.extend((depth+1,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    delta,jcp = model["phase_readout"](root)
    assert arb("1.14102143533") < delta < arb("1.14102143535")
    print("global smallest gain and selected delta",root[0].str(15),delta.str(15),
          "J",jcp.str(15),"gain/contact rectangles",boxes)


if __name__ == "__main__":
    verify()
