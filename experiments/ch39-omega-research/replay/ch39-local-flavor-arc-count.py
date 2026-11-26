"""Real nonnegative contact terms on incident shell arcs at three specified zeros."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

native = run_path(str(Path(__file__).with_name("ch39-native-flavor-eigenaxes.py")))
basis, fit, p = (native[k] for k in ("basis","fit","p"))
shell = fit["shell"]


def components(h):
    return (h[0,1].real, -h[0,1].imag, (h[0,0].real-h[1,1].real)/2)


def arc_projectors():
    arcs = {}
    for point in native["geo"]["sites"]:
        distance = sum((a-b)**2 for a,b in zip(point,p))
        if 0 < distance <= 804:
            n = shell["cross"](p,point)
            col = basis.transpose()*acb_mat([[v/native["norm"](n)] for v in n])
            assert distance not in arcs
            arcs[distance] = col*col.transpose()
    assert sorted(arcs) == [32,152,216,244,296,504,544,632,724,804]
    return arcs


def central_actions():
    actions = []
    for xy in fit["enclose_roots"]():
        w = basis.transpose()*acb_mat(fit["rephased_matrix"](xy))*basis
        tr = w[0,0]+w[1,1]
        hd = (abs(tr)/2).acos()
        centered = w*(tr.conjugate()/abs(tr))
        traceless = (centered-centered.conjugate().transpose())*(hd/(2*acb(0,1)*hd.sin()))
        actions.append((traceless,tr.arg()))
    return actions


def verify(index, outer, inner, expected_gain, actions, arcs):
    t = native["plane_response"](acb.zeta_zero(index).imag)
    tv = components(t)
    gains = []
    for k, angle in actions:
        kv = components(k)
        for projector in arcs.values():
            pv = components(projector)
            triple = shell["dot"](shell["cross"](tv,pv),kv)
            assert 0 not in triple
        gain = kv[1]/tv[1]
        if gain > 0:
            gains.append((gain,k,angle))
    gain,k,angle = min(gains,key=lambda x:float(x[0].mid()))
    # The signed logarithm-branch certificate gives the absolute minimum over
    # every phase winding (including the nearest opposite-signed branch).
    assert abs(gain-arb(expected_gain)) < arb("1e-8")
    assert all(gain < v for v,_,_ in gains if v is not gain)
    rem = k-t*gain
    a,_,b = components(rem)
    def solve_pair(d,e):
        x,_,z = components(arcs[d])
        y,_,w = components(arcs[e])
        pair = arb_mat([[x,y],[z,w]]).solve(arb_mat([[a],[b]]))
        return pair[0,0],pair[1,0]
    if index == 11737:
        low = solve_pair(32,152)
        assert any(v < 0 for v in low)
    h,j = solve_pair(outer,inner)
    assert h > 0 and j > 0
    iso = angle+((rem[0,0]+rem[1,1]).real-h-j)/2
    iso -= 2*arb.pi()*(iso/(2*arb.pi())).floor()
    assert 0 < iso and iso < 2*arb.pi()
    print(index,"arc lengths squared",outer,inner,"strengths",h.str(10),j.str(10),
          "isotropic phase",iso.str(10))


if __name__ == "__main__":
    arcs, actions = arc_projectors(), central_actions()
    for case in ((11431,32,152,".72539605969"),(11439,32,152,"3.18344268088"),
                 (11737,216,152,"10.5835553016")):
        verify(*case,actions,arcs)
