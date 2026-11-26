"""Stiffness on the twelve actual clock edges, selected by the Gaussian determinant."""
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat

src=runpy.run_path(str(Path(__file__).with_name("ch39-single-arc-source.py")))
s=src["shell"]
J,W=src["J"],src["W"]
block=src["block"]
edges=[(i,j,d) for i,j,d,_ in s["edges"] if d in (724,804)]
p,q=(s[key]*src["ref"] for key in ("a","b"))


def evaluate(p,q):
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    full=m.inv()
    c=block(m,J,J).inv()
    b=block(m,W,J)
    k=block(m,W,W)-b*c*b.conjugate().transpose()
    resistances={}
    for i,j,d in edges:
        phase=-m[i,j]
        resistances[i,j]=(full[i,i]+full[j,j]-phase.conjugate()*full[i,j]-phase*full[j,i]).real
    # h_i=1, h_j=-conj(phase), so the reduced perturbation vector is h_W - B C^-1 h_J.
    i,j=1,6
    v=-b*block(c,range(len(J)),(J.index(i),))
    v[W.index(j),0]+=m[i,j].conjugate()
    return resistances,k,v,c[J.index(i),J.index(i)].real


def check(p,q,radius):
    res,k,v,g=evaluate(p+arb(0,radius),q+arb(0,radius))
    win={(1,6),(2,8),(5,7)}
    assert max(res[a].upper() for a in win)<min(res[a].lower() for a in res.keys()-win)
    eig,a,b=src["fourier_columns"](k,v)
    for f in src["all_faces"]():
        z=src["value"](src["loop_poly"](f,k,v),arb(0,arb(".001")))
        assert z.real>0 or not z.imag.contains(0)
    assert all(abs(z)>arb(".02") for z in a+b)
    lo,hi=src["jarlskog"](eig,a,b,arb(".001"))
    assert lo>0
    return res[1,6],min(res[a].lower() for a in res.keys()-win)-res[1,6],lo,hi


if __name__=="__main__":
    r,gap,lo,hi=check(p,q,arb(".00001"))
    print("winning resistance, gap, Jarlskog enclosure",r.str(8),gap.lower().str(12),
          lo.lower().str(12),hi.upper().str(12))
