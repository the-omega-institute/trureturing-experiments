"""Uniform rank-one transfer obstruction on the robust flux-three torus at zero energy."""
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb

work=Path(__file__).resolve().parent
tri=runpy.run_path(str(work/"ch39-schur-triangulation-density.py"))
model=tri["model"]
W,J=model["wall"],model["inside"]
shell=model["shell"]
block=model["block"]
omega=acb(-arb(1)/2,arb(3).sqrt()/2)
c=tri["coefficients"]
modes=[(0,0),(-1,0),(0,-1)]
outer=(11,12,18)
inner=(6,7,8)
rows=[[sum(c[k][W.index(6),W.index(j)]*omega**(l*i) for i,j in enumerate(outer)) for k in modes]
      for l in range(3)]


def bounds_on_defects():
    m=shell["matrix"](0,{})
    inv=block(m,J,J).inv()
    part=[block(m,W,J),acb_mat(7,12),acb_mat(7,12)]
    driven=[]
    for i,j,d,z in shell["edges"]:
        if d in (724,804):
            driven.append((i,j))
            row,col=W.index(j),J.index(i)
            part[1+(d==804)][row,col]=part[0][row,col]
            part[0][row,col]=0
    part=[u*inv for u in part]
    maximum=arb(0)
    for i,j in driven:
        # A clocked edge vector has a unit wall component at its inner-ring endpoint.
        def norm(sites):
            return [
                sum(abs(p[W.index(k),J.index(i)])**2 for k in sites).sqrt() for p in part]
        inner_bound=1+sum(norm(inner))
        outer_bound=sum(norm(outer))
        g=inv[J.index(i),J.index(i)].real
        maximum=max(maximum,(inner_bound*outer_bound/g).upper())
    return maximum


def cover():
    correction=bounds_on_defects()
    assert correction<arb("2.8")
    # A retained small triangle along the outer edge must have v > t_face-pi > -arg(y).
    # Its real part is positive, hence Im(y*loop)>0. These redundant conditions still
    # rule out spurious boxes near zeros of the inner edge X.
    tests=tri["tests"]+[tri["terms"](tri["loop"](f),tri["y"]) for f in [(6,12,11),(7,12,11)]]
    todo=[(0,0,0)]
    accepted=0
    level_max=0
    while todo:
        level,i,j=todo.pop()
        if any(tri["test_value"](t,level,i,j)<0 for t in tests):
            continue
        ph=[1]
        for k in (i,j):
            co,si=tri["trig"](level,2*k+1,1)
            ph.append(acb(co,-si))
        ev=[abs(sum(w*v for w,v in zip(row,ph))) for row in rows]
        # Singular values interlace across one rank-one update; the second keeps
        # at least the smallest undeformed singular value, with arbitrary complex alignment.
        if min(z.lower() for z in ev)>arb(".095")*(max(z.upper() for z in ev)+correction):
            accepted+=1
            level_max=max(level_max,level)
            continue
        assert level<18
        todo.extend((level+1,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    return correction,accepted,level_max


if __name__=="__main__":
    correction,count,depth=cover()
    print("uniform operator perturbation upper, certified clock boxes and maximum depth",
          correction.str(12),count,depth)
    print("bare middle/heavy ratio > 0.095; required endpoint condition products",
          (arb(".095")/(arb(".00316")/arb(".8573"))).str(12),
          (arb(".095")/(arb(".000271")/arb(".01392"))).str(12))
