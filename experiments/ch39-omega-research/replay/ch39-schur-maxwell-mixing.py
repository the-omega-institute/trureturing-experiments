"""Three-point mixing bound for degree-three Schur faces at all spectral parameters."""
from pathlib import Path
import runpy
from flint import arb,acb_mat

work=Path(__file__).resolve().parent
s=runpy.run_path(str(work/"ch39-shell-schur-flux.py"))
calc=runpy.run_path(str(work/"ch39-wall-certified-moments.py"))
h=calc["h"]
C=h["C"]
shell=s["shell"]
inside=s["inside"]
block=s["block"]
pi=arb.pi()


def matrix_positive_parameter(u):
    # Positive Taylor coefficients of cosh and sinh(v)/v give monotone endpoint bounds.
    matrix=acb_mat(19,19)
    admittance={}
    for i,v in enumerate(shell["V"]):
        matrix[i,i]=u+v
    for i,j,d,z in shell["edges"]:
        if d not in admittance:
            argument=u*shell["rho"]*shell["theta"](d)/4
            lo,hi=argument.lower().max(0),argument.upper()
            denom=lo.hypgeom_0f1(arb(3)/2).union(hi.hypgeom_0f1(arb(3)/2))
            cosh=lo.hypgeom_0f1(arb(1)/2).union(hi.hypgeom_0f1(arb(1)/2))
            admittance[d]=1/denom,cosh/denom
        off,diag=admittance[d]
        matrix[i,i]+=diag
        matrix[j,j]+=diag
        matrix[i,j]=-off*z
        matrix[j,i]=matrix[i,j].conjugate()
    return matrix


def base_ok(low,high):
    m=matrix_positive_parameter(arb(low).union(high))
    interior=block(m,inside,inside)
    try:
        # Clock factors never enter a row connecting the base ring to the interior.
        y=-m[11,12]+(block(m,(11,),inside)*interior.solve(block(m,inside,(12,))))[0,0]
        if y.real>0 and arb(".034")<y.arg()<arb(".038"):
            return 1
    except ZeroDivisionError:
        pass
    assert high-low>arb(".0001")
    mid=(arb(low)+high)/2
    return base_ok(low,mid)+base_ok(mid,high)


if __name__=="__main__":
    assert not any(d in (724,804) and (i in (11,12) or j in (11,12))
                   for i,j,d,_ in shell["edges"])
    cells=base_ok(0,6400)
    direct=next(z.arg() for i,j,d,z in shell["edges"] if (i,j)==(11,12))
    edge_error=(arb(12)/(6400+4-11)).asin()
    assert arb(".034")<direct-edge_error and direct+edge_error<arb(".038")
    print("base phase in (.034,.038), compact interval cells",cells)

    faces=s["geometry"]["big"]
    geo=h["geo"]
    orig=[geo["geometry"]["lifted"](f,geo["negative"]) for f in faces]
    base_angle=arb(".034").union(arb(".038"))
    x=(pi-arb(".038")).union(pi)
    areas=[geo["area"](f) for f in faces]
    density_diff=[]
    for f,area,old in zip(faces,areas,orig):
        angle=2*pi-x-base_angle if len(f)==4 else x if 14 in f else 3*base_angle
        density_diff.append((angle-old)/area)
    magnetic_span=max(v.upper() for v in density_diff)-min(v.lower() for v in density_diff)
    diameter=max(abs(i-j).upper() for i in C for j in C)
    oscillation=magnetic_span*diameter**2/4
    # Fixed outward-rounded reference boxes, verified when the reference integral is built.
    reference=calc["reference_boxes"]()
    inflation=arb(0,(oscillation/2).upper()).exp()
    moment_bound=calc["coherence_upper"]([box*inflation for box in reference])
    assert moment_bound<arb(".165")
    # Use a rounded bound so the reported CKM limits reproduce independently of hardware.
    mu=arb(".165")/(1-arb(".165"))
    def direction(r):
        return (2*(1+mu)).sqrt()*mu*r/(1-2*mu-(1+mu)*r)
    up=arb(".00316")/arb(".8573")
    down=arb(".000271")/arb(".01392")
    error=direction(up)+direction(down)
    print("log-metric oscillation bound, moment contrast bound",oscillation.str(12),moment_bound.str(12))
    print("same-site Vcb and different-site Vtb bounds",error.str(12),(mu+arb(2).sqrt()*error).str(12))
    assert error<arb(".013") and mu+arb(2).sqrt()*error<arb(".217")
