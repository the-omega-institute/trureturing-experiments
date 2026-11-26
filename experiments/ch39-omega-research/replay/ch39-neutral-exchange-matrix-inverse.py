"""Separate zero labels even with an arbitrary gapped positive matrix spectral measure."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

mixed = run_path(str(Path(__file__).with_name("ch39-mixed-spectral-witness.py")))
shape, poles, mat, q = (mixed[key] for key in ("shape","poles","mat","q"))


def coefficients(energy, b, d, z, probe, observed, start):
    _,A,x = probe
    s = x*energy
    # The rank-one hypercharge change drops out of the XX entry of the
    # inverse graph exchange. The breaking charge is (3,-1,2).
    effective = observed.inv()[2,2]/(x*start)
    return d+(b+4/z)/s-effective, A[2,2]


def verify():
    time=poles["study"]["nominal"][2]
    goals=[poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11432)]
    times = (q("1/2"),q(1))
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        c1,c2,_,n1,_ = poles["reference_readings"](logs)
        e,b,d = n1**2/scale,c/rho,eta/rho
        n2 = shape["geometry"]["second_neutral_mass"](scale,z,b,d)
        probes=[shape["probe"](t,e,z,b,d,masses[0]) for t in times]
        k=1/rho
        h=mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
        start=mat["massive_exchange"](e,z,b,d)[2,2]
        observed = []
        for P,A,x in probes:
            raw=P*k+mat["hypercharge"]*h
            observed.append((raw.inv()-A).inv()*(1/(x*k*start)))
        err=arb(0,q("1e-15").upper())
        masses_box=[m*(1+err) for m in (c1,c2,n1,n2)]
        l,zz,_=poles["charged_geometry"](*masses_box[:2])
        bb,dd=shape["geometry"]["neutral_shape"](l,zz,masses_box[2:])
        ee=masses_box[2]**2/l
        triplet,_=shape["geometry"]["physical_masses"](bb,arb(1),zz)
        readings=[shape["probe"](t,ee,zz,bb,dd,triplet) for t in times]
        origin=mat["massive_exchange"](ee,zz,bb,dd)[2,2]
        one,two=(r[2] for r in readings)
        threshold=4*triplet**2
        cap=((one-1)/(two-1)*(threshold+two*ee)/(threshold+one*ee)).upper()
        u=arb(0,q("1e-14").upper())
        box=[arb_mat([[R[i,j]*(1+u) for j in range(3)] for i in range(3)]) for R in observed]
        (a1,b1),(a2,b2)=(coefficients(ee,bb,dd,zz,r,obs,origin) for r,obs in zip(readings,box,strict=True))
        k30,k32=[1/(poles["offset"]+g/(15*l)) for g in goals]
        slope=b1-cap*b2
        margin=a1-cap*a2-k30*slope
        high=a2-k32*b2
        print(case,"slope",slope.lower().str(8),
              "low-index margin",margin.lower().str(9),"high-index upper bound",high.upper().str(9),flush=True)
        assert cap < q(".406427") and slope > q(".035")
        assert margin > q(".000025") and high < q("-.00016")
        for (P,A,x),measurement in zip(probes,observed,strict=True):
            algebra = d+(b+4/z)/(x*e)-k*A[2,2]
            direct = measurement.inv()[2,2]/(x*start)
            assert direct.overlaps(algebra) and abs(direct-algebra)<q("1e-6")



if __name__ == "__main__":
    verify()
