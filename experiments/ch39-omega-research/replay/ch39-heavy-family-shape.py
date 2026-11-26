"""Compare two subthreshold Euclidean polarizations for arbitrary positive heavy family spectra."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

finite = run_path(str(Path(__file__).with_name("ch39-subthreshold-triplet-polarization.py")))
mat, poles, geometry, q = (finite[k] for k in ("matrix","poles","geometry","q"))


def probe(multiplier, e, z, b, d, triplet):
    x = multiplier.exp()
    p = multiplier/(16*arb.pi()**2)
    correction = 48*(finite["scalar_integral"](x*e/triplet**2) -
                     finite["scalar_integral"](e/triplet**2))/(16*arb.pi()**2)
    A = mat["charge_trace"]()*p
    A[0,0] += correction
    return mat["massive_exchange"](x*e,z,b,d), A, x


def inverse_response(k, h, data):
    P, A, _ = data
    # F = P_XX + k v' (A^{-1} - kP - h yy')^{-1} v .
    v = arb_mat([[P[i,2]] for i in range(3)])
    kernel = A.inv()-P*k-mat["hypercharge"]*h
    assert mat["positive"](kernel)
    W = kernel.inv()
    vertex = W*v
    f = P[2,2]+k*(v.transpose()*vertex)[0,0]
    derivative = -(v.transpose()*vertex+k*vertex.transpose()*P*vertex)[0,0]/f**2
    return 1/f, derivative


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"]
    goals=[target(n,time) for n in (11430,11431,11432)]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        c1,c2,_,n1,_ = poles["reference_readings"](logs)
        e,b,d = n1**2/scale,c/rho,eta/rho
        probes=[probe(t,e,z,b,d,masses[0]) for t in (q("1/2"),q(1))]
        h=mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
        k=1/rho
        base = mat["massive_exchange"](e,z,b,d)[2,2]
        amplitudes = []
        for P,A,x in probes:
            T = P*k+mat["hypercharge"]*h
            amplitudes.append((T.inv()-A).inv()[2,2]/(k*x*base))
        n2 = geometry["second_neutral_mass"](scale,z,b,d)
        mass_error = arb(0,q("1e-10").upper())
        data=[m*(1+mass_error) for m in (c1,c2,n1,n2)]
        l,zz,_=poles["charged_geometry"](*data[:2])
        bb,dd=geometry["neutral_shape"](l,zz,data[2:])
        ee=data[2]**2/l
        triplet,_=geometry["physical_masses"](bb,arb(1),zz)
        readings=[probe(t,ee,zz,bb,dd,triplet) for t in (q("1/2"),q(1))]
        for P,_,_ in readings:
            assert mat["positive"](P)
        one,two=[row[2] for row in readings]
        threshold=4*triplet**2
        cap=(one-1)/(two-1)*(threshold+two*ee)/(threshold+one*ee)
        assert q(".4") < cap < q(".41")
        cap=cap.upper()
        origin=mat["massive_exchange"](ee,zz,bb,dd)[2,2]
        amp_error=arb(0,q("2e-8").upper())
        refs=[1/(amp*(1+amp_error)*x*origin) for amp,(_,_,x) in zip(amplitudes,readings,strict=True)]
        k30,k32=[1/(poles["offset"]+g/(15*l)) for g in (goals[0],goals[2])]
        # The normalized exchange increases with k and h on the positive branch.
        negative=inverse_response(k32,q(0),readings[1])[0]-refs[1]
        assert negative < 0
        k_range=arb(0).union(k30)
        margins,direction=[],[]
        for j in range(8):
            h_lo,h_hi=q(j)/8,q(j+1)/8
            h_range=h_lo.union(h_hi)
            derivative=[inverse_response(k_range,h_range,r)[1] for r in readings]
            direction.append((derivative[0]-cap*derivative[1]).upper())
            g1=inverse_response(k30,h_hi,readings[0])[0]-refs[0]
            g2=inverse_response(k30,h_lo,readings[1])[0]-refs[1]
            margins.append((g1-cap*g2).lower())
        print(case,"spectral ratio cap",cap.str(11),"direction upper",max(direction).str(9),
              "lower-label contrast",min(margins).str(9),"higher-label excess",negative.upper().str(9),flush=True)
        assert max(direction) < q("-.034") and min(margins) > q("7.2e-7")
        assert negative < q("-.00016")


if __name__ == "__main__":
    verify()
