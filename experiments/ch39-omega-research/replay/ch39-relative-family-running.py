"""Recover the site normalization from a traceless current at two one-loop scales."""
from pathlib import Path
from runpy import run_path

from flint import arb

geometry = run_path(str(Path(__file__).with_name("ch39-isospectral-resonance-masses.py")))
poles, q, unit = (geometry[key] for key in ("poles", "q", "unit"))
running = poles["window"]["running"]
step = 1/(32*arb.pi()**2)


def normalized_charges(energy, b, d):
    delta = b-d*energy
    family = -1/(d+geometry["green"](energy,True)*delta**2/4).sqrt()
    alpha = family*poles["colour"](energy,unit)*delta/16
    beta = -family*poles["right"](energy,unit)*delta/8
    return alpha, beta, family


def recover_polarization(ratios, charges):
    yy, yv, vv = running["polarization"](*charges)
    k = yv**2/(yy*vv)
    assert 0 < k < 1
    s1, s2 = ((1-1/r)/(i*step) for i,r in enumerate(ratios, 1))
    change = s2/s1-1
    assert s1 > 0 and change > 0
    def equation(u):
        return k*u-change*(1-2*u)*(1-(1-k)*u)
    u = poles["roots"](equation,0,".5")
    amplitude = s1-(1-2*u)*(s2-s1)
    assert amplitude > 0
    return vv/amplitude, (u/(step*yy)).sqrt()


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"]
    targets = [target(n,time) for n in (11430,11431,11432)]
    poles["verify_ratio_monotonicity"]()
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        model,_,_ = poles["retune"]["predictions"](logs,time)
        scale,z,rho,eta,c,masses = model
        c1,c2,_,n1,_ = poles["reference_readings"](logs)
        n2 = geometry["second_neutral_mass"](scale,z,c/rho,eta/rho)
        onset = time+n1.log()
        gy = running["hypercharge_window"](onset).mid()
        physical = poles["study"]["retuned"]["current"]["light_current"](logs)[2]
        predicted = normalized_charges(n1**2/scale,c/rho,eta/rho)
        assert all((v/rho.sqrt()).overlaps(w) for v,w in zip(predicted,physical,strict=True))
        exact = [(running["run_window"](t,gy,physical)[4][2]/physical[2])**2
                 for t in (q("1/4"),q("1/2"))]
        print("group",case,"gY",gy.str(8),"ratios",*(r.str(13) for r in exact))
        mass_error = arb(0,q("1e-8").upper())
        obs = [m*(1+mass_error) for m in (c1,c2,n1,n2)]
        fitted_scale,fitted_z,_ = poles["charged_geometry"](*obs[:2])
        b,d = geometry["neutral_shape"](fitted_scale,fitted_z,obs[2:])
        triplet,_ = geometry["physical_masses"](b,arb(1),fitted_z)
        assert (triplet*fitted_scale.sqrt()/obs[2]).log() > q("1/2")
        charges = normalized_charges(obs[2]**2/fitted_scale,b,d)
        # The two-scale inverse is exact for any positive hypercharge coupling
        # on the positive kinetic branch. Curvature can be smaller than the
        # input error; the inequality below covers that regime as well.
        fitted_rho,gauge_y = recover_polarization(exact, predicted)
        assert abs(fitted_rho-rho) < q("1e-10") and abs(gauge_y-gy) < q("1e-7")
        # Retain one massive neutral vector and the massless hypercharge field.
        # An external lepton flavour difference cancels its hypercharge vertex.
        # A common Euclidean subtraction mass fixed to M_N gives A/A0=2/(1+x/R).
        xs = [q("1/2").exp(), q(1).exp()]
        amplitudes = [2/(1+x/r) for x,r in zip(xs,exact,strict=True)]
        amplitude_error = arb(0,q("1.5e-7").upper())
        for x,amp,truth in zip(xs,amplitudes,exact,strict=True):
            recovered = x/(2/(amp*(1+amplitude_error))-1)
            assert abs(recovered/truth-1) < q("2.5e-7")
        print("normalized Euclidean exchange", *(v.str(12) for v in amplitudes))
        yy,yv,vv = running["polarization"](*charges)
        ratio_error = arb(0,q("2.5e-7").upper())
        data = [r*(1+ratio_error) for r in exact]
        s1,s2 = [(1-1/r)/(i*step) for i,r in enumerate(data,1)]
        amplitude = (2*s1-s2).lower().union(s1.upper())
        assert amplitude > 0
        rho_two = vv/amplitude
        # One upper-scale ratio is enough if the independent bound g_Y <= 1 is supplied.
        one_mix = (2*step*yv**2/(1-2*step*yy)).upper()
        rho_one = (vv+arb(0).union(one_mix))/s2
        for label, inferred in (("two",rho_two),("one with gY bound",rho_one)):
            interval = 15*fitted_scale*(inferred-poles["offset"])
            print(label,"target",interval.lower().str(12),interval.upper().str(12))
            assert targets[0] > interval > targets[2] and interval.contains(targets[1])
        # At the retained upper-scale datum there is also an exact neighbour
        # with the same masses if hypercharge is independently unconstrained.
        yy,yv,vv = running["polarization"](*predicted)
        alternate_rho = poles["offset"]+targets[0]/(15*scale)
        measured_s2 = (1-1/exact[1])/(2*step)
        excess = alternate_rho*measured_s2-vv
        assert excess > 0
        y_sq = excess/(2*step*(yv**2+yy*excess))
        assert 2 < y_sq.sqrt() < 3 and 2*step*yy*y_sq < 1
        alternate_charges = tuple(v/alternate_rho.sqrt() for v in predicted)
        family = running["run_window"](q("1/2"),y_sq.sqrt(),alternate_charges)[4][2]
        rebuilt = (family/alternate_charges[2])**2
        assert rebuilt.overlaps(exact[1]) and abs(rebuilt-exact[1]) < q("1e-12")
        print("single-scale 11430 witness g_Y",y_sq.sqrt().str(12))


if __name__ == "__main__":
    verify()
