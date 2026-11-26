"""Include every tree-level family-current pole in the two-scale Euclidean test."""
from pathlib import Path
from runpy import run_path

from flint import arb

flow = run_path(str(Path(__file__).with_name("ch39-relative-family-running.py")))
geometry, poles, q = (flow[key] for key in ("geometry", "poles", "q"))
green = geometry["green"]
xs = (q("1/2").exp(), q(1).exp())


def components(energy, z, b, d):
    weight = 1/(d + green(energy,True)*(b-d*energy)**2/4)

    def exchange(x):
        s = energy*x
        return 1/(b+d*s+4/(z+green(-s)))

    whole = exchange(arb(1))
    light = weight/(2*energy)
    baseline = whole/light
    # Normalize all tree amplitudes to the first pole's value at M_N.
    heavy = [exchange(x)/light - 2/(1+x) for x in xs]
    assert baseline > 1 and min(heavy) > 0
    return baseline, heavy


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11431,11432)]
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,_),_,_ = poles["retune"]["predictions"](logs,time)
        c1,c2,_,n1,_ = poles["reference_readings"](logs)
        e = n1**2/scale
        b,d = c/rho,eta/rho
        n2 = geometry["second_neutral_mass"](scale,z,b,d)
        total,heavy = components(e,z,b,d)
        gy = flow["running"]["hypercharge_window"](time+n1.log()).mid()
        currents = poles["study"]["retuned"]["current"]["light_current"](logs)[2]
        ratios = [(flow["running"]["run_window"](t,gy,currents)[4][2]/currents[2])**2
                  for t in (q("1/4"),q("1/2"))]
        amplitudes = [(h+2/(1+x/r))/total for x,h,r in zip(xs,heavy,ratios,strict=True)]
        diffs = [(h+2/(1+x))/(total*2/(1+x))-1 for x,h in zip(xs,heavy,strict=True)]
        naive = [x/(2/a-1) for x,a in zip(xs,amplitudes,strict=True)]
        print(case,"heavy fraction at M_N", (1-1/total).str(12),
              "relative tree shape change",*(v.str(12) for v in diffs))
        print("naive corrected-current ratios", *(v.str(12) for v in naive))
        assert q(".0927") < 1-1/total < q(".0934")
        assert q(".0289") < diffs[0] < q(".0292")
        assert q(".0756") < diffs[1] < q(".0761")
        error = arb(0,q("1e-10").upper())
        data = [m*(1+error) for m in (c1,c2,n1,n2)]
        fitted_scale,fitted_z,_ = poles["charged_geometry"](*data[:2])
        b,d = geometry["neutral_shape"](fitted_scale,fitted_z,data[2:])
        fitted_e = data[2]**2/fitted_scale
        baseline,background = components(fitted_e,fitted_z,b,d)
        amps = [a*(1+arb(0,q("1e-7").upper())) for a in amplitudes]
        coeff = [x/(2/(baseline*a-h)-1) for x,a,h in zip(xs,amps,background,strict=True)]
        s1,s2 = [(1-1/r)/(i*flow["step"]) for i,r in enumerate(coeff,1)]
        strengths = (2*s1-s2).lower().union(s1.upper())
        charges = flow["normalized_charges"](fitted_e,b,d)
        vv = flow["running"]["polarization"](*charges)[2]
        answer = 15*fitted_scale*(vv/strengths-poles["offset"])
        assert targets[0]>answer>targets[2] and answer.contains(targets[1])
        print("target",answer.lower().str(9),answer.upper().str(9),
              targets[0]>answer>targets[2],flush=True)


if __name__ == "__main__":
    verify()
