"""测量最轻中性质量和颜色三重态，排除所有较高零点的混合局部矢量谱。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat, arb_poly

family = run_path(str(Path(__file__).with_name("ch39-retuned-family-spectral-obstruction.py")))
poles, mat, q = (family[key] for key in ("poles", "mat", "q"))
edges = (
    (("9.360590","9.061786"), ("5.184614","5.744860")),
    (("9.368670","9.069269"), ("5.189996","5.750391")),
    (("9.403785","9.101795"), ("5.213382","5.774430")),
)


def verify():
    time = poles["study"]["nominal"][2]
    goal = poles["study"]["comp"]["model"]["target"](11432,time)
    *_, colour, right = poles["neutral"]["graph_constants"](port_builder=poles["neutral"]["pole_polynomials"])
    (pc,vc),(pr,vr) = ([arb_poly(p) for p in pair] for pair in (colour,right))
    for i, logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        b0,d0,m = c/rho,eta/rho,masses[0]
        neutral = poles["reference_readings"](logs)[3]
        e = neutral**2/scale
        h = mat["flow"]["running"]["hypercharge_window"](time+neutral.log()).mid()**2
        k = 1/(poles["offset"]+goal/(15*scale))
        refs = family["model"]["retuned"]["readings"](e,z,b0,d0,m,1/rho,h,q(0))
        origin = family["bare_data"](e,z,q(0),q(0),k,h)[0]
        for mass_err, pair in zip((q("1e-10"),q("1e-6")),edges[i],strict=True):
            cut = e*(1-mass_err)**2
            col, weak = family["model"]["ports_at"](cut)
            y = -4/(z+9*col+weak)
            assert cut < 4*m*m and y > 0
            for log_x, ref, edge in zip((q(".5"),q(1)),refs,pair,strict=True):
                x = log_x.exp()
                s = x*e
                r = ref*(1+q("2e-8"))
                fixed, base, leg = family["bare_data"](s,z,q(0),q(0),k,h)
                _, slope, _ = family["bare_data"](s,z,q(0),q(0),q(1),q(0))
                two = arb_mat([[slope[a,b] for b in range(2)] for a in range(2)])
                assert two[0,0] > 0 and two.det() > 0
                assert e-r*s < 0 and 1-r+(e-r*s)/cut < 0

                def slack(change):
                    b = b0*(1+change)**2
                    d = (b-y)/cut
                    ell0 = b+d*e+origin
                    ell = b+d*s+fixed
                    pol = family["correction"](e,x,m**2*(b/b0))
                    dressed = family["dressed_term"](k,s,base,leg,log_x,pol)
                    assert d > 0 and ell > dressed and ell0 > 0
                    return ell0-r*(ell-dressed)

                center = q(edge)/10**6
                radius = q("1e-12")
                assert slack(center-radius)>0>slack(center+radius)
                # The zero-hidden-field construction has exactly one massless
                # pole and its first positive tree pole is the prescribed cut.
                box = (center-radius).union(center+radius)
                b = b0*(1+box)**2
                d = (b-y)/cut
                old = ((vc*z-pc*(9/(24*mat["unit"])))*vr-pr*vc*(1/(4*mat["unit"])))*arb_poly([b,-d*mat["unit"]])+4*vc*vr
                assert old[0]==0 and not old[1].contains(0) and old(cut/mat["unit"]).contains(0)
                positive = arb_poly(old.coeffs()[1:])
                variations = family["model"]["ordered"]["variations_in"]
                assert variations(positive,q(0),cut*q("1.001")/mat["unit"]) == 1
            print("case",i,"neutral mass relative error",mass_err.str(3),
                  "individual triplet thresholds (ppm)",*pair,flush=True)


if __name__ == "__main__":
    verify()
