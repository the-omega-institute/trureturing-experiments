"""用两处保极点插值约束任意多个只与族流混合的局部矢量。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

model = run_path(str(Path(__file__).with_name("ch39-retuned-dark-current-match.py")))
mat, poles, q = (model[key] for key in ("mat", "poles", "q"))
scalar_integral = model["retuned"]["shape"]["finite"]["scalar_integral"]


def bare_data(s, z, b0, d0, k, h):
    c, r = model["ports_at"](-s)
    side = z+9*c+r
    direction = arb_mat([[3*c],[-r],[0]])
    bare = (arb_mat([[c,0,0],[0,r,0],[0,0,0]])-model["tensor"](direction)*(1/side))*k*s
    bare += mat["hypercharge"]*(h-k/24)
    leg = arb_mat([[-6*c/side],[2*r/side],[1]])
    return b0+d0*s+4/side, bare, leg


def correction(e, x, mass_squared):
    return 48*(scalar_integral(x*e/mass_squared)-scalar_integral(e/mass_squared))/(16*arb.pi()**2)


def dressed_term(k, s, bare, leg, log_x, extra):
    trace = mat["charge_trace"]()*log_x*(1/(16*arb.pi()**2))
    trace[0,0] += extra
    inverse = trace.inv()-bare
    assert mat["positive"](inverse)
    return k*s*model["scalar"](leg,inverse.inv(),leg)


def check_case(logs, label, target):
    time = poles["study"]["nominal"][2]
    (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
    b0, d0, m = c/rho,eta/rho,masses[0]
    n1 = poles["reference_readings"](logs)[3]
    e = n1**2/scale
    second = mat["geometry"]["second_neutral_mass"](scale,z,b0,d0)**2/scale
    h = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
    k = 1/(poles["offset"]+target/(15*scale))
    refs = model["retuned"]["readings"](e,z,b0,d0,m,1/rho,h,q(0))
    ell_e = bare_data(e,z,b0,d0,k,h)[0]
    e_tol = q("2e-8")
    readings = [refs[0]*(1+arb(0,e_tol.upper())), refs[1]*(1+q("1e-7"))]
    data = []
    for i, log_x in enumerate((q(".5"), q(1))):
        x = log_x.exp()
        s = x*e
        ell, bare, leg = bare_data(s,z,b0,d0,k,h)
        # Positive branch uniformly up to the massless triplet limit.
        full = 16*log_x/(16*arb.pi()**2)
        dressed_term(k,s,bare,leg,log_x,full)
        initial = correction(e,x,m**2)
        data.append((ell,bare,leg,log_x,x,s,initial))
    def kernel(s, w, ratio):
        def component(v):
            return (v+e)*(v+second)/(e*second*(1+w*v/second))
        return component(s)-component(e)/ratio
    # r1 in its permitted interval. For all shrinkages the known triplet
    # polarization grows, so T1 cannot exceed its value at the reference mass.
    ell,bare,leg,log_x,x,s,initial = data[0]
    t1_max = ell-dressed_term(k,s,bare,leg,log_x,initial)-ell_e/readings[0]
    upper = correction(e,x,m**2*(1-q(".00005")))
    t1_min = ell-dressed_term(k,s,bare,leg,log_x,upper)-ell_e/readings[0]
    ratio = q("3.455")
    for i in range(2048):
        w = (q(i)/2048).union(q(i+1)/2048)
        g1 = kernel(s,w,readings[0])
        assert g1 > 0 and t1_max < b0*q(".00005")*g1
        assert b0*q(".000024")*g1 < t1_min
        g2 = kernel(data[1][5],w,readings[1])
        assert g2 > ratio*g1
    # b/b0 lies in [1-0.00005,1]. Bound each resolvent at both ends,
    # using monotonicity of the once-subtracted triplet correction.
    limits = []
    for (ell,bare,leg,log_x,x,s,initial),r in zip(data,readings,strict=True):
        upper = correction(e,x,m**2*(1-q(".00005")))
        correction_box = initial.union(upper)
        limits.append(ell-dressed_term(k,s,bare,leg,log_x,correction_box)-ell_e/r)
    assert limits[1] < ratio*limits[0]
    print("label",label,"first family reading within 20 ppb forces second above 100 ppb excess",flush=True)


if __name__ == "__main__":
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11430,time)
    for logs in poles["retune"]["log_sets"]:
        check_case(logs,11430,target)
