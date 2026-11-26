"""Sum light-field logarithmic polarization with all neutral graph exchanges."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

flow = run_path(str(Path(__file__).with_name("ch39-relative-family-running.py")))
geometry, poles, q = (flow[key] for key in ("geometry", "poles", "q"))
unit = poles["unit"]
y = arb_mat([[q("1/6")], [q("1/2")], [0]])
vev = arb_mat([[3], [-1], [2]])
hypercharge = y*y.transpose()


def check_massless_residue():
    _,_,(pc,vc),(pr,vr) = poles["neutral"]["graph_constants"](
        port_builder=poles["neutral"]["pole_polynomials"])
    assert 15*pc(0) == vc.derivative()(0) and 14*pr(0) == vr.derivative()(0)
    residue = arb_mat([[q("1/360"),0,0],[0,q("1/56"),0],[0,0,0]])
    broken = residue*vev
    massless = residue-broken*broken.transpose()*(1/(vev.transpose()*broken)[0,0])
    for value in (massless-hypercharge*q("1/24")).entries():
        assert value.contains(0) and abs(value) < q("1e-45")


def charge_trace():
    matrix = arb_mat(3,3)
    fields = poles["study"]["retuned"]["current"]["weyl_charges"]
    for multiplicity,(d,r,sign) in zip((6,3,3,2,1),fields.values(),strict=True):
        for x in (1,1,-2):
            charge = arb_mat([[d],[r],[sign*x]])
            matrix += (charge*charge.transpose()) * (q("2/3")*multiplicity)
    higgs = arb_mat([[0],[1],[0]])
    matrix += (higgs*higgs.transpose())*q("2/3")
    assert ((y.transpose()*matrix*y)[0,0]-q("41/6")).contains(0)
    assert positive(matrix)
    return matrix


def positive(matrix):
    return (matrix[0,0] > 0 and arb_mat([[matrix[i,j] for j in range(2)] for i in range(2)]).det()>0
            and matrix.det()>0)


def massive_exchange(s, z, b, d):
    bare = arb_mat(3,3)
    bare[0,0] = poles["colour"](-s,unit)/24
    bare[1,1] = poles["right"](-s,unit)/4
    bare[2,2] = 1/(b+d*s)
    coupling = bare*vev
    graph = bare - coupling*coupling.transpose()*(1/(z+(vev.transpose()*coupling)[0,0]))
    result = graph*s - hypercharge*(q("1/24"))
    return result


def enhancement(P, inverse_trace, k, h, p):
    column = arb_mat([[P[i,2]] for i in range(3)])
    kinetic = inverse_trace-(P*k+hypercharge*h)*p
    assert positive(kinetic)
    return P[2,2]+p*k*(column.transpose()*kinetic.inv()*column)[0,0]


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11431,11432)]
    check_massless_residue()
    trace = charge_trace()
    inverse_trace = trace.inv()
    x = q(1).exp()
    p = 2*flow["step"]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,_),_,_ = poles["retune"]["predictions"](logs,time)
        c1,c2,_,n1,_ = poles["reference_readings"](logs)
        e = n1**2/scale
        b,d = c/rho,eta/rho
        n2 = geometry["second_neutral_mass"](scale,z,b,d)
        G0,G1 = (massive_exchange(s,z,b,d) for s in (e,x*e))
        start = G0[2,2]
        full_tree = 1/(b+d*e+4/(z+geometry["green"](-e)))
        assert abs(start/e-full_tree) < q("1e-15")
        assert positive(G0) and positive(G1)
        gy = flow["running"]["hypercharge_window"](time+n1.log()).mid()
        T = G1*(1/rho)+hypercharge*gy**2
        dressed = (T.inv()-trace*p).inv()
        ref = rho*dressed[2,2]/(x*start)
        check = enhancement(G1,inverse_trace,1/rho,gy**2,p)/(x*start)
        assert abs(check-ref)<q("1e-12") and check.overlaps(ref)
        print(case,"full normalized amplitude",ref.str(13),"tree", (G1[2,2]/(x*start)).str(13))
        err = arb(0,q("1e-10").upper())
        data = [m*(1+err) for m in (c1,c2,n1,n2)]
        fitted_scale,zz,_ = poles["charged_geometry"](*data[:2])
        bb,dd = geometry["neutral_shape"](fitted_scale,zz,data[2:])
        ee = data[2]**2/fitted_scale
        triplet,_ = geometry["physical_masses"](bb,arb(1),zz)
        assert (triplet/ee.sqrt()).log() > q("1/2")
        base = massive_exchange(ee,zz,bb,dd)[2,2]*x
        matrix = massive_exchange(ee*x,zz,bb,dd)
        assert positive(matrix) and base > 0
        ks = [1/(poles["offset"]+n/(15*fitted_scale)) for n in targets]
        low = enhancement(matrix,inverse_trace,ks[0],q(1),p)/base
        high = enhancement(matrix,inverse_trace,ks[2],q(0),p)/base
        print("outward lower margins (ppm)", ((1-low/ref)*10**6).lower().str(10),
              ((high/ref-1)*10**6).lower().str(10),flush=True)
        assert low < ref*(1-q("2.05e-6")) and high > ref*(1+q("1.40e-6"))


if __name__ == "__main__":
    verify()
