"""用共同有理平方动量与减法点源比校准中性形状及超荷。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

cal = run_path(str(Path(__file__).with_name("ch39-calibrated-mixed-current.py")))
q, poles, mat = (cal[key] for key in ("q", "poles", "mat"))
uncertain = cal["uncertain"]
# All models use the same externally specified momentum unit.
anchor = q("1/450")
logs_x = (q("5/3").log(), q("8/3").log())
shape_source = cal["inv_basis"]*arb_mat([[3],[-1],[2]])
hyper_source = cal["calibration"]


def collect(logs, time, targets, mass_error, shape_error, source=cal["mixed"]):
    (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
    c1 = poles["reference_readings"](logs)[0]
    b,d,k0 = c/rho,eta/rho,1/rho
    e = anchor/scale
    h0 = mat["flow"]["running"]["hypercharge_window"](time+anchor.log()/2).mid()**2
    ref0 = cal["projected_tree"](e,z,b,d)+cal["hypercharge"]*(h0/k0)

    def ratio(source):
        return cal["form"](source,ref0,source)/ref0[0,0]

    # The octet-triplet splitting grows as z^(21/5), while the
    # first charged squared mass in graph units decreases with z.
    m3, m8 = [v*scale.sqrt()*uncertain(mass_error) for v in masses]
    splitting = m8**2-m3**2
    gamma = splitting/(c1*uncertain(mass_error))**2
    gamma0 = (masses[1]**2-masses[0]**2)*scale/c1**2
    zz = z*q(1).union((gamma/gamma0)**q("5/21"))
    gap_squared = ((2+21*(56*zz).log())/5).exp()
    length = splitting/gap_squared
    ee = anchor/length
    triplet = m3/length.sqrt()
    test_b = q("1/5")
    unit_triplet,_ = mat["geometry"]["physical_masses"](test_b,q(1),zz)
    bb = test_b*(triplet/unit_triplet)**2

    c_site = poles["colour"](-ee,mat["unit"])/24
    r_site = poles["right"](-ee,mat["unit"])/4
    sites = 9*c_site+r_site
    side = zz+sites
    observed = ratio(shape_source)*uncertain(shape_error)
    ell = (observed-4*(zz/side)**2)*side/(zz*sites)
    dd = (ell-bb-4/side)/ee
    assert dd > 0
    origin = cal["projected_tree"](ee,zz,bb,dd)
    t = (ratio(hyper_source)*uncertain(q("1e-6"))*origin[0,0]
         -cal["form"](hyper_source,origin,hyper_source))/cal["form"](hyper_source,cal["hypercharge"],hyper_source)
    assert t > 0 and dd.contains(d) and t.contains(h0/k0)
    origin += cal["hypercharge"]*t
    records = cal["probe_records"]((e,z,b,d,masses[0],h0,k0,ref0),
                                   (ee,zz,bb,dd,triplet,t,origin,length),targets,logs_x,source)
    return records, dd


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11432)]
    for i,logs in enumerate(poles["retune"]["log_sets"]):
        (records,k0,ks,cap,t),d = collect(logs,time,targets,q("1e-11"),q("1e-10"))
        count, high = cal["identify"](records,k0,ks,cap)
        print("case",i,"shape",d.str(10),"gY^2/k",t.str(9),"lower boxes",count,
              "upper excess ppm",(high*10**6).lower().str(9),flush=True)


if __name__ == "__main__":
    verify()
