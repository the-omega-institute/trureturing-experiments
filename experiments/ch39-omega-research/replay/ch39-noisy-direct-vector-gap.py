"""在有限精度的极化轻子率上排除任意多个质量至少为两倍三重态的直接矢量。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

root = Path(__file__).parent
noisy = run_path(str(root/"ch39-noisy-contact-direction.py"))
gap = run_path(str(root/"ch39-direct-vector-gap.py"))
poles, mat, contact, scheme, light, q = (noisy[key] for key in ("poles","mat","contacts","scheme","light","q"))
window = scheme["uncertain"]


def verify_case(logs, time, targets):
    (scale,z0,rho,raw_d,raw_b,masses), _, _ = poles["retune"]["predictions"](logs,time)
    c1,c2 = poles["reference_readings"](logs)[:2]
    mass_window = window(q("1e-11"))
    length,z,_ = poles["charged_geometry"](c1*mass_window, c2*mass_window)
    triplet = masses[0]*(scale/length).sqrt()*window(q("1e-7"))
    test_b = q("1/5")
    base_mass,_ = mat["geometry"]["physical_masses"](test_b,q(1),z)
    b = test_b*(triplet/base_mass)**2
    d, k0 = raw_d/rho, 1/rho
    e, e0 = scheme["anchor"]/length, scheme["anchor"]/scale
    low,high = [1/(poles["offset"]+t/(15*length)) for t in targets]
    diff = k0-low
    assert diff > 0 and high > k0
    graphs = light["graphs"](logs,time)
    log_x = scheme["logs_x"][0]
    _,polar,x = scheme["cal"]["shape"]["probe"](log_x,e,z,b,d,triplet)
    _,polar0,_ = scheme["cal"]["shape"]["probe"](log_x,e0,z0,raw_b/rho,d,masses[0])
    side1 = contact["axes"](e,z)[0]
    project = arb_mat([[3,-1,2*(1-z/side1)],[0,0,1]])
    masses_part, sides, p = [], [], []

    def projected_tree(energy, shape, family_b, family_d):
        side, _, _, _, beta = contact["axes"](energy,shape)
        spatial = side-shape
        coefficient = 2*(shape/side-z/side1)
        family = arb_mat([[coefficient], [1]])
        aligned = arb_mat([[beta*spatial**2,0],[0,0]])
        return aligned + family*family.transpose()*(energy/(family_b+family_d*energy+4/side))

    for i,mult in enumerate((q(1),x)):
        energy = mult*e
        side, _, _, _, beta = contact["axes"](energy,z)
        tree = projected_tree(energy,z,b,d)
        center = projected_tree(mult*e0,z0,raw_b/rho,d)
        graph = graphs[i]
        data_error = noisy["measurement_error"](graph, q("1e-11"))
        eye = scheme["cal"]["eye"]
        if i == 0:
            observed = graph+data_error
            error = data_error
            inverse = noisy["dual_basis"](energy,z)
            aligned = inverse*observed*inverse.transpose()
            assert aligned[1,1] < high*beta
        else:
            varied = graph+data_error
            error = (eye+varied*polar).inv()*(data_error-varied*(polar-polar0)*graph)*(eye+polar0*graph).inv()
        delta = project*error*project.transpose()+(center-tree)*k0
        uncertain = arb_mat([[arb(0,(abs(delta[r,j])/diff).upper()) for j in range(2)] for r in range(2)])
        masses_part.append(tree+uncertain)
        sides.append(side)
        p.append(energy/(b+d*energy+4/side))

    a = x*(4*triplet**2+e)/(4*triplet**2+x*e)
    intercept = masses_part[1]-masses_part[0]*a
    wedge = 2*z*(sides[0]-sides[1])/(sides[0]*sides[1])
    ratio_sq = (p[1]/p[0])**2
    slope = arb_mat([[ratio_sq*wedge**2,ratio_sq*wedge],[ratio_sq*wedge,ratio_sq-a]])
    quadratic,linear,constant = -a*ratio_sq*wedge**2,gap["mixed"](intercept,slope),intercept.det()
    discriminant = linear**2-4*quadratic*constant
    print("negative determinant coefficients", quadratic.str(10), discriminant.str(10), flush=True)
    assert quadratic < 0 and discriminant < 0


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11432)]
    assert poles["verify_ratio_monotonicity"]() > 0
    for logs in poles["retune"]["log_sets"]:
        verify_case(logs,time,targets)


if __name__ == "__main__":
    verify()
