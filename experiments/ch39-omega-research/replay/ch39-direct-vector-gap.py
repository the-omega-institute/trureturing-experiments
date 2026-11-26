"""限制直接交换的正谱质量，覆盖全部正零点，并夹住下邻编号可达的谱隙。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

vectors = run_path(str(Path(__file__).with_name("ch39-three-direct-vectors.py")))
poles, mat, contact, scheme, q = (vectors[k] for k in ("poles","mat","contact","scheme","q"))
projection = arb_mat([[3,-1,0],[0,0,1]])


def mixed(a, b):
    # 二维行列式的极化。
    return a[0,0]*b[1,1]+a[1,1]*b[0,0]-2*a[0,1]*b[0,1]


def det_polynomial(left, right, matrices, rays, coeffs):
    base = matrices[0]*left+matrices[1]*right
    (f1,s1),(f2,s2) = coeffs
    r1,r2 = rays[0]*(left*f1),rays[1]*(right*f2)
    b1,b2 = mixed(base,r1),mixed(base,r2)
    c = base.det()
    return c, c*(s1+s2)+b1+b2, c*s1*s2+b1*s2+b2*s1+mixed(r1,r2)


def verify_case(logs, time, target):
    (length,z,rho,raw_d,raw_b,masses), _, _ = poles["retune"]["predictions"](logs,time)
    e, x = scheme["anchor"]/length, q("5/3")
    k,b,d = 1/rho,raw_b/rho,raw_d/rho
    t = rho/(poles["offset"]+target/(15*length))
    assert 0 < t < 1
    masses_sq = masses[0]**2
    mats, rays, coeffs, p, all_mats, all_rays = [], [], [], [], [], []
    for en in (e,x*e):
        side, _, v, _, _ = contact["axes"](en,z)
        ell = b+d*en+4/side
        tree = mat["massive_exchange"](en,z,b,d)
        p.append(en/ell)
        all_mats.append(tree)
        all_rays.append(v*v.transpose())
        mats.append(projection*tree*projection.transpose())
        rays.append(projection*v*v.transpose()*projection.transpose())
        coeffs.append((t*d*(en/ell)**2,(1-t)*d*en/ell))
    rate = p[1]/p[0]
    assert rate > 1

    def min_ratio(gap):
        squared = gap**2*masses_sq
        return x*(e+squared)/(x*e+squared)

    # 全部较低编号：D'<=D 时族流比至多为 rate。
    # D'>D 时放大第二个秩一项至切线给出的最大值。
    lower = min_ratio(q("2"))
    assert rate < lower
    base = mats[1]-mats[0]*lower
    slope = rays[1]*rate**2-rays[0]*lower
    quadratic, linear, constant = slope.det(), mixed(base,slope), base.det()
    assert quadratic < 0 and linear**2-4*quadratic*constant < 0
    print("all lower labels, family ratio and required ratio", rate.str(10), lower.str(10),flush=True)

    # 对 11430 保留两个秩一项间的精确有理关系。
    sep = q("1.2")
    lower = min_ratio(q(".644"))
    for left,right,side in ((-lower,q(1),"low"),(x,q(-1),"high")):
        c1,c2,c3 = det_polynomial(left,right,mats,rays,coeffs)
        value = c1+c2*sep+c3*sep**2
        derivative = c2+2*c3*sep
        assert c3 < 0 and value < 0
        if side == "low":
            assert derivative > 0
        else:
            assert derivative < 0
    # 选一个具体的重配，构造位于质量下界之上的正谱。
    delta = q("1")
    shifts = [tree+ray*(f*delta/(1+s*delta)) for tree,ray,(f,s) in
              zip(all_mats,all_rays,coeffs,strict=True)]
    assert mat["positive"](x*shifts[0]-shifts[1])
    assert mat["positive"](shifts[1]-min_ratio(q(".640"))*shifts[0])
    ratios = [mass/masses[0] for mass in vectors["effective_masses"](*shifts,e,x)]
    print("retuned relative masses", *(r.str(11) for r in ratios), flush=True)
    assert min(ratios) > q(".640")


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11430,time)
    for logs in poles["retune"]["log_sets"]:
        verify_case(logs,time,target)


if __name__ == "__main__":
    verify()
