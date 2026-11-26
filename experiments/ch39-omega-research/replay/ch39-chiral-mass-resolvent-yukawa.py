"""以已有手征质量映射与 Higgs 背景构造有限逆矩阵顶角并核验重阈值。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


filtered = run_path(str(Path(__file__).with_name("ch39-rational-yukawa-shell-filter.py")))
model = filtered["model"]
adjoint = model["adjoint"]


def norm_squared(matrix):
    return sum(abs(x)**2 for x in matrix.entries())


def quartic(matrix):
    return norm_squared(adjoint(matrix)*matrix)


def graph_response_check():
    # 非对角的有理质量映射，在原基底独立核验逆矩阵与核投影公式。
    mass = acb_mat([[1,2,0,1],[0,1,1,2]])
    squ = mass*adjoint(mass)
    eye = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    two = acb_mat([[1,0],[0,1]])
    p = eye-adjoint(mass)*squ.inv()*mass
    def zero(matrix):
        return all(x.contains(0) for x in matrix.entries())
    assert zero(p*p-p) and zero(mass*p)
    s = arb(5)/7
    f = eye-adjoint(mass)*(s*two+squ).inv()*mass
    assert zero((s*eye+adjoint(mass)*mass)*f-s*eye)
    assert zero(f*p-p)


def verify_case(number, coefficients, vacuum):
    parts, matrices = filtered["source"](number,model["column"](coefficients),model["column"](vacuum))
    v = arb("246.22")
    s = v*v
    weights = {q:[s/(s+mass**2) for mass in vals[1]] for q,vals in parts.items()}
    light, norms, quadratic = filtered["observables"](parts,matrices,weights)
    sextic = arb(0)
    quartic_remainder = arb(0)
    radial_norms = []
    for sector,(left,right) in enumerate(filtered["sectors"]):
        base = matrices[sector]
        colors = (3,3,1)[sector]
        ml,mr = parts[left][1],parts[right][1]
        fl,fr = weights[left],weights[right]
        eff = acb_mat([[fl[i]*fr[j]*base[i,j] for j in range(19)] for i in range(19)])
        assert all((eff[i,j]-base[i,j]).contains(0) for i in range(3) for j in range(3))
        assert all((light[sector][i,j]-base[i,j]).contains(0) for i in range(3) for j in range(3))
        excess = quartic(eff)-quartic(light[sector])
        assert excess > 0
        quartic_remainder += colors*excess
        radial_norms.append(sum(abs(eff[i,j]*(5-2*fl[i]-2*fr[j]))**2
                                for i in range(19) for j in range(19) if i>=3 or j>=3).sqrt())
        sextic += colors*sum(abs(base[i,j])**2/(ml[i]**2+mr[j]**2)
                             for i in range(19) for j in range(19) if (i<3)!=(j<3))
    excess_trace = v**2*quadratic + v**4*quartic_remainder/4
    assert arb("6.19e-10") < excess_trace < arb("6.37e-10")
    assert abs(excess_trace/(v**6*sextic)-1) < arb("1e-16")
    assert max(n.upper() for n in radial_norms) < arb("1.84e-21")
    # 完整 35 维 Dirac 矩阵相对于轻重块对角模型的谱扰动。
    delta = v/arb(2).sqrt()*max(n.upper() for n in norms)
    assert delta < arb("1.07e-19")
    for (left,right), y in zip(filtered["sectors"],light):
        sig = sorted(((v/arb(2).sqrt())*a.real.sqrt()
                      for a in (y*adjoint(y)).eig(algorithm="rump")), key=float)
        heavy_min = min(x.lower() for q in (left,right) for x in parts[q][1][3:])
        gap = min(x.lower() for x in (2*sig[0],sig[1]-sig[0],sig[2]-sig[1],heavy_min-sig[2]))
        assert gap > arb(".00097") and delta/sig[0] < arb("2.2e-16")
        assert 4*delta/(gap-delta) < arb("4.5e-16")
    print(number, "六次背景系数 GeV^-2、物理真空重层迹 GeV^4、含重模质量与径向顶角范数",
          *(v.str(10) for v in (sextic,excess_trace,max(n.upper() for n in norms),
                                 max(n.upper() for n in radial_norms))))
    print("相对六次主项", (excess_trace/(v**6*sextic)).str(12))


def verify():
    ctx.prec = 448
    graph_response_check()
    products = filtered["products"]
    for number,w,z in ((22027,model["WITNESS"],model["HIGGS"]),
                       (28325,products["SECOND_WITNESS"],products["SECOND_HIGGS"])):
        verify_case(number,w,z)


if __name__ == "__main__":
    verify()
