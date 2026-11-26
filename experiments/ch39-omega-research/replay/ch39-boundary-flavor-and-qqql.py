"""局部 Yukawa 变化幅度与三边界 QQQL 接触项在轻场上的投影。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

light = run_path(str(Path(__file__).with_name("ch39-boundary-light-yukawa.py")))
model = light["boundary_model"]
adjoint = light["adjoint"]


def flavor_bound(left, right, higgs):
    _, masses = light["squared_yukawa"](left, right, higgs)
    absolute = acb_mat(19,19)
    for i in range(19):
        absolute[i,i] = abs(higgs[i,0])

    def compressed_norm(profile):
        gram = adjoint(profile)*profile
        weighted = gram.solve(adjoint(profile)*absolute*profile)
        return max(v.real.upper() for v in weighted.eig(multiple=True))

    # Weighted Cauchy-Schwarz controls every complex site perturbation with
    # |delta_i| <= epsilon. Weyl then bounds both extreme singular values.
    factor = (compressed_norm(left)*compressed_norm(right)).sqrt()
    ratio = arb(1)/100
    return (masses[0]-ratio*masses[2])/((1+ratio)*factor)


def contact_bounds(profiles, boundary, height):
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    flat = light["extension"](0,height,boundary)
    for i in range(19):
        for j in range(3):
            reference = identity[boundary.index(i),j] if i in boundary else acb(1)/3
            assert abs(flat[i,j]-reference) < arb("1e-45")
    # Check the unit couplings from every interior site to every boundary site.
    # At zero flux this makes each interior harmonic measure exactly one third.
    matrix = model["laplacian"](0,height)
    assert all((matrix[i,j]+1).is_zero() for i in range(19) if i not in boundary for j in boundary)
    cap = arb(19)/3
    diagonal = 1/(3*cap*cap)
    delta = (2*(1-arb(".97"))).sqrt()
    mixed = (1/(2*cap).sqrt()-delta)/(3*cap*cap.sqrt())
    assert mixed > arb(".00074") and diagonal > arb(".00831")
    print("任意相位下相同夸克及重叠至少 0.97 的混合夸克投影下界", diagonal, mixed)
    # If xx*y vanishes at every boundary site, the two nonzero traces have
    # disjoint supports. Flat harmonic extensions of (1,0,0) and (0,1,1)
    # attain the maximum overlap for this constraint.
    a = flat*acb_mat([[1],[0],[0]])
    b = flat*acb_mat([[0],[1],[1]])
    overlap = abs((adjoint(a)*b)[0,0])/((adjoint(a)*a)[0,0].real*(adjoint(b)*b)[0,0].real).sqrt()
    threshold = 16*arb(2).sqrt()/(5*arb(41).sqrt())
    assert abs(overlap-threshold) < arb("1e-45")
    print("全轻轻子投影为零时两夸克波函数的尖锐最大重叠",threshold.str(12))

    omega = acb(0,2*arb.pi()/3).exp()
    def norm_squared(charge, k):
        f = acb_mat([[omega**(i*k)/arb(3).sqrt()] for i in range(3)])
        wave = profiles[charge]*f
        return sum(abs(wave[i,0])**2 for i in range(19))
    lepton = norm_squared(-3,0)
    for k in range(3):
        amplitude = 1/(3*(norm_squared(1,k)**3*lepton).sqrt())
        assert arb(".0515") < amplitude < arb(".144")
        print("三份同一夸克特征标与循环不变轻子的边界重叠",k,amplitude.str(12))


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(11431).imag
    boundary = model["boundary_sites"]()
    profiles = {q:light["extension"](q,height,boundary) for q in (1,-4,2,-3,6)}
    higgs = light["extension"](3,height,[model["soft"]["pole"]])
    higgs = higgs/sum(abs(higgs[i,0])**2 for i in range(19)).sqrt()
    for l,r,threshold in ((1,-4,".8355"),(1,2,".8214"),(-3,6,".7866")):
        eps = flavor_bound(profiles[l],profiles[r],higgs if r==-4 else higgs.conjugate())
        assert eps > arb(threshold)
        print("最小/最大 Yukawa 至多百分之一所需的相对站点变化",l,r,eps.str(12))
    contact_bounds(profiles,boundary,height)


if __name__ == "__main__":
    verify()
