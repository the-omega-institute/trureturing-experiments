"""十九层的非负非降味响应与第三代 CKM 元的全域不相容性。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify(root):
    family = run_path(str(root/"ch39-three-gram-source-flavor.py"))
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 320
    frame, _ = family["family_space"](native["gram_cholesky"](root))
    adj = family["adj"]
    positions = [shell["soft"]["pole"], *shell["boundary_sites"]()]
    assert positions == [14,8,7,6]
    height = acb.zeta_zero(2_306_888_100_085).imag
    projected, lengths = [], []
    for charge in (-4,2):
        levels, vectors = shell["laplacian"](charge,height).eig(right=True,algorithm="rump")
        ordered = sorted(range(19),key=lambda i:float(levels[i].real.mid()))
        assert all(levels[i].imag.contains(0) for i in ordered)
        assert levels[ordered[-2]].real < levels[ordered[-1]].real
        assert all(levels[i].real<levels[ordered[-1]].real for i in ordered[:-1])
        highest = ordered[-1]
        norm = sum(abs(vectors[i,highest])**2 for i in range(19)).sqrt()
        four = acb_mat([[vectors[i,highest]/norm] for i in positions])
        v = adj(frame)*four
        g = (adj(v)*v)[0,0].real.sqrt()
        assert g>0
        projected.append(v*(1/g))
        lengths.append(g)
    cosine = abs((adj(projected[0])*projected[1])[0,0])
    angle = cosine.acos()
    radius = arb("1.001")/arb(".999")
    mass_ratios = (radius*arb(".00311")/arb(".8616"),radius*arb(".000268")/arb(".01401"))
    u,d = [r/g for r,g in zip(mass_ratios,lengths,strict=True)]
    physical_gap = angle-u.asin()-d.asin()
    assert physical_gap>0
    ceiling = physical_gap.cos()
    s23 = ckm["centers"][1]+ckm["sine_radius"]
    s13 = ckm["centers"][2]+ckm["sine_radius"]
    reference = ((1-s23**2)*(1-s13**2)).sqrt()
    assert ceiling < arb(".597803") < arb(".999115") < reference
    print("两荷最高投影模长、归一投影夹角及内积模",lengths[0].str(14),lengths[1].str(14),
          angle.str(14),cosine.str(14))
    print("两个质量比的上界、顶底角的下界、|Vtb|上界、参考下界",
          *[v.str(12) for v in (*mass_ratios,physical_gap,ceiling,reference)])
    # 即使谱响应非单调，也可用末层相对该部门最大读数的比例衡量必需的抑制。
    delta = angle-reference.acos()
    threshold = (u**2+d**2+2*u*d*delta.cos())/delta.sin()**2
    assert delta>0 and delta<arb.pi()/2
    assert max(u,d)**2<threshold<arb(".001174")
    for value, other in ((u,d),(d,u)):
        limit = (other/(delta-value.asin()).sin())**2
        assert limit < (arb(".000932") if value is u else arb(".0000351"))
    print("两末层/峰值比例中较小者的必要上界",threshold.str(14))
    # 若共同三代投影改变 eps，则每条已投影壳向量的位移至多 eps。
    eps = arb(".331")
    tilts = [((eps/g).asin() + (r/(g-eps)).asin())
             for r,g in zip(mass_ratios,lengths,strict=True)]
    gap = angle - sum(tilts)
    assert gap>0 and gap.cos()<arb(".99904")<reference
    print("共同味投影算子范数移动至多 0.331 后 |Vtb| 上界",gap.cos().upper().str(14))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
