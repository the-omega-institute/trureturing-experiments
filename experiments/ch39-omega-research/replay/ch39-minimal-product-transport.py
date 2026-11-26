"""十九点壳上三条腿各自作一次独立输运时，所需乘积顶角项数的下界与见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx

study = run_path(str(Path(__file__).with_name("ch39-three-leg-yukawa.py")))
masks = study["masks"]

SECOND_WITNESS = (
    ("16.41466301022515", "-45.80705362235868"),
    ("-7.057476964906035", "155.15532079553898"),
    ("22.27096654945317", "-78.39661555773696"),
    ("-1399.3822679241896", "3088.6066482527713"),
    ("7.038282244951122", "-129.56456086865026"),
    ("-1070.6332933976212", "2773.762613663774"),
    ("1415.5858979836041", "-2945.45849283614"),
    ("-1683.038021853561", "4175.417182457006"),
)
SECOND_HIGGS = (
    ("0.9999999921065358", "0.0"),
    ("2.0687997963158922e-05", "0.0001234671374720156"),
    ("6.8750214778259555e-06", "-8.217996414874848e-06"),
)


def lift_groups(leg):
    # w⊗w 的相同单项式共享一列。leg=None 时三条腿均独立。
    groups = {}
    for i,a in enumerate(masks):
        for j,b in enumerate(masks):
            if leg is None:
                key = tuple(x+y for x,y in zip(a,b))
            else:
                rest = tuple(sorted(x[:leg]+x[leg+1:] for x in (a,b)))
                key = a[leg]+b[leg], rest
            groups.setdefault(key, []).append((i,j))
    assert len(groups) == (27 if leg is None else 30)
    return list(groups.values())


def lifted(form, groups):
    # P†(form⊗I)P，只求和第二个索引相同的项。
    return acb_mat([[sum(form[i,j] for i,a in left for j,b in right if a==b)
                     for right in groups] for left in groups])


def rank_two(w):
    # 两切片 A,B。BA⁻¹ 有两个不同本征值，谱投影分解给两个秩一切片。
    first, second = [acb_mat([[w[4*a+2*b+c,0] for c in range(2)] for b in range(2)])
                     for a in range(2)]
    assert abs(first.det()) > 0
    ratio = second*first.inv()
    discriminant = (ratio[0,0]+ratio[1,1])**2 - 4*ratio.det()
    assert abs(discriminant) > 0


def verify():
    refs = run_path(str(Path(__file__).with_name("ch39-shell-flavor-matching-windows.py")))
    mixing = run_path(str(Path(__file__).with_name("ch39-finite-orthogonal-ckm.py")))
    ctx.prec = 448
    target = [y for g in ("up","down","leptons") for y in refs[g]]
    sites = (22027,28325)
    families = ((None,".267"),(0,".583"),(1,".906"),(2,".456"))
    for n in sites:
        blocks, _, _, _ = study["tensors"](n)
        for leg, limit in families:
            groups = lift_groups(leg)
            for block in blocks:
                up, denom = study["gram_pair"](block)
                study["positive"](lifted(arb(limit)**2*denom - up, groups))
        first = study["column"](study["WITNESS"])
        second = study["column"](SECOND_WITNESS)
        w = first if n == 22027 else second
        rank_two(w)
        z = study["column"](study["HIGGS"] if n == 22027 else SECOND_HIGGS)
        z = z/(study["adjoint"](z)*z)[0,0].real.sqrt()
        matrices = study["mass_matrices"](blocks,w,z)
        masses = [s for y in matrices for s in sorted(
            (r.real.sqrt() for r in (y*study["adjoint"](y)).eig(algorithm="rump")), key=float)]
        error = max(abs(y/ref-1).upper() for y,ref in zip(masses,target))
        assert error < arb("3.3e-8")
        angles = mixing["observables"](*matrices[:2])[4:]
        goals = [arb(x).sin() for x in (".22704",".04275",".00377")]+[arb("1.139")]
        assert all(abs(v-ref)<arb("1e-8") for v,ref in zip(angles[:4],goals))
        print(n, "两乘积项的质量最大相对误差、三个混合角正弦和 CP 相位",
              error.str(9), *(a.str(9) for a in angles[:4]))
        other = second if n == 22027 else first
        relative = []
        for block in blocks:
            u,d = study["gram_pair"](block)
            relative.append(((study["adjoint"](other)*u*other)[0,0].real /
                             (study["adjoint"](other)*d*other)[0,0].real).sqrt())
        cap = max(v.upper() for v in relative)
        assert cap < arb(".571" if n == 22027 else ".757")
        print(n,"用另一枚系数时，联合质量判据上界",cap.str(9))


if __name__ == "__main__":
    verify()
