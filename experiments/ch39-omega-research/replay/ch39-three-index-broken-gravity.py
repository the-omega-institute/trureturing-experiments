"""固定五个正截止节点，在相邻三个实际零点共享电弱背景引力权。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def gravity_columns(root, numbers, cap=1):
    single, lib, fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"](root,numbers)
    r = (arb("246.22")/arb("1.22089e19"))**2/2
    yt, g1, g2 = map(arb, (".8616", ".467453", ".63811"))
    kw, kz = g2**2/2, (g2**2+3*g1**2/5)/2
    light = [(-59,arb(2)/3,arb(0)), (-12,1,yt**2*r), (6,-1,kw*r), (3,-arb(1)/2,kz*r)]
    sites = list(map(arb, ("4e-15","3.5e-11","1.23e-6",".5","1")))
    all_kw = []
    for rows in fields:
        moved = light+rows[1:49]+rows[65:]
        for g,a,z in rows[49:65]:
            moved.extend(((g-4,a-arb(2)/3,z), (4,arb(2)/3,z+(-cap*r).union(cap*r))))
        all_kw.append(arb_mat([single["moments"](moved,lib["expint"],x) for x in sites]).transpose())
    return all_kw, sites, fields[0][49:65], lib["expint"]


def verify():
    n = 2_306_888_100_085
    root = Path(__file__).resolve().parent
    spectra, _, _, _ = gravity_columns(root, (n,n+1,n-1))
    target = arb("1e-122")
    for fixed in (0,1):
        varied = 1-fixed
        a = spectra[0]
        diff = [b-a for b in spectra[1:]]
        constraints = [[1]*5,[a[fixed,j] for j in range(5)]]
        constraints += [[b[fixed,j] for j in range(5)] for b in diff]
        mat = arb_mat(constraints+[[a[varied,j] for j in range(5)]])
        print("固定读数、行列式", fixed,mat.det().str(12),flush=True)
        assert not mat.det().contains(0)
        rhs = arb_mat([[1],[1 if fixed==0 else target],[0],[0],[target if varied==1 else 1]])
        zero = mat.solve(rhs)
        slope = mat.solve(arb_mat([[0],[0],[0],[0],[1]]))
        # t 是中心编号上未固定读数相对目标的偏差。
        def response(b,col):
            return sum(b[varied,j]*col[j,0] for j in range(5))
        offsets = [arb(0)]+[response(b,zero) for b in diff]
        rates = [arb(1)]+[1+response(b,slope) for b in diff]
        print("三条仿射响应",*[v.str(16) for pair in zip(offsets,rates) for v in pair],flush=True)
        # 对于任意 t，两行的有符号线性组合消去 t，给出最大偏差下界。
        i,j = 1,2
        assert rates[i] > 0 and rates[j] > 0
        bound = abs(rates[j]*offsets[i]-rates[i]*offsets[j])/(rates[i]+rates[j])
        t = -(offsets[i]+offsets[j])/(rates[i]+rates[j])
        errors = [b+s*t for b,s in zip(offsets,rates,strict=True)]
        assert abs(errors[0]) < bound
        if fixed==0:
            assert arb("1.85650964e-8") < bound < arb("1.85650966e-8")
        else:
            assert arb(".00122900360") < bound < arb(".00122900361")
        print("邻点对给出的锐下界、中心误差、达成三项误差",bound.str(12),
              t.str(12),*[v.str(12) for v in errors],flush=True)
        w = zero+slope*t
        print("最优五权",*[v.str(10) for v in w.entries()],flush=True)


if __name__=="__main__":
    verify()
