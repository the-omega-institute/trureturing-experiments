"""从完整训练盒生成第五壳轻夸克扰动证书的有理输入表。"""
import argparse
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def build_table(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    z, error, sectors, _, _ = lib("fourth-third-pivot-joint")["verify"](root, return_transport=True)
    ctx.prec = 384
    assert error < arb("1.28e-10")
    scale = 10**14
    def ceil(value):
        return int((scale*value).upper().ceil().unique_fmpz())
    def floor(value):
        return int((scale*value).lower().floor().unique_fmpz())
    def grid(value):
        q = int((scale*value.mid()+arb(1)/2).floor().unique_fmpz())
        assert abs(scale*value-q)<1
        return q

    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    adj = lib("rank-two-generator-up-quark")["adj"]
    sites = (shell["soft"]["pole"], *shell["boundary_sites"]())
    assert sites==(14,8,7,6)
    rows = [("scale",scale),("shell",2_306_888_100_090),("sites",*sites)]
    for sector, (charge, op) in enumerate(((-4,adj(z**4)),(2,z**2))):
        src, _, training, _, _ = sectors[sector]
        dst = spectrum(shell, charge, 2_306_888_100_090)[1]
        # 用每列中第一个模大于 0.1 的分量固定该列的相位。
        phases = []
        for j in range(19):
            k = next(i for i in range(19) if abs(dst[i,j])>arb(".1"))
            assert all(abs(dst[i,j])<arb(".1") for i in range(k))
            v = dst[k,j]
            phases.append(v.conjugate()/abs(v))
        dst = acb_mat([[dst[i,j]*phases[j] for j in range(19)] for i in range(19)])
        ov = adj(dst)*op*src
        rows.append(("train",charge,*map(floor,training)))
        for i in range(19):
            prob = [abs(ov[i,j])**2 for j in range(19)]
            assert sum(prob).contains(1)
            rows.append(("prob",charge,i,*map(ceil,prob)))
        for k in sites:
            rc = [a for j in range(19) for a in (grid(dst[k,j].real),grid(dst[k,j].imag))]
            rows.append(("basis",charge,k,*rc))
    return "".join("\t".join(map(str,row))+"\n" for row in rows)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true",help="由原始训练盒写出有限输入表")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    path = root/"ch39-fifth-quark-distance-cells.tsv"
    data = build_table(root)
    if args.write:
        path.write_text(data)
    else:
        assert path.read_text()==data
    print("有理输入表逐项认证完毕")
