"""按 02823ab9 的协议，在第九壳读取已由第七、八壳联合选定的同一个算子。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    star = lib("two-shell-symmetric-star-gates")
    box, _ = star["experiment"](root, fifth, columns, levels)
    (geometry,) = star["prepare"](root, lib, fifth, columns, levels, (2_306_888_100_094,))
    fixed, blocks, matching, weights, heat, field = geometry
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    midpoint = [q.mid() for q in box]
    center, _ = gravity(*star["project"](midpoint, fixed, blocks, matching),
                        levels, weights, heat, field, True)
    _, jac = gravity(*star["project"](box, fixed, blocks, matching),
                     levels, weights, heat, field, False)
    # 中心取有质量原核；同一原盒上的导数把所有联合真根传播到第九壳。
    values = [base + center[i,0] + sum(jac[i,j]*(box[j]-q) for j,q in enumerate(midpoint))
              for i,base in enumerate((arb(1),arb("1e-122")))]
    print("第九壳的牛顿、真空项及其有向端点",
          *[(v.str(10),v.lower().str(15),v.upper().str(15)) for v in values],flush=True)
    k,w = values
    assert arb("-1.440750") < k < arb("-1.433173")
    assert arb("-.016891") < w < arb("-.016474")
    return values


def verify(root):
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
