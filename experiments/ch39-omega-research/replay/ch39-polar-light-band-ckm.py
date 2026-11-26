"""把原壳最低三条本征线用投影的极分解等距送至同一算术三源空间。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    four = run_path(str(root / "ch39-four-gram-triplet-ckm.py"))
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 320
    gram = native["gram_cholesky"](root)
    adj = flavor["adj"]
    sites = [shell["soft"]["pole"], *shell["boundary_sites"]()]
    center = 2_306_888_100_085
    for offset in (0, -1, 1):
        height = acb.zeta_zero(center + offset).imag
        restricted = []
        for charge in (-4, 2):
            levels, vectors = shell["laplacian"](charge, height).eig(right=True, algorithm="rump")
            order = sorted(range(19), key=lambda j: float(levels[j].real.mid()))
            assert all(v.imag.contains(0) for v in levels)
            assert all(levels[i].real < levels[j].real for i, j in zip(order, order[1:]))
            light = order[:3]
            norms = [sum(abs(vectors[k,j])**2 for k in range(19)).sqrt() for j in light]
            restricted.append(acb_mat([[vectors[i,j]/n for j,n in zip(light,norms,strict=True)]
                                       for i in sites]))
        for omit in range(4):
            normal = acb_mat([[v] for v in four["normal"](gram, omit)])
            projection = flavor["identity"](4) - normal*adj(normal)
            transported = []
            for low in restricted:
                z = projection * low
                eigenvalues, subspaces = flavor["projections"](adj(z) * z)
                inv_sqrt = sum((p * (1/v.sqrt()) for v,p in zip(eigenvalues,subspaces,strict=True)),
                               acb_mat(3,3))
                q = z * inv_sqrt
                assert (adj(q)*q).contains(flavor["identity"](3))
                transported.append(q)
            v = adj(transported[0]) * transported[1]
            assert (adj(v)*v).contains(flavor["identity"](3))
            mod = [abs(x) for x in v.entries()]
            j = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
            nearest = min(abs(t-ckm["centers"][2]).lower() for t in mod)
            print("编号偏移、省略源、按能序的 |Vub|、J、任意重排距 |Vub| 中心的下界",
                  offset,omit,mod[2].str(12),j.str(12),nearest.str(12),flush=True)
            assert min(x.lower() for x in mod) > arb(".025")
            assert nearest > arb(".02") > ckm["sine_radius"]
            assert (j > 0 if offset == -1 else j < 0)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
