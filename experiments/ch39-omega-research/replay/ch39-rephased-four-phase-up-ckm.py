"""四相位构造的上型读数，在所有共同三平面上与中心 CKM 的必要矩相比较。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def verify(root):
    construct = run_path(str(root / "ch39-rephased-involution-up-centers.py"))
    witness = construct["verify"](root, return_matrices=True)
    adj = construct["adj"]
    four_bounds = witness["levels"]
    model = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    data = run_path(str(root / "ch39-up-unrestricted-character-floor.py"))
    labels, target = data["shell_data"](root)
    _, source = spectrum(shell, -4, model["N"]+1)
    center = witness["up_center"]
    err = witness["up_error"]
    native = model["read"](source, target, labels, center)
    sites = [14, 8, 7, 6]
    cells = acb_mat([[target[i,j] for j in range(19)] for i in sites])
    weighted = acb_mat([[cells[i,j] * native[j] for j in range(19)] for i in range(4)])
    raw = weighted * adj(cells)
    omega = acb(-1,arb(3).sqrt())/2
    characters = []
    for z in (omega, omega.conjugate()):
        vector = acb_mat([[0], [1], [z], [z**2]])*(1/arb(3).sqrt())
        mean = (adj(vector)*raw*vector)[0,0].real
        delta = max(a.upper() for a in labels).sqrt()*err
        assert mean.sqrt()>delta
        bound = ((mean.sqrt()-delta)**2).union((mean.sqrt()+delta)**2)
        j, = [j for j,b in enumerate(four_bounds) if b.overlaps(bound)]
        characters.append(j)
    assert len(set(characters))==2
    assert characters == [1, 3]
    reference = [four_bounds[j] for j in (1, 3, 0, 2)]
    lib = run_path(str(root / "ch39-fixed-reflection-quark-spectral-neighbourhood.py"))
    pencil = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    lines, geom, mixing, masses, a, b, _ = pencil["setup"](root)
    probabilities = lib["probability_box"](reference, arb(0), masses[0])
    print("共同弱轴投到上型三条质量本征线的概率", [p.str(12) for p in probabilities], flush=True)
    for first, second in ((a, b), (b, a)):
        count, depth, gap = lib["exclude_chart"](lines, geom, mixing, first, second, probabilities)
        print("铅笔坐标图的区间数、最大深度、概率距离下界", count, depth, gap.str(12), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
