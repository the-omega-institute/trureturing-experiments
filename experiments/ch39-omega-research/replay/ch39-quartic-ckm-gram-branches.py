"""从质量留数与三阶实 Gram 行列式计算两个四次见证的带号 J。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat


def verify(root):
    quartic = run_path(str(root/"ch39-common-quartic-six-masses.py"))
    identity = acb_mat([[1, 0], [0, 1]])
    for number, coeff, box, nodes, rows, masses, _, _ in quartic["certified_cases"](root, radius_text="1e-25"):
        residues, _, operators = quartic["read"](coeff, box, nodes, rows, masses, matrices=True)
        blocks = [acb_mat([[op[i, j] for j in range(2)] for i in range(2)]) for op in operators]
        (ua, ub), (da, db) = [(op[2, 2].real, op[3, 3].real) for op in operators]
        det = ((blocks[0]-identity*ua)*(da-db)+(blocks[1]-identity*da)*(ub-ua)).det()
        assert det.imag.contains(0) and not det.real.contains(0)
        rho = [(r[0]+r[1]) for r in residues]
        pa, pb = [residues[0][i].intersection(residues[1][i]) for i in (2, 3)]
        length = rho[0].intersection(rho[1])
        assert pa > 0 and pb > 0 and (pa+pb+length).contains(1)
        pauli = [[b[0, 1].real, -b[0, 1].imag, (b[0, 0].real-b[1, 1].real)/2] for b in blocks]
        uv = [[sum(x*y for x, y in zip(v, w, strict=True)) for w in pauli] for v in pauli]
        projection = [uv[i][i].sqrt()*(r[1]-r[0]) for i, r in enumerate(residues)]
        gram = arb_mat([[*uv[i], projection[i]] for i in range(2)]+[[*projection, length**2]]).det()
        assert gram > 0
        vandermonde = arb(1)
        for a, b, c in masses:
            vandermonde *= (b-a)*(c-b)*(c-a)
        # 镜像正支指向 a×b，一般因子式和目标的递增谱序固定 J 的符号。
        positive_orientation = -pa*pb*det.real*gram.sqrt()/vandermonde
        expected = {2_306_888_100_085: ("-.000531999235997127", "-.000531999235997125"),
                    2_306_888_100_086: (".0000122467675281111", ".0000122467675281113")}
        low, high = map(arb, expected[number])
        assert low < positive_orientation < high
        print("编号、Gram 行列式、det D、负与正镜像支的 J", number,
              *[v.str(15) for v in (gram, det.real, -positive_orientation, positive_orientation)], flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
