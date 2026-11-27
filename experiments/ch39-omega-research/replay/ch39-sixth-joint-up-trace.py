"""在第六壳联合输运六条上型谱线，排除末步射影距离 0.0700 内的上夸克窗。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


LAYERS = (1, 7, 10, 12, 16, 18)
SITES = (14, 8, 7, 6)


def adjoint(matrix):
    return matrix.conjugate().transpose()


def trace_directions(overlap, spectrum):
    count = overlap.nrows()
    weighted = acb_mat([[overlap[i, j]*value for j, value in enumerate(spectrum)] for i in range(count)])
    compression = weighted*adjoint(overlap)
    # 中点矩阵只选择试探方向。随后正交化精确的二进有理列，并在原盒上估值。
    centre = acb_mat([[compression[i, j].mid() for j in range(count)] for i in range(count)])
    eigenvalues, candidates = centre.eig(right=True, algorithm="rump")
    indices = sorted(range(count), key=lambda j: eigenvalues[j].real.mid(), reverse=True)[:3]
    directions = []
    for j in indices:
        vector = acb_mat([[candidates[i, j].mid()] for i in range(count)])
        for q in directions:
            vector = vector-q*(adjoint(q)*vector)[0, 0]
        norm = (adjoint(vector)*vector)[0, 0].real.sqrt()
        assert norm > 0
        directions.append(vector/norm)
    return acb_mat([[q[i, 0].conjugate() for i in range(count)] for q in directions])


def row_floors(overlap, resolvents, shifts, fidelity):
    prob = arb_mat([[abs(z)**2 for z in row] for row in overlap.tolist()])
    den = prob*resolvents
    return [max([arb(0)]+[(fidelity/den[i, j]-s).lower() for j, s in enumerate(shifts)])
            for i in range(overlap.nrows())]


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    w, dw, sectors, _, _ = lib("fifth-third-pivot-joint")["verify"](root, return_transport=True)
    ctx.prec = 384
    assert dw < arb("2.81e-9")
    shell = lib("three-ray-chiral-boundary")
    assert (shell["soft"]["pole"], *shell["boundary_sites"]()) == SITES
    dst = lib("adjacent-pvm-transport")["spectrum"](shell, -4, 2_306_888_100_091)[1]
    src, _, spectrum, _, _ = sectors[0]
    overlap = adjoint(dst)*adjoint(w**4)*src
    distance = 2*(arb(".0700")+2*dw)
    assert distance < 1 and min(spectrum) > 0
    fidelity = (1-distance**2/2)**2
    shifts = [arb(2)**k*m/4 for k in range(-38, 10) for m in (4, 5, 6, 7)]
    resolvents = arb_mat([[1/(v+s) for s in shifts] for v in spectrum])
    lower = row_floors(overlap, resolvents, shifts, fidelity)
    selected = acb_mat([overlap.tolist()[i-1] for i in LAYERS])
    test_rows = trace_directions(selected, spectrum)*selected
    joint = sum(row_floors(test_rows, resolvents, shifts, fidelity))
    excess = (joint-sum(lower[i-1] for i in LAYERS)).lower()
    assert excess > arb(".00660")
    four = acb_mat([[dst[i, j] for j in range(dst.ncols())] for i in SITES])
    base = lib("common-unitary-quark-right-trace")["weighted"](four, lower)
    floor = arb(".00000633")
    identity = acb_mat([[int(i == j) for j in range(4)] for i in range(4)])
    for i in LAYERS:
        col = acb_mat([[four[k, i-1]] for k in range(4)])
        test = base+excess*col*adjoint(col)-floor*identity
        for k in range(1, 5):
            minor = acb_mat([[test[a, b] for b in range(k)] for a in range(k)]).det()
            assert minor.imag.contains(0) and minor.real > 0, (i, k, minor)
    assert floor > arb(".00000629")
    print("六层超额、六个顶点的统一严格质量下界、原上窗窗顶",
          excess.str(14), floor.str(9), ".00000629", flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
