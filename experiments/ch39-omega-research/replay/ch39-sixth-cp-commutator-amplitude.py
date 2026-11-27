"""从第六壳实际四站块计算共同三代平面上 CP 交换子行列式的精确极值区间。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def adj(h):
    return h.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    (tu, td), ((ua, ub), (da, db)) = lib("sixth-heavy-frontier")["blocks"](root)
    eye = acb_mat([[1, 0], [0, 1]])
    h = (tu-eye*ua)*(da-db)+(td-eye*da)*(ub-ua)
    comm = tu*td-td*tu
    # 无迹二阶反厄米矩阵的两个奇异值均为 sqrt(det(comm))。
    determinant = h.det()
    assert determinant.imag.contains(0) and comm.det().imag.contains(0)
    norm = comm.det().real.sqrt()
    amplitude = abs(determinant.real)*norm/27
    assert arb("1.15448710e-10") < amplitude < arb("1.15455756e-10")
    # 两个极值法向量直接压回四站矩阵，独立核验两种符号。
    for sign in (-1, 1):
        eigproj = (eye+comm*(acb(0, sign)/norm))/2
        v = [eigproj[i, 0]/eigproj[0, 0].real.sqrt() for i in range(2)]
        sqrt2, sqrt6 = arb(2).sqrt(), arb(6).sqrt()
        frame = acb_mat([[-v[1].conjugate(), 2*v[0]/sqrt6, 0],
                         [v[0].conjugate(), 2*v[1]/sqrt6, 0],
                         [0, -1/sqrt6, 1/sqrt2], [0, -1/sqrt6, -1/sqrt2]])
        compressed = []
        for block, a, b in ((tu, ua, ub), (td, da, db)):
            raw = acb_mat([[block[0, 0], block[0, 1], 0, 0],
                           [block[1, 0], block[1, 1], 0, 0], [0, 0, a, 0], [0, 0, 0, b]])
            compressed.append(adj(frame)*raw*frame)
        left, right = compressed
        direct = (left*right-right*left).det()
        assert direct.real.contains(0) and direct.imag.overlaps(sign*amplitude)
        assert sign*direct.imag > 0
        print("极值平面直接压缩得到的虚部", sign, direct.imag.str(10), flush=True)
    mix = acb_mat(lib("shell-clock-ckm-blocks")["standard_matrix"](*[arb(v).sin() for v in
                  (".22704", ".04275", ".00377")], arb("1.139")))
    jarlskog = (mix[0, 0]*mix[1, 1]*mix[0, 1].conjugate()*mix[1, 0].conjugate()).imag
    mass = ((".00000615", ".00311", ".8616"), (".0000135", ".000268", ".01401"))
    diag = lib("complex-four-site-unit-quark-amplitudes")["diag"]
    target_up = adj(mix)*diag(list(map(arb, mass[0])))*mix
    target_down = diag(list(map(arb, mass[1])))
    target = (target_up*target_down-target_down*target_up).det()
    vandermonde = arb(1)
    for row in mass:
        a, b, c = map(arb, row)
        vandermonde *= (b-a)*(c-b)*(c-a)
    assert target.real.contains(0) and abs(target.imag).overlaps(2*jarlskog*vandermonde)
    assert arb("7.20658690126e-15") < target.imag < arb("7.20658690127e-15")
    assert arb("16019.88728") < amplitude/target.imag < arb("16020.86506")
    print("det D、平凡交换子范数、全平面锐幅、中心目标虚部、锐幅/目标绝对值",
          *[(v.lower().str(14), v.upper().str(14)) for v in
            (determinant.real, norm, amplitude, target.imag, amplitude/abs(target.imag))], flush=True)
    return amplitude, target.imag


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
