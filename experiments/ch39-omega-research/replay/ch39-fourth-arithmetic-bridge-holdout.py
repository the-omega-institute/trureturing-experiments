"""从两步联合证书验证三条算术源所选秩一桥的第四壳读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root, rule=None, *, charge_steps=None):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    x, y, incoming, readout, cutoffs, pivot, thermal = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    if charge_steps is None:
        step = lib("fourth-arithmetic-rank-one-bridge")["verify"](root) if rule is None else rule(x,y)
        up, down, gravity = adj(step**2), step, step**3
    else:
        up, down, gravity = charge_steps(root, y)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    _, labels, _ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root, ((n+1, ".0294"),), ".1", (), ((6,7,10,13), (6,7,10,12)))
    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    sites = (shell["soft"]["pole"], *shell["boundary_sites"]())
    assert sites == (14,8,7,6)
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    transports = (up*adj(y**2)*adj(x**2), down*y*x)
    for charge, op, limits in zip((-4,2), transports, windows, strict=True):
        source, target = [spectrum(shell, charge, i)[1] for i in (n+1,n+4)]
        amplitudes = rank["read"](source, target, labels[n+1,charge], op)
        four = acb_mat([[target[i,j] for j in range(19)] for i in sites])
        matrix = acb_mat([[four[i,j]*amplitudes[j] for j in range(19)] for i in range(4)])*adj(four)
        levels = sorted((z.real for z in matrix.eig(algorithm="rump")), key=lambda s:s.mid())
        assert all(a < b for a,b in zip(levels,levels[1:]))
        verdict = []
        for i,(lo,hi) in enumerate(limits):
            if levels[i] > arb(hi) or levels[i+1] < arb(lo):
                verdict.append("排除")
            else:
                assert levels[i] < arb(lo) < arb(hi) < levels[i+1]
                verdict.append("整窗可达")
        print(charge, [v.str(14) for v in levels], verdict, flush=True)

    low, high, highest = cutoffs
    moved = readout * gravity * incoming
    amplitudes = [low + (high-low)*abs(moved[i,0])**2 + (highest-low)*abs(moved[i,1])**2 for i in range(19)]
    moments, heat, field = thermal
    p0,p1 = [(arb(s)/arb("1.22089e19"))**2 for s in ("91.1876","1.176397e14")]
    assert all(a > 0 and a*a > p0 for a in amplitudes)
    kernels = [moments(field,heat,a*a) for a in amplitudes]
    weak = lib("full-neutral-free-cutoff")["weak_column"]
    columns = (*zip(*kernels), [weak(a*a,p0,p1) for a in amplitudes])
    values = [sum(u*v for u,v in zip(pivot,col,strict=True)) for col in columns]
    assert all(not val.overlaps(arb(goal)) for val,goal in zip(values,("1","1e-122","0"),strict=True))
    print("荷六 K、W、D", [v.str(15) for v in values], flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
