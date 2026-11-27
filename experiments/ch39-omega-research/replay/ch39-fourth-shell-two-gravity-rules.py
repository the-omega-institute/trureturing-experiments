"""从两壳联合证书全相干延续到 N+4，检验两种第三步的荷六读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def shell_data(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    upstream = lib("one-generator-gravity-full-ckm")
    joint = lib("next-common-ckm-gravity")
    cx, rx, cy, ry, _ = joint["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root, refine_node=n+1)[2][1]
    assert all(v is ordered[0] for v in ordered[:17])
    low, high, highest = [ordered[j].sqrt() for j in (0,17,18)]
    assert 0 < low < high < highest

    moments, heat, fields = lib("three-index-free-atoms")["broken_fields"](root, [n+4])
    native = lib("native-rh-ray-pvm")
    _, pivot_levels, crosses, _ = next(row for row in native["shell_weights"](root, n+4) if row[0] == 6)
    weights = [w.real for w in crosses[0]]
    assert sum(weights).contains(1)
    ctx.prec = 384
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    shell = lib("three-ray-chiral-boundary")
    _, source = spectrum(shell, 6, n+1)
    levels, target = spectrum(shell, 6, n+4)
    assert all(a.overlaps(b) for a,b in zip(levels, pivot_levels, strict=True))
    source = acb_mat([[source[i,17], source[i,18]] for i in range(19)])

    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    local = lib("next-common-ckm")["local"]
    eye = rank["eye"](19)
    def generator(center, radius, directions, name):
        box = [v + arb(0,radius.upper()) for v in center]
        return local(seed(root, rank["eye"], name), box, directions, eye, False)[0]

    x = generator(cx, rx, upstream["DIRECTIONS"], "ch39-full-ckm-gravity-cayley.tsv")
    y = generator(cy, ry, joint["DIRECTIONS"], "ch39-next-common-ckm-gravity-cayley.tsv")
    incoming = y**3 * x**3 * source
    readout = target.conjugate().transpose()
    return x, y, incoming, readout, (low,high,highest), weights, (moments,heat,fields[0])


def verify(root):
    x,y,incoming,readout,(low,high,highest),weights,(moments,heat,field) = shell_data(root)
    z = y * x.conjugate().transpose() * y
    p0, p1 = [(arb(s)/arb("1.22089e19"))**2 for s in ("91.1876", "1.176397e14")]
    weak = run_path(str(root / "ch39-full-neutral-free-cutoff.py"))["weak_column"]
    goals = (arb(1), arb("1e-122"), arb(0))
    results = {}
    for name, step in (("YX*Y",z), ("Y",y)):
        ov = readout * step**3 * incoming
        amplitudes = [low + (high-low)*abs(ov[i,0])**2 + (highest-low)*abs(ov[i,1])**2
                      for i in range(19)]
        assert all(v>0 and v*v>p0 for v in amplitudes)
        thermal = [moments(field, heat, v*v) for v in amplitudes]
        columns = (*zip(*thermal), [weak(v*v, p0, p1) for v in amplitudes])
        result = [sum(w*b for w,b in zip(weights, col, strict=True)) for col in columns]
        print("第三步、最小振幅下界、最大振幅上界、K、W、D", name,
              min(v.lower() for v in amplitudes).str(15),max(v.upper() for v in amplitudes).str(15),
              *[v.str(18) for v in result], flush=True)
        assert all(not v.overlaps(goal) for v,goal in zip(result,goals,strict=True))
        results[name] = result
    return results


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
