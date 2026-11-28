"""四个单位置例外使匹配两个引力矩的强逆耦合在整段能标上为正。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def verify():
    def load(name):
        return run_path(str(ROOT / f"ch39-{name}.py"))

    edges, rows, heat, _ = load("two-copy-gram-cutoff-two-exceptions")["model_inputs"]()
    moments = load("signed-cutoff-gravity")["moments"]
    gauge = load("signed-gravity-unification")
    model, _, _, initial, _ = load("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()
    fields = gauge["gauge_weights"](model)
    fixed = {16: 1, 17: 1, 18: 1, 323: 0}
    pi = arb.pi()

    def evaluate(c, t, *, derivative=False, kinetic=False):
        out = [arb(0), arb(0)]
        jac = arb_mat(2, 2)
        at_z = list(initial)
        for i, (a, b, weight) in enumerate(edges):
            u = fixed.get(i, t)
            v = (1-u)*a+u*b
            x = c*v
            if derivative:
                rates = [(-z/x).exp() for _, _, z in rows]
                dk = weight*sum(alpha*y for (_,alpha,_),y in zip(rows,rates,strict=True))/(2*pi)
                dw = -weight*x*sum(g*y for (g,_,_),y in zip(rows,rates,strict=True))/(4*pi)
                for j, val in enumerate((dk,dw)):
                    jac[j,0] += val*v
                    if i not in fixed:
                        jac[j,1] += val*c*(b-a)
            elif kinetic:
                js = [heat(1,z/x) for _,_,z in rows[1:]]
                for j in range(3):
                    at_z[j] -= weight*sum(row[j]*h for row,h in zip(fields,js,strict=True))/(4*pi)
            else:
                k,w = moments(rows,heat,x)
                out[0] += weight*k
                out[1] += weight*w
        if derivative:
            return jac
        if kinetic:
            return at_z
        return arb_mat([[out[0]-1],[out[1]]])

    center = (arb(".019541502206").mid(), arb(".437018430143").mid())
    radius = arb(2)**-18
    scale, mix = [q+arb(0,radius) for q in center]
    assert 0 < scale < 1 and 0 < mix < 1
    inverse = arb_mat([[arb(".01956").mid(),arb(".51398").mid()],
                       [arb("-.000003896").mid(),arb("3.11636").mid()]])
    assert not inverse.det().contains(0)
    defect = arb_mat([[1,0],[0,1]])-inverse*evaluate(scale,mix,derivative=True)
    contraction = max(sum(abs(defect[i,j]) for j in range(2)).upper() for i in range(2))
    tau = arb(0,arb("5e-7").upper())
    shift = inverse*(evaluate(*center)-arb_mat([[0],[tau]]))
    move = max(abs(z).upper() for z in shift.entries())
    assert contraction < arb(".33") and move < arb(".00000156")
    assert move+contraction*radius < radius
    assert tau.contains(arb("1e-122"))
    intercept = evaluate(scale,mix,kinetic=True)
    assert model["beta"][0] > 0 and model["beta"][2] == -7
    planck = intercept[2]+7*(arb("1.22089e19")/model["mz"]).log()/(2*pi)
    assert intercept[0] < 0 and intercept[2] > 2687 and planck < 2747
    print("压缩率、中心最大预条件位移",contraction.str(12),move.str(12))
    print("电弱标尺超荷及色支、普朗克标尺色支",
          intercept[0].str(12),intercept[2].str(12),planck.str(12))


if __name__ == "__main__":
    verify()
