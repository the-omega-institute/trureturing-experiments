"""同一第四步精确满足完整夸克、第三枢轴四式与全程正动能。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


NUMBER = 2_306_888_100_090
QUATERNION = (-219150165, -404816044, -944357695, 999999448143)
DIRECTIONS = ((14,14,"i"), (8,7,"i"), (4,0,"r"), (7,5,"i"), (17,8,"r"), (14,4,"r"),
              (10,5,"r"), (9,6,"i"), (6,5,"r"), (18,14,"r"), (14,6,"i"), (11,5,"r"),
              (18,8,"i"), (9,7,"r"), (12,8,"i"), (8,6,"r"))


def setup(root, lib):
    gauge = lib("all-field-signed-crossing")
    data = gauge["joint_system"](root, "1.176397e14", number=NUMBER, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, NUMBER)
    weights = next(row[2][2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER) if row[0] == 6)
    weights = [v.real for v in weights]
    assert sum(weights).contains(1)
    rho = (arb("246.22") / arb("1.22089e19"))**2 / 2
    neutral = lib("full-neutral-free-cutoff")
    species, _ = neutral["response_rows"](data, masses, signs, rho, arb(".42995401872"))
    return data, species, weights


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    data, species, weights = setup(root, lib)
    gauge = lib("all-field-signed-crossing")
    means = lib("fourth-third-pivot-cyclic-five")["angular_means"]
    terms = [(kind,z,(mix-col,qq-8*col/3)) for kind,z,(qq,mix,col) in species]
    p0, p1 = data["p0"], data["p1"]
    log_ratio = (p1/p0).log()
    a1, a2, a3 = data["initial"]

    def endpoints(overlap, levels):
        low, h, k = levels
        anchors = [(low+(h-low)*abs(overlap[i,0])**2+(k-low)*abs(overlap[i,1])**2).mid()**2
                   for i in range(19)]
        ctx.prec = 192
        rule = gauge["angular_rule"](19, 52)
        center = [a2-a3, 5*a1/3+a2-8*a3/3]
        for x, c in zip(anchors, weights, strict=True):
            for kind,z,b in terms:
                value = c*gauge["angular_exposure"](kind,z,x,p0,p1,log_ratio,rule)/(4*arb.pi())
                for i in range(2):
                    center[i] -= b[i]*value

        def slopes(x):
            angle = [a-b for a,b in zip(means(p0/x), means(p1/x), strict=True)]
            result = [arb(0), arb(0)]
            for kind,z,b in terms:
                factor = -(-z/x).exp()*angle[kind]/(4*arb.pi())
                for i in range(2):
                    result[i] += factor*b[i]
            return result

        def evaluate(amps):
            values = center.copy()
            derivatives = [[], []]
            for a, ref, c in zip(amps, anchors, weights, strict=True):
                x = a*a
                delta = (x/ref).log()
                for i, grad in enumerate(slopes(ref.union(x))):
                    values[i] += c*grad*delta
                    derivatives[i].append(c*grad*2/a)
            return values, derivatives

        return evaluate

    sites = lib("fifth-flavor-cutoff-chart")["solve"](
        root, "ch39-fifth-flavor-four-positive", QUATERNION, DIRECTIONS, endpoints=endpoints)
    lower = min(x.lower() for x in sites)
    assert lower > arb("9.95e-11")
    running = lib("full-neutral-running")
    weak = running["weak_bound"](
        data, [p0, *sites], [arb(0), *weights], lib("neutral-axial-top"), lib("neutral-bosonic"))
    _, floors = lib("fourth-third-pivot-neutral-positive-path")["path"](
        data, species, weights, [s.mid() for s in sites], sites, weak, lib)
    weak_column = lib("full-neutral-free-cutoff")["weak_column"]
    weak_end = sum(c*weak_column(s,p0,p1)
                   for c,s in zip(weights,sites,strict=True))
    print("最小截止平方、全程弱差绝对上界、色与中性行列式下界、端点弱差",
          *[e.str(11) for e in (lower, weak.upper(), floors[2], floors[3], weak_end)], flush=True)
    assert arb("-.074256578") < weak_end < arb("-.074256576")
    assert floors[2] > arb("2.04041") and floors[3] > arb("1618.95")
    return sites


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
