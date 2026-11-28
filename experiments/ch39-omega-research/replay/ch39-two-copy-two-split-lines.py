"""保持 357 条的原格点，以两条端点分拆和两条单点选位给出正的规范交汇。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


ROOT = Path(__file__).resolve().parent
SPLIT = (0, 16)
MOVING = (17, 338)
SEED = (".04570626066740993", ".33879910987562295",
        ".9686248814813744", ".14392524328127743")


def verify():
    def load(name):
        return run_path(str(ROOT / f"ch39-{name}.py"))

    edges, rows, heat, _ = load("two-copy-gram-cutoff-two-exceptions")["model_inputs"]()
    gravity = load("signed-cutoff-gravity")
    gauge = load("signed-gravity-unification")
    lib = load("chiral-minimal-neutral-gravity")
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    labels = "".join(load("two-copy-four-split-lines")["LABELS"])
    fields = gauge["gauge_weights"](model)
    scale = arb(7)/1000
    pi = arb.pi()
    logscale = (arb(10**16)/model["mz"]).log()
    beta = [b/(2*pi) for b in model["beta"]]
    baseline = arb_mat([[0], [0], *[[a-b*logscale] for a,b in zip(initial,beta,strict=True)]])

    def response(i, t):
        a, b, w = edges[i]
        x = scale*((1-t)*a+t*b)
        assert 0 < x <= 1
        column = gauge["threshold_column"](rows, gravity["moments"], heat, fields, x)
        return arb_mat([[w*v] for v in column[1:]])

    def derivative(i, t):
        a, b, w = edges[i]
        x = scale*((1-t)*a+t*b)
        exps = [(-z/x).exp() for _,_,z in rows]
        # d[x E₂(z/x)]/dx=e^(-z/x)，d[x² E₃(z/x)]/dx=x e^(-z/x)。
        k = sum(alpha*e for (_,alpha,_),e in zip(rows,exps,strict=True))/(2*pi)
        vacuum = -x*sum(g*e for (g,_,_),e in zip(rows,exps,strict=True))/(4*pi)
        shifts = [-sum(coeff[j]*e for coeff,e in zip(fields,exps[1:],strict=True))/(4*pi*x)
                  for j in range(3)]
        return arb_mat([[w*scale*(b-a)*v] for v in (k,vacuum,*shifts)])

    directions = []
    assert len(edges) == len(labels) == 361
    for i, digit in enumerate(labels):
        if i in SPLIT:
            right = response(i, arb(1))
            baseline += right
            directions.append(response(i, arb(0))-right)
        elif i not in MOVING:
            baseline += response(i, arb(int(digit,12))/11)

    def readings(params):
        f0, f16, t17, t338 = params
        return (baseline + directions[0]*f0 + directions[1]*f16
                + response(17,t17) + response(338,t338))

    project = arb_mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,-1,0],[0,0,1,0,-1]])
    target = arb_mat([[1],[arb("1e-122")],[0],[0]])

    def jacobian(params):
        columns = [*directions, *(derivative(i,t) for i,t in zip(MOVING,params[2:],strict=True))]
        return project*arb_mat([v.entries() for v in columns]).transpose()

    # 浮点初值只用于定位。预条件矩阵取逐项精确的二进制中点。
    center = [arb(v) for v in SEED]
    step = jacobian(center).solve(project*readings(center)-target)
    center = [(v-step[i,0]).mid() for i,v in enumerate(center)]
    J0 = jacobian(center)
    midpoint = arb_mat([[J0[i,j].mid() for j in range(4)] for i in range(4)])
    inverse = midpoint.inv()
    preconditioner = arb_mat([[inverse[i,j].mid() for j in range(4)] for i in range(4)])
    assert not preconditioner.det().contains(0)
    radius = arb(2)**-36
    box = [m + arb(0,radius.upper()) for m in center]
    assert all(0 < x < 1 for x in box)
    correction = preconditioner*(project*readings(center)-target)
    displacement = max(abs(z).upper() for z in correction.entries())
    identity = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    error = identity-preconditioner*jacobian(box)
    contraction = max(sum(abs(error[i,j]) for j in range(4)).upper() for i in range(4))
    assert contraction < 1 and displacement+radius*contraction < radius, (contraction, displacement, correction)
    # Banach 定理给盒内唯一精确解，且其距中心不超过 displacement/(1-contraction)。
    error_radius = displacement/(1-contraction)
    root = [m + arb(0,error_radius.upper()) for m in center]
    values = readings(root)
    common = values[4,0]
    assert beta[2] < beta[1] < 0 < beta[0] and logscale > 0
    minimum = common + beta[2]*logscale
    assert contraction < arb("1.062e-8") and error_radius < arb("5.584e-16")
    assert arb("37.5627922203") < common < arb("37.5627922205")
    assert arb("1.5461721635") < minimum < arb("1.5461721637")
    print("半宽、预条件残差、雅可比收缩界、根的半宽", radius.str(5),displacement.str(12),
          contraction.str(12),error_radius.str(12))
    print("两条线的左端点份额与另两条线的位置", *[z.str(17) for z in root])
    print("共同逆耦合和全程最小值",common.str(17),minimum.str(17))

    # 在盒内的同一个中心解处冻结所有选择，仅更换重场质量所用的零点。
    frozen = []
    for i, ((a, b, weight), digit) in enumerate(zip(edges, labels, strict=True)):
        if i in SPLIT:
            share = root[SPLIT.index(i)]
            choices = [(arb(0), share), (arb(1), 1-share)]
        elif i in MOVING:
            choices = [(root[2+MOVING.index(i)], arb(1))]
        else:
            choices = [(arb(int(digit, 12))/11, arb(1))]
        for t, share in choices:
            frozen.append((scale*((1-t)*a+t*b), weight*share))

    # 沿四个匹配方程隐式调节四个参数，求公共逆耦合对共同标尺的局部响应。
    scale_column = arb_mat([[0] for _ in range(5)])
    for x, weight in frozen:
        rates = [(-z/x).exp() for _, _, z in rows]
        dk = x*sum(alpha*e for (_,alpha,_),e in zip(rows,rates,strict=True))/(2*pi)
        dw = -x*x*sum(g*e for (g,_,_),e in zip(rows,rates,strict=True))/(4*pi)
        du = [-sum(b[j]*e for b,e in zip(fields,rates[1:],strict=True))/(4*pi) for j in range(3)]
        scale_column += arb_mat([[weight*f] for f in (dk,dw,*du)])
    response_matrix = arb_mat([v.entries() for v in
                               [*directions, *(derivative(i,t) for i,t in zip(MOVING,root[2:],strict=True))]]).transpose()
    system = project*response_matrix
    assert not system.det().contains(0)
    parameter_slopes = system.solve(-project*scale_column)
    common_slope = (scale_column+response_matrix*parameter_slopes)[4,0]
    assert arb("-.832249") < common_slope < arb("-.832248")
    print("沿精确交汇解支的公共逆耦合对 log c 的导数", common_slope.str(14))

    codes = (2_345_665, 2_345_667)
    _, _, neighbours = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT, codes, tight_neutral=True)
    for code, neighbour in zip(codes, neighbours, strict=True):
        result = arb_mat([[0], [0], *[[a-b*logscale] for a,b in zip(initial,beta,strict=True)]])
        for x, weight in frozen:
            column = gauge["threshold_column"](neighbour, gravity["moments"], heat, fields, x)
            result += arb_mat([[weight*v] for v in column[1:]])
        k, vacuum, u1, u2, u3 = result.entries()
        differences = (k-1, vacuum-arb("1e-122"), u1-u3, u2-u3)
        assert all(not v.contains(0) for v in differences)
        print("冻结两分拆解后的邻点编号，牛顿、真空、超荷减色、弱减色残差",
              code, *[v.str(13) for v in differences])


if __name__ == "__main__":
    verify()
