"""沿连续比较能标，从低能初值严格下包第四壳两项绝对规范残差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def angular_bands(parts):
    def primitive(t):
        return (3*t*t-2*t**3, 3*t-6*t*t+4*t**3,
                (5*t+6*t*t-4*t**3)/7, (9*t+6*t*t-4*t**3)/11)

    for i in range(parts):
        left, right = arb(i)/(2*parts), arb(i+1)/(2*parts)
        areas = [2*(b-a) for a,b in zip(primitive(left),primitive(right),strict=True)]
        assert all(v > 0 for v in areas)
        yield left*(1-left), right*(1-right), areas


def evaluate(root, cutoffs, weights):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    number = 2_306_888_100_089
    data = lib("all-field-signed-crossing")["joint_system"](
        root, "1.176397e14", order=4, last=8, number=number)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, number)
    ctx.prec = 192
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    species, _ = lib("full-neutral-free-cutoff")["response_rows"](
        data, masses, signs, rho, arb(".42995401872"))
    assert len(species)==119 and len(cutoffs)==len(weights)==19
    assert all(0 < x < 1 for x in cutoffs) and sum(weights).contains(1)

    # 对每行预先算出保留符号的 c_i exp(-z/x_i)。
    rows = []
    beta_norm = [arb(0), arb(0)]
    for kind, mass, (qq, mixed, colour) in species:
        assert mass.is_zero() or mass > 0
        beta = [qq-arb(8)*colour/3, mixed-colour]
        beta_norm = [a+abs(b) for a,b in zip(beta_norm,beta,strict=True)]
        coeff = [c*(-mass/x).exp() for x,c in zip(cutoffs,weights,strict=True)]
        rows.append((kind, mass, beta, coeff))

    a1,a2,a3 = data["initial"]
    current = [5*a1/3+a2-8*a3/3, a2-a3]
    floors = [v.lower() for v in current]
    bands = list(angular_bands(24))
    p = data["p0"]
    step = arb(2).log()
    norm = 4*arb.pi()
    while p < 2**16:
        # 每个矩形 [p,2p] x [t_0,t_1] 上 y=p' t(1-t) 的两个端点。
        integral = [arb(0), arb(0)]
        for hlo,hhi,areas in bands:
            lo,hi = p*hlo,2*p*hhi
            exponentials = [(-hi/x).exp().union((-lo/x).exp()) for x in cutoffs]
            for kind,mass,beta,coeff in rows:
                factor = arb(1) if mass.is_zero() else (lo/(lo+mass)).union(hi/(hi+mass))
                response = areas[kind]*factor*sum(a*b for a,b in zip(coeff,exponentials,strict=True))
                for j in range(2):
                    integral[j] += beta[j]*response
        for j in range(2):
            delta = (-step*integral[j]/norm).lower()
            current[j] = (current[j]+delta).lower()
            floors[j] = min(floors[j], current[j])
        p *= 2

    # h(t)>=t/2 到半角。两端积分给出 |dR/dlog p|<=12 x_max TV ||beta||_1/(4 pi p)。
    variation = sum(abs(c) for c in weights)
    tail = [12*max(x.upper() for x in cutoffs)*variation*b/(norm*p) for b in beta_norm]
    floors = [min(f,(v-t).lower()) for f,v,t in zip(floors,current,tail,strict=True)]
    print("两项规范残差的全程下界", *[v.str(14) for v in floors], flush=True)
    print("末端下界与剩余总变化上界", *[v.str(12) for v in current+tail], flush=True)
    assert floors[0] > arb("103.70") and floors[1] > arb("19.75")
    return floors


def verify(root):
    full = run_path(str(root/"ch39-fourth-full-ckm-gravity.py"))
    return evaluate(root, *full["verify"](root, return_cutoffs=True))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
