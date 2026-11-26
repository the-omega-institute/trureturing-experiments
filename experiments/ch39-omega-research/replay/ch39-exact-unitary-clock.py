"""在不使用 11431 的六段时间表中精确复制它的规范均值和末态幺正矩阵。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

free = run_path(str(Path(__file__).with_name("ch39-center-free-three-pulses.py")))
ordered, study, box, flavor = (free[key] for key in ("ordered", "study", "box", "flavor"))
labels = (11381, 11397, 11430, 11432, 11303, 11576)
seed = (".0025244007063487", ".0270792666395162", ".2280674739172651",
        ".6314691050814896", ".0546549231024736", ".0562048305529069")


def components(matrix):
    return [matrix[0,0].imag, matrix[0,1].real, matrix[0,1].imag]


def verify():
    ctx.prec = 240
    heights = [acb.zeta_zero(n).imag for n in (*labels,11431)]
    xs = [(heights[-1]**2+study["offset"])/(h**2+study["offset"])-1 for h in heights[:-1]]
    blocks = [study["native"]["plane_response"](h) for h in heights]
    traces = [(b-blocks[-1]).trace().real for b in blocks[:-1]]
    identity = flavor["identity"]
    common = (flavor["arcs"][32], flavor["arcs"][152], identity)
    centers = box["centers"]
    coeff_radius = arb(2)**-28
    cube = [v+arb(0,coeff_radius) for v in centers]

    def evaluate(times, coeff, *, coeff_derivatives=False):
        gates = []
        actions = []
        sources = []
        for block, duration in zip(blocks, [*times,arb(1)]):
            basis = (block, *common)
            action = sum((v*b for v,b in zip(coeff,basis)), acb_mat(2,2))
            actions.append(action)
            gates.append(ordered["evolution"](block,duration,coeff))
            if coeff_derivatives:
                modes = free["three"]["spectral_modes"](block,coeff)
                derivatives = []
                for direction in basis:
                    total = acb_mat(2,2)
                    for j,(a,proj_a) in enumerate(modes):
                        for k,(b,proj_b) in enumerate(modes):
                            diff = duration*(a-b)/2
                            sinc = arb(1) if j==k else diff.sin()/diff
                            total += (acb(0,duration)*acb(0,duration*(a+b)/2).exp()*sinc
                                      *proj_a*direction*proj_b)
                    derivatives.append(total)
                sources.append(derivatives)
        target = gates.pop()
        target_inverse = target.conjugate().transpose()
        before = []
        prod = identity
        for gate in gates:
            before.append(prod)
            prod = gate*prod
        after = [None]*len(gates)
        suffix = target_inverse
        for j in reversed(range(len(gates))):
            after[j] = suffix
            suffix = suffix*gates[j]
        relative = target_inverse*prod
        values = arb_mat([[v] for v in (sum(times)-1, sum(t*x for t,x in zip(times,xs)),
                                         sum(t*d for t,d in zip(times,traces)), *components(relative))])
        columns = []
        for j,gate in enumerate(gates):
            derivative = acb(0,1)*after[j]*actions[j]*gate*before[j]
            columns.append([1,xs[j],traces[j],*components(derivative)])
        jacobian = arb_mat(columns).transpose()
        if not coeff_derivatives:
            return values,jacobian,relative
        derivatives = []
        full_derivatives = []
        for i in range(4):
            total = sources[-1][i].conjugate().transpose()*prod
            for j in range(len(gates)):
                total += after[j]*sources[j][i]*before[j]
            derivatives.append([0,0,0,*components(total)])
            full_derivatives.append(total)
        return arb_mat(derivatives).transpose(), full_derivatives

    point = list(map(arb,seed))
    for _ in range(3):
        values,jac,_ = evaluate(point,centers)
        step = jac.solve(values)
        point = [(t-step[j,0]).mid() for j,t in enumerate(point)]
    residual,jac,_ = evaluate(point,centers)
    inverse = jac.mid().inv().mid()
    assert 0 not in inverse.det()
    # 在系数的整个闭盒上控制固定点残量，保留纯态和候选态对参数的共同依赖。
    sensitivity = inverse*evaluate(point,cube,coeff_derivatives=True)[0]
    forces = [abs((inverse*residual)[i,0]) +
              coeff_radius*sum(abs(sensitivity[i,j]) for j in range(4)) for i in range(6)]
    radii = [arb("1e-13")+2*v.upper() for v in forces]
    times = [v+arb(0,r) for v,r in zip(point,radii)]
    _,jac,relative = evaluate(times,cube)
    identity6 = arb_mat([[int(i==j) for j in range(6)] for i in range(6)])
    contraction = identity6-inverse*jac
    ratios = []
    for i in range(6):
        change = sum(abs(contraction[i,j])*radii[j] for j in range(6))
        ratios.append((forces[i]+change)/radii[i])
    print("Central times",*[v.str(21) for v in point])
    print("Tube radii",*[v.upper().str(9) for v in radii])
    print("Weighted contraction and displacement",*[v.upper().str(8) for v in ratios])
    assert all(v < 1 for v in ratios)
    assert all(t > 0 for t in times) and relative[0,0].real > 0
    # 保持六段相对比例，再将总作用时长减半，检验额外的末态读数。
    _,_,half_relative = evaluate([t/2 for t in times],cube)
    # evaluate 的参照时长为一，将相对矩阵换到半时长的纯态参照。
    center_half = ordered["evolution"](blocks[-1],arb(1)/2,cube)
    half_relative = center_half*half_relative
    error = (sum(abs(v)**2 for v in (half_relative-identity).entries())).sqrt()
    print("Half-time Frobenius separation",error.str(15))
    budget = run_path(str(Path(__file__).with_name("ch39-center-free-neighborhood-time.py")))
    _,scale,eps = budget["model"]()
    cost,*_ = budget["support"]["prev"]["setup"]()
    port = scale*sum(t*cost(x)[0] for t,x in zip(times,xs))/eps
    print("Port shift over original half-width",port.str(15))
    assert error > arb(".000491") and port > arb("538.28")
    # 用同一张整数表覆盖整个系数盒；近似见证独立于系数取值。
    length = 10**12
    ticks = [int((v*length+arb(1)/2).floor().unique_fmpz()) for v in point]
    ticks[3] += length-sum(ticks)
    fixed = [arb(n)/length for n in ticks]
    fixed_gauge = abs(sum(t*x for t,x in zip(fixed,xs)))
    _,_,mid_matrix = evaluate(fixed,centers)
    _,matrix_derivatives = evaluate(fixed,cube,coeff_derivatives=True)
    def frobenius(matrix):
        return sum(abs(v).upper()**2 for v in matrix.entries()).sqrt()
    fixed_error = frobenius(mid_matrix-identity) + coeff_radius*sum(map(frobenius,matrix_derivatives))
    print("Fixed integer clock",ticks,"gauge error upper",fixed_gauge.upper().str(10),
          "matrix error upper",fixed_error.upper().str(10))
    assert fixed_gauge < arb("1e-13") and fixed_error < arb("1e-9")
    _, derivatives = evaluate(times,cube,coeff_derivatives=True)
    slope = derivatives[1][0,1].imag
    assert slope < arb("-.000342")
    separation = -2*coeff_radius*slope
    assert separation > arb("2.54e-12")
    print("Two endpoint calibrations, imaginary off-diagonal difference above",separation.lower().str(12))
    return times


if __name__ == "__main__":
    verify()
