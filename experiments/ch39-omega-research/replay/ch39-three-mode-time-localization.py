"""算术时间前缀对三模支撑概率与共同逆高度尺度的全谱限制。"""
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx


root = Path(__file__).parent
clock = run_path(str(root/"ch39-arithmetic-height-bandwidth.py"))
mask = run_path(str(root/"ch39-three-mode-arithmetic-mass.py"))


def partial_fractions(degree):
    rates = [j*j for j in range(1,degree+1)]
    return [(rate,prod(1-arb(other)/rate for other in rates if other != rate))
            for rate in rates]


def residues_of(numerator, poles):
    return [numerator(-arb(200)/rate)/denominator for rate,denominator in poles]


def window_polynomial(moments, x):
    centered = arb_mat([[moments[i+j+2]-2*x*moments[i+j+1]+x*x*moments[i+j]
                         for j in range(20)] for i in range(20)])
    basis = arb_mat([[x**i] for i in range(20)])
    coefficients = centered.solve(basis)
    coefficients = coefficients/(basis.transpose()*coefficients)[0,0]
    return arb_poly([coefficients[i,0].mid() for i in range(20)])


def weight_cap(times, cap, degree, targets, error=0):
    moments = clock["cleared_moments"](times,degree)
    order = degree//2
    gram = arb_mat([[moments[i+j] for j in range(order+1)] for i in range(order+1)])
    basis = arb_mat([[x**i for x in targets] for i in range(order+1)])
    q = arb_poly([1])
    for j in range(1,degree+1):
        q *= arb_poly([1,arb(j*j)/200])
    values = arb_mat([[q(x).sqrt()] for x in targets])
    inverse_basis = gram.solve(basis)
    coefficients = inverse_basis*(basis.transpose()*inverse_basis).solve(values)
    polynomial = arb_poly([coefficients[j,0].mid() for j in range(order+1)])
    power = arb_poly([polynomial[j]*cap**j for j in range(order+1)])
    normalization = 1+arb("1e-9")
    difference = normalization*power**2-arb_poly([q[j]*cap**j for j in range(degree+1)])
    for i in range(degree+1):
        bernstein = sum(difference[j]*comb(i,j)/comb(degree,j) for j in range(i+1))
        assert bernstein > 0
    numerator = normalization*polynomial**2
    nominal = sum(numerator[i]*moments[i] for i in range(numerator.degree()+1))
    if error:
        return nominal+error*sum(abs(c) for c in residues_of(numerator,partial_fractions(degree)))
    return nominal


def required_windows(times):
    moments = clock["cleared_moments"](times,41)
    poles = partial_fractions(41)
    result = []
    for n, tol in ((1,"2e-16"), (2,"3.2e-12")):
        height = acb.zeta_zero(n).imag
        x = 200/height**2
        polynomial = window_polynomial(moments,x)
        numerator = polynomial**2*(arb_poly([(x*arb(tol))**2])-arb_poly([-x,1])**2)
        # 分子至多四十次。用四十一个不同的极点展开后，各系数直接乘原始读数。
        residues = residues_of(numerator,poles)
        signal = sum(c*t for c,t in zip(residues,times[1:],strict=True))
        margin = signal/sum(abs(c) for c in residues)
        assert margin > arb(2)**-176
        print("窗口",n,"允许的逐项绝对误差裕量",margin.str(12))
        result.append((height,arb(1) + arb(0,tol)))
    return result


def excluded_phase_box():
    # 极限荷一磁壳取两相位 (i,-1)，三道界的两侧惯性不变。
    window = mask["window"]
    levels = ((15-arb(189).sqrt())/2,(15+arb(189).sqrt())/2,(15+arb(261).sqrt())/2)
    phases = {724:acb(0,1),804:acb(-1)}
    for sign in (-1,1):
        counts = [1+window["supplier"](1,t+sign*arb(".01"))(phases) for t in levels]
        assert counts == [0,7,12] and counts[0]+counts[2]-counts[1] == 5
    assert arb(8).sqrt()/1000 + arb("3e-8") < arb(".01")


def verify(hits):
    ctx.prec,ctx.cap = 1536,3
    excluded_phase_box()
    times = clock["times"]()
    cap = 200/acb.zeta_zero(200).imag**2
    for degree, targets, ceiling in ((2,[cap],".58565"),(20,[cap],".16508"),(40,[0,cap],".12182")):
        weight = weight_cap(times,cap,degree,targets)
        assert weight < arb(ceiling)
        print(degree,"个原单位时间读数容许的三模支撑概率上界",weight.str(11))
    (h1,window1),(h2,window2) = required_windows(times)
    ctx.prec = 128
    for candidate in hits:
        second = candidate * (h2/h1) * (window1/window2).sqrt()
        assert second.zeta_nzeros().unique_fmpz() is not None
    floor = acb.zeta_zero(8193).imag / h1 * window1.lower().sqrt()
    assert floor > arb("588.4188")
    print("同时满足前 41 项所需的共同逆高度因子下界",floor.str(13))


if __name__ == "__main__":
    verify(mask["verify"]())
