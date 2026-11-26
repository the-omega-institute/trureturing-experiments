"""径向单位 d/12 图上，弱 Casimir 动能项给出的一参数有限阈值交点。

固定 p=1、R=theta(804)/theta(12)、第 11431 个零点的连接和三代轻场。
仅在每个图上较高的十八个物质模取 Z=1+c C2(SU(2))；Q、L、H 的
C2=3/4，单态为零。图质量项维持原值。轻层的 SM 输入单独计一次。
这是一项另选的树级 Wilson 系数。用一圈自能单独确定 c 再插入一圈
极化，会得到更高阶贡献。

对剩余对角规范群，固定背景的图谱投影作用在站点空间，和表示生成元
对易。协变四维动能使同一个 Z 乘逆传播子的动量项、三点顶点以及
标量海鸥项。场正则化给出标量质量平方 lambda/Z；固定一阶费米质量块
且左右动能同乘 Z 时，费米质量平方为 lambda/Z²。也单列以 Weyl 指数
计数的二阶调节器 lambda/Z。物理矢量样重层须计入共轭伴侣。

若 Z 同乘整个二次算子，规范化后图质量也恢复。固定背景相对行列式中
log det Z 相消。协变测度可能改写局域匹配项（例如 F²）；给定 M_Z
初始耦合后，这种约定不会产生非局域阈值。此处 Z 是正实常量，未做
手征相位旋转。完整站点规范场/链场动力学需有自己的协变作用量。
"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq_mpoly_ctx

here = Path(__file__).parent
radial = run_path(str(here/"ch39-native-radial-weyl-exclusion.py"))
scalar = run_path(str(here/"ch39-native-shell-higgs-budget.py"))["smooth"]
sheet = run_path(str(here/"ch39-smooth-two-stiffness-crossing.py"))["sheet"]


def number(string):
    value = Fraction(str(string))
    return arb(value.numerator)/value.denominator


def ward_check():
    # One Cartesian component; sum the quadratic identity over spacetime.
    algebra = fmpq_mpoly_ctx.get(("Z", "k", "q", "M"))
    z, k, q, mass = algebra.gens()
    inverse = lambda momentum: z*momentum**2+mass**2
    assert q*z*(2*k+q) == inverse(k+q)-inverse(k)
    assert inverse(k+2*q)-2*inverse(k+q)+inverse(k) == 2*z*q**2
    assert z*(k+q)+mass-(z*k+mass) == z*q  # coefficient of gamma_mu
    assert z*k**2+z*mass**2 == z*(k**2+mass**2)
    assert z*k+z*mass == z*(k+mass)


def projection(values):
    return (115*values[0]-333*values[1]+218*values[2])/218


def certify(masses, radius, power, partners, c_bounds, scale_bounds, common_bounds):
    # partners=1 counts one Weyl index per high mode; partners=2 counts a Dirac pair.
    rows = [[number(v) for v in row] for row in
            (("1/15",1,"2/3"), ("8/15",0,"1/3"), ("2/15",0,"1/3"),
             ("1/5","1/3",0), ("2/5",0,0))]
    higgs = (number("1/10"),number("1/6"),0)
    gauge = (0,2,3)
    beta = (number("4.1"),number("-19/6"),-7)
    assert projection(beta).contains(0)
    em, weak, mz, pi = number("127.930"), number(".23129"), number("91.1876"), arb.pi()
    initial = (3*(1-weak)*em/5, weak*em, 1/number(".1180"))
    vector = sum(sheet(m,radius**2,True) for m in masses[0])

    def intercept(c):
        z = 1+3*c/4
        higgs_sum = sum(scalar(m/z,radius**2) for m in masses[4])
        fermions = []
        for j, row in enumerate(masses[1:]):
            divisor = z**power if j in (0,3) else 1
            fermions.append(sum(sheet(m/divisor,radius**2,False) for m in row))
        return [initial[i]+7*gauge[i]*vector/(4*pi)-higgs[i]*higgs_sum/(2*pi)
                -3*partners*sum(f*row[i] for f,row in zip(fermions,rows))/(2*pi)
                for i in range(3)]

    lo, hi = arb(0), arb(100)
    assert projection(intercept(lo)) < 0 < projection(intercept(hi))
    for _ in range(43):
        mid = (lo+hi)/2
        value = projection(intercept(mid))
        if value < 0:
            lo = mid
        else:
            assert value > 0
            hi = mid

    c = lo.union(hi)+arb(0,number("1e-8").upper())
    values = intercept((lo+hi)/2)
    scale_log = 2*pi*(values[0]-values[2])/(beta[0]-beta[2])
    logs = scale_log+arb(0,number("1e-5").upper())
    # Each reference subtraction is bounded by q0²/(8 min m²).
    # 150 bounds the absolute coefficients summed over the eighteen levels for each A_i.
    for i in range(3):
        assert 18*(7*gauge[i]/(4*pi)+higgs[i]/(2*pi)
                   +3*partners*sum(row[i] for row in rows)/(2*pi)) < 150
    min_mass2 = min(m.lower() for row in masses for m in row)/(1+3*c.upper()/4)**power
    error = 150*(radius*(-logs.lower()).exp())**2/(8*min_mass2)
    assert 0 < error < number("1e-16")
    uncertainty = arb(0,error.upper())
    # The light running cancels in P; only reference subtractions depend on L.
    assert projection([a+uncertainty for a in intercept(c.lower())]) < 0
    assert projection([a+uncertainty for a in intercept(c.upper())]) > 0
    # Positive matter indices make each intercept nonincreasing in c.
    estimates = [a.union(b)+uncertainty for a,b in zip(intercept(c.lower()),intercept(c.upper()))]
    def difference(log):
        return estimates[0]-estimates[2]-(beta[0]-beta[2])*log/(2*pi)
    assert difference(logs.lower()) > 0 > difference(logs.upper())
    # Opposite signs on all four faces give simultaneous P=0 and A1=A3.
    scale = mz*logs.exp()
    common = estimates[0]-beta[0]*logs/(2*pi)
    for value, bounds in ((c,c_bounds), (scale,scale_bounds), (common,common_bounds)):
        assert number(bounds[0]) < value < number(bounds[1])
    print("p=1, fermion mass-squared power, Weyl multiplicity", power, partners,
          "c, scale GeV, common inverse", c.str(13), scale.str(13), common.str(13))


def verify():
    ctx.prec = 256
    ward_check()
    blocks, unweighted = radial["distance_blocks"]()
    logs = radial["log_masses"](blocks, unweighted, arb(1))
    masses = [[(2*v).exp() for v in row] for row in logs]
    geo = radial["soft"]["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    certify(masses,radius,1,1,("65.00307","65.00309"),("4.6987e14","4.6989e14"),("38.4159","38.4161"))
    certify(masses,radius,2,1,("8.27139","8.27140"),("5.0233e14","5.0235e14"),("38.4123","38.4125"))
    certify(masses,radius,2,2,("4.56179","4.56181"),("5.0756e14","5.0758e14"),("38.1843","38.1845"))


if __name__ == "__main__":
    verify()
