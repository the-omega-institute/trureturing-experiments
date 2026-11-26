"""四个算术 Gram 源嵌入极点与三射线，由原磁拉普拉斯谱投影给截止。"""
from pathlib import Path
from runpy import run_path
import sys
from math import comb, factorial

from flint import acb, arb, ctx


def gram_cholesky(root):
    record = run_path(str(root/"ch39-finite-rh-pivot-mass.py"))
    table = [[arb(0) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(i+1):
            lo,hi = record["read"](i,j)
            g = (arb(lo.numerator)/lo.denominator).union(arb(hi.numerator)/hi.denominator)
            value = g-sum(table[i][k]*table[j][k] for k in range(j))
            if i==j:
                assert value>0
                table[i][j] = value.sqrt()
            else:
                table[i][j] = value/table[j][j]
    return table


def readout(eigenvectors, order, sites, coefficients, pivot):
    weights = []
    for j in order:
        norm = sum(abs(eigenvectors[k,j])**2 for k in range(19))
        incoming = sum(coefficients[k]*eigenvectors[site,j] for k,site in enumerate(sites))
        weights.append(incoming*eigenvectors[sites[pivot],j].conjugate()/(norm*coefficients[pivot]))
    assert sum(weights).contains(1)
    return weights


def weak_series(order=32):
    # E1 的对数部分沿用已证的参考值；幂级数各项以 beta 积分精确求角矩。
    r = (arb("246.22")/arb("91.1876"))**2/2
    w = arb(".63811")**2*r/2
    z = w+3*arb(".467453")**2*r/10
    t = arb(".8616")**2*r
    assert max(w,z,t)+arb(1)/4<4

    def moment(n, mass, linear=False, h_power=0, average_v=False, omit_j0=False):
        total = arb(0)
        for j in range(int(omit_j0),n+1):
            a = j+h_power+(n-j if linear else 0)
            b = j+h_power
            coefficient = arb(comb(n,j)*factorial(a)*factorial(b))/factorial(a+b+1)
            if average_v:
                coefficient /= j+1
            total += coefficient*mass**(n-j)
        return total

    coefficients = []
    for n in range(1,order+1):
        charged = w*moment(n,w,omit_j0=True)
        neutral = z*moment(n,z,linear=True,omit_j0=True)
        axial = t*moment(n,t,omit_j0=True)
        vector = 6*(moment(n,t,h_power=1)-moment(n,arb(0),h_power=1))
        scalar = moment(n,z,linear=True,h_power=1,average_v=True)-moment(n,w,h_power=1,average_v=True)
        term = (-(2*charged+neutral)/4+scalar/8+3*axial/8-vector/24)/arb.pi()
        coefficients.append((-1)**n*term/(n*factorial(n)))
    remainder = 10*arb(4)**(order+1)/((order+1)*factorial(order+1))
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p_ratio = (arb("91.1876")/arb("1.176397e14"))**2

    def at(c):
        # x=c p0, c>=1. The finite-cutoff high endpoint costs at most 2x/p1.
        assert c>arb(".99")
        inverse = 1/c
        correction = arb(0)
        for a in reversed(coefficients):
            correction = (correction+a)*inverse
        return reference+correction+arb(0,(remainder+2*c*p_ratio).upper())
    assert arb("-.048342336655422")<at(arb(1))<arb("-.048342336655420")
    return at


def shell_weights(root, number=2_306_888_100_085):
    """逐荷产生 (q, 有序19谱值, 三组复谱权, 镜向三组复谱权)。"""
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    gram = gram_cholesky(root)
    pole,*rays = [shell["soft"]["pole"],*shell["boundary_sites"]()]
    assert [pole,*rays]==[14,8,7,6]
    height = acb.zeta_zero(number).imag
    for q in (1,-4,2,-3,6):
        matrix = shell["laplacian"](q,height)
        vals,vecs = matrix.eig(right=True,algorithm="rump")
        assert all(v.imag.contains(0) for v in vals)
        order = sorted(range(19),key=lambda k:float(vals[k].real))
        levels = [vals[k].real for k in order]
        assert levels[0]>0 and all(a<b for a,b in zip(levels,levels[1:]))
        direct = [readout(vecs,order,[pole,*rays],gram[j],j) for j in (1,2,3)]
        mirror = [readout(vecs,order,[pole,rays[0],rays[2],rays[1]],gram[j],j) for j in (1,2,3)]
        yield q,levels,direct,mirror


def bracket_balance(weak_at, ratios, weights):
    def response(scale):
        return sum(c*weak_at(scale*r) for c,r in zip(weights,ratios,strict=True))
    low,high = arb(1),arb(2)
    assert response(low)>0
    while response(high)>0:
        low,high = high,2*high
        assert high<4096
    assert response(high)<0
    for _ in range(24):
        center=(low+high)/2
        value=response(center)
        if value>0:
            low=center
        else:
            assert value<0
            high=center
    mass=arb("91.1876")*low.union(high).sqrt()
    return mass


def verify(details):
    root = Path(__file__).resolve().parent
    ctx.prec = 256
    weak_at = weak_series()
    weak_cost = arb("1.86547379546839").union(arb("1.86547379546841"))
    for q,levels,direct,mirrored in shell_weights(root):
        lower_endpoint = [weak_at(arb(1) if k==0 else v/levels[0]) for k,v in enumerate(levels)]
        for pivot in (1,2,3):
            weights,mirror = direct[pivot-1],mirrored[pivot-1]
            for c,m in zip(weights,mirror,strict=True):
                assert c.real.overlaps(m.real)
            real = [v.real for v in weights]
            # 零权允许由循环特征标强制；绝对值写法也覆盖区间含零。
            negative = (sum(abs(v) for v in real)-1)/2
            tails = [sum(real[j:]) for j in range(1,19)]
            depth = -min(v.lower() for v in tails)
            upper_min = min(v.upper() for v in tails)
            depth = (-upper_min).union(depth)
            assert depth>0
            if pivot==1 and q!=6:
                assert depth<weak_cost
            else:
                assert depth>weak_cost
            if pivot==1 and q==-4:
                assert negative>weak_cost
            if pivot==3:
                assert negative>arb("43.42808")
            else:
                assert negative<arb("43.42807")
            # 附加的共同质量标尺 x_j=lambda_j/lambda_max，采用逐节点大截止余量。
            total_variation = sum(abs(v) for v in real)
            tiny = arb("2.398e-21")+arb("7.468e-36")*levels[-1]/levels[0]
            error = (total_variation*tiny).upper()
            assert error<arb("2e-18")
            weak_at_floor = sum(c*d for c,d in zip(real,lower_endpoint,strict=True))
            if pivot==1:
                assert weak_at_floor<arb("-.03")
            else:
                assert weak_at_floor>arb(".25")
            minimum_at = min(range(1,19),key=lambda j:float(tails[j-1]))+1
            print(q,pivot,"负量、最大轮廓负深、最深后缀起点",negative.str(9),
                  depth.str(9),minimum_at,"最低截止取 Z 质量时的弱差",weak_at_floor.str(10),flush=True)
            if q==1 and pivot==2:
                assert arb("-6.094")<weights[13].imag<arb("-6.091")
                assert arb("6.091")<mirror[13].imag<arb("6.094")
            if pivot>1:
                ratios=[arb(1),*(e/levels[0] for e in levels[1:])]
                mass=bracket_balance(weak_at,ratios,real)
                assert mass>arb("91.1876")
                print("同一比例标尺的相消点最低截止 GeV",q,pivot,mass.str(12))
            else:
                # 令 u=p0/x_min，先对 [1/4096,1] 严格覆盖，再用大截止统一余量。
                ratios=[arb(1),*(e/levels[0] for e in levels[1:])]
                ends=[arb(1)/4096,(1+arb(1)/4096)/2,arb(1)]
                for a,b in zip(ends,ends[1:]):
                    upper=sum(c*weak_at((r/b).union(r/a)) for c,r in zip(real,ratios,strict=True))
                    assert upper<arb("-.003")
                p0=(arb("91.1876")/arb("1.22089e19"))**2
                error=total_variation*(arb("2.398e-21")+arb("7.468e-36")/(4096*p0))
                assert arb("-.07425")+error<arb("-.073")
            if details:
                for j,(e,c) in enumerate(zip(levels,weights,strict=True),1):
                    print("层、谱值、复交叉权",j,e.str(12),c.str(12))


if __name__=="__main__":
    verify("--levels" in sys.argv[1:])
