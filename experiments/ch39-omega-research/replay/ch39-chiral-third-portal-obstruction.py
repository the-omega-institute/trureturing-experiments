"""检验原窗口之后首枚三代候选，并在所有非负谱多项式上比较质量阈值。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
common = run_path(str(ROOT / "ch39-chiral-common-portal-scale.py"))
fixed, chiral, two = common["fixed"], common["chiral"], common["two"]


def layer_table(read, model, height):
    rows = read(height)
    return {q: chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}


def verify():
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    def passes(layers):
        return all(layers[q][2]<9<layers[q][3] for q in model["charges"])
    for n,zero in enumerate(acb.zeta_zeros(12093,115),12093):
        layers = layer_table(read,model,zero.imag)
        if n==12207:
            assert passes(layers)
        else:
            assert any(layers[q][2]>9 or layers[q][3]<9 for q in model["charges"])

    scale = arb("3.725e14")
    unit = arb("2.435e18")/724**2
    targets, spectra, low_moments = [], [], []
    for n in (11661,11934,12207):
        levels = layer_table(read,model,acb.zeta_zero(n).imag)
        assert passes(levels)
        _,_,alpha = fixed["crossing_enclosure"](model,initial,levels,False,724,"29.040",".1184")
        assert arb(".117")<alpha<arb(".119")
        table = [((unit/model["mz"]).log()+l,w,B,f,s) for l,w,B,f,s
                 in two["thresholds"](model,levels,False)]
        target,moments = common["portal_at"](model,initial,(levels,1/alpha,table),scale,steps=65536)
        targets.append(target)
        low_moments.append(moments)
        active = [unit*h for h in levels[3][3:] if unit*h<scale]
        assert len(active)==14 and all(unit*h<scale or unit*h>scale for h in levels[3][3:])
        spectra.append([((mass/(80*unit))**2,mass**2*(1+2*(scale/mass).log())/(16*arb.pi()**2))
                        for mass in active])
    a,b,c = targets
    x,y,z = low_moments
    fit = arb_mat([[x[0],x[4]],[y[0],y[4]]]).solve(arb_mat([[a],[b]])).entries()
    # 消去共同常数项，避免把同一系数的不确定度重复计入第三枚。
    r = z[0]/x[0]
    residual = c-r*a-fit[1]*(z[4]-r*x[4])
    assert -arb("8.4e21")<residual<-arb("7.9e21")

    first,last = spectra[0],spectra[2]
    def moment(s,k):
        return sum(w*u**k for u,w in s)
    best = moment(last,5)/moment(first,5)
    for k in range(512):
        if k!=5:
            assert moment(last,k)/moment(first,k)>best
    top1 = first[-1][0]
    top3,weight = last[-1]
    assert all(u<top1 for u,_ in first[:-1]) and top3>top1
    tail = weight*(top3/top1)**512/sum(w*(u/top1)**512 for u,w in first)
    assert tail>best
    # 最小单项比即任意非负多项式的锐下界，由第五次单项取到。
    margin = best-c/a
    assert margin>arb("1e-5")
    error = (best*a-c)/(1+best)
    print("四次两项门户在第三枚的质量余额/GeV²",residual.str(8))
    print("全部非负谱门户的锐矩比、目标比、二谱最大残差下界/GeV²",
          best.str(12),(c/a).str(12),error.lower().str(12))

    # 允许有符号系数时，直接检查所得多项式对每个实参数的值。
    matrix = arb_mat([row[:3] for row in low_moments])
    constant,linear,quadratic = matrix.solve(arb_mat([[v] for v in targets])).entries()
    minimum = constant-linear**2/(4*quadratic)
    assert quadratic>0 and minimum>arb(".006") and linear<0
    assert max(constant.upper(),(constant+linear+quadratic).upper())<arb(".04")
    print("二次门户的系数与全球最小值",*[v.str(8) for v in (constant,linear,quadratic,minimum)])


if __name__=="__main__":
    verify()
