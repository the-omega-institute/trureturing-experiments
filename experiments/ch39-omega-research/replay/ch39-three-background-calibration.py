"""三组相邻背景的共同交汇、全程正性及单条独立线的分类。"""
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
SHELL = (2306888100085, 2306888100086, 2306888100087)
HEAVY = (2345666, 2345668, 2345670)
# 分母均为 10^11。各背景依次取 K−1、W−目标、A₁−A₃、A₂−A₃、A₃(M_Z)。
LINE_323 = ((-100000000000,-100000000000,-11908287254,34545312970,0),
            (-35137932052,-100000000000,896711509,-2586511152,0),
            (-99796078122,-100000000000,34491393033,-100000000000,0))
LINE_0 = ((-9236592476,-100000000000,5954064,-15763022,0),
          (742566174,-100000000000,-26802,68924,0),
          (1952164595,40854269349,-516375759,63152557,-644063941))
THIRD = ((2875417518,34873570380,-79557105379,8072027442,-100000000000),)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def primal():
    certificate = load("one-background-specific-line")["solve_layout"]
    examples = json.loads((ROOT/"ch39-three-background-layouts.json").read_text())
    for exponent, lower in [(15, ("2.39281",)*3), (16, ("-0.725557",)*3)]:
        layout = examples[str(exponent)]
        common = [arb(k)/11000 for row in layout["base"] for k in row]
        assert len(common)==361
        labels = [[arb(k)/11000, *common[1:]] for k in layout["independent_base"]]
        extra = [(n,i,arb(k)/11000) for n,i,k in layout["extra"]]
        print("共同交汇能标 GeV",10**exponent,flush=True)
        certificate(labels, extra, 0, lower, SHELL, HEAVY, False, reference_gev=10**exponent)
    # 取两份概率分布的逐线凸组合，读数对概率分布、交汇等式对 log E 均仿射。
    theta = arb("5.852").log()/arb(10).log()
    lower = (1-theta)*arb("2.39281")-theta*arb("0.725557")
    assert 0<theta<1 and lower>arb("0.000073")


def projection_terms(sectors,fields,beta,initial,inputs,kernel,directions,meeting,scale_millionths):
    """返回投影的基线、能标斜率及已扣热核尾界的逐线数据。"""
    scales = [7000]*len(sectors) if scale_millionths is None else scale_millionths
    order = kernel["ORDER"]
    tail = 2*arb(7)**(order+1)/factorial(order+1)
    assert arb(7)/(order+2) < arb(1)/2
    constant, slope = arb(0), arb(0)
    tables = []
    for (edges,series,bg),integers,micro in zip(sectors,directions,scales,strict=True):
        assert micro>0
        scale = arb(micro)/10**6
        k,h,u,v,s = [arb(z)/10**11 for z in integers]
        row = [k,h,u,v,s-u-v]
        coeff = arb_mat([row])*series
        tables.append([(scale*a,scale*(b-a),w*coeff) for a,b,w in edges])
        assert len(edges)==361 and all(0<a and (a<b or a==b) and b<=1 for a,b,w in edges)
        floor = min((scale*a).lower() for a,b,w in edges)
        assert max(z.upper() for g,j,z in bg)<min(arb(1),7*floor)
        size = sum(abs(2*k*j)+abs(h*g)+sum(abs(d*b) for d,b in zip(row[2:],bs,strict=True))
                   for (g,j,z),bs in zip(bg,fields,strict=True))/(4*arb.pi())
        constant -= tail*size*sum(abs(w) for a,b,w in edges)
        constant += -k-h*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])+s*initial[2]
        rate = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
        if meeting is None:
            constant -= (rate*(arb("1.22089e19")/inputs["mz"]).log()).max(0)
        else:
            constant -= rate*(arb(meeting)/inputs["mz"]).log()
        slope -= rate

    return constant,slope,list(zip(*tables,strict=True))


def projection(sectors,fields,beta,initial,inputs,kernel,directions,meeting,independent,
               *, scale_millionths=None):
    constant,slope,rows = projection_terms(sectors,fields,beta,initial,inputs,kernel,directions,
                                           meeting,scale_millionths)
    minimum = load("shared-continuous-position-gap")["projected_minimum"]
    for i,terms in enumerate(rows):
        groups = [[term] for term in terms] if i in independent else [terms]
        for group in groups:
            constant += minimum(kernel["basis"],kernel["S"],group,arb("1e-6"))
    return constant.lower(),slope


def dual():
    source = load("shared-twelve-node-gauge-gap")
    kernel = load("continuous-split-roots")
    _,lib,backgrounds = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,HEAVY,tight_neutral=True)
    inputs,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](inputs)]
    shell = load("native-rh-ray-pvm")["shell_weights"]
    beta = [b/(2*arb.pi()) for b in inputs["beta"]]
    sectors = [(source["paired_lines"](ROOT,shell,n),kernel["series"](bg,fields),bg)
               for n,bg in zip(SHELL,backgrounds,strict=True)]
    bound,_ = projection(sectors,fields,beta,initial,inputs,kernel,LINE_323,None,{323})
    norm = sum(abs(arb(x)/10**11) for row in LINE_323 for x in row[2:4])
    print("仅第 323 线独立，任意三处能标的六项规范差之最大绝对值下界",(bound/norm).str(12),flush=True)
    assert bound/norm > arb("2.80331")
    meeting = 6200000000000000
    bound,slope = projection(sectors,fields,beta,initial,inputs,kernel,LINE_0,meeting,{0})
    assert slope>0 and LINE_0[2][4]<0 and LINE_0[0][4]==LINE_0[1][4]==0
    upper = arb(meeting)*(-bound/slope).exp()
    assert upper < arb("6.025e15")
    print("仅第 0 线独立，三背景保持全程正的共同交汇能标之必要上界 GeV",upper.str(13),flush=True)
    bound,slope = projection(sectors[2:],fields,beta,initial,inputs,kernel,THIRD,10**16,set())
    assert bound > arb("0.33179") and slope > 0
    ceiling = arb(10**16)*(-bound/slope).exp()
    assert ceiling < arb("7.83e15")
    print("第三背景在 10^16 GeV 交汇时强逆耦合 A3(M_Z) 的上界",(-bound).str(13),flush=True)
    print("第三背景保持全程正的交汇能标之必要上界 GeV",ceiling.str(13),flush=True)


if __name__=="__main__":
    primal()
    dual()
