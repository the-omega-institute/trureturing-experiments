"""固定第五壳源见证，只调第16线的一个份额及共同标尺，精确恢复引力双矩。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb,arb_mat

ROOT = Path(__file__).resolve().parent
SOURCE = 2345698
HOLDOUT = tuple(range(2345700,2345723,2))
FORWARD = tuple(range(2345724,2345747,2))


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def setup():
    frozen = load("fifth-shell-frozen-witnesses")
    layout = frozen["data"]("ch39-fifth-shell-quadratic-layout.json")
    assert layout["heavy"]==SOURCE and layout["extra"][0]==[16,11000]
    shell = 2306888100089
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](ROOT,load("native-rh-ray-pvm")["shell_weights"],shell)
    labels = [arb(k)/11000 for row in layout["base"] for k in row]
    extra = [(None,i,arb(k)/11000) for i,k in layout["extra"]]
    shares,_ = load("one-background-specific-line")["solve_layout"](
        [labels],extra,set(),["11.5"],[shell],[SOURCE],False,reference_gev=10**12)
    kernel = load("continuous-split-roots")
    scale = arb(7)/1000
    def point(i,t):
        a,b,w = edges[i]
        return w*kernel["basis"](scale*((1-t)*a+t*b))
    base = [point(i,t) for i,t in enumerate(labels)]
    moments = sum(base[1:],base[0])
    changes = [point(i,t)-base[i] for who,i,t in extra]
    for share,delta in zip(shares,changes,strict=True):
        moments += share*delta
    return edges,kernel,moments,changes[0],shares[0]


def dilate(row,s,derivative=False):
    """在原基底中表示改变标尺比后的核，或其对标尺比的导数。"""
    if derivative:
        return load("fifth-shell-rescaled-gravity")["derivative"](row,s)
    entries = row.entries()
    result = [value*s**(2-i) for i,value in enumerate(entries[:-1])]
    result[2] += entries[-1]*s.log()
    return arb_mat([[*result,entries[-1]]])


def bounds(s,d,coeff,moments,change,tail,delta_tail):
    """两残差及其雅可比，s 为标尺比，d 为第16线额外份额的增量。"""
    m = moments+d*change
    f,jac = [],[]
    for n,(row,target) in enumerate(zip(coeff,[arb(1),arb("1e-122")],strict=True)):
        value = dilate(row,s)
        deriv = dilate(row,s,True)
        f.append((value*m)[0,0]-target+arb(0,tail[n].upper()))
        jac.append([(deriv*m)[0,0]+arb(0,(2*tail[n]).upper()),
                    (value*change)[0,0]+arb(0,delta_tail[n].upper())])
    return arb_mat([[z] for z in f]),arb_mat(jac)


def slopes(matrix):
    return [sum(abs(matrix[i,j]) for j in range(2)).upper() for i in range(2)]


def repair_box(args, share):
    # 每个源份额由原四式的区间解给出，固定盒中的 Newton 中心只用二进制中点。
    center = [arb(1),arb(0)]
    for _ in range(3):
        val,jac = bounds(*center,*args)
        move = jac.solve(val)
        center = [(x-v).mid() for x,v in zip(center,move.entries(),strict=True)]
    radius = arb("2e-9")
    half = arb(0,radius)
    boxes = [x+half for x in center]
    assert boxes[0]>arb(".99") and boxes[0]<arb("1.01") and 0<share+boxes[1]<1
    at,center_jac = bounds(*center,*args)
    _,jac = bounds(*boxes,*args)
    inv = arb_mat(2,2,[x.mid() for x in center_jac.inv().entries()])
    contraction = max(slopes(arb_mat([[1,0],[0,1]])-inv*jac))
    error = max(abs(v).upper() for v in (inv*at).entries())
    assert contraction<1 and error+contraction*radius<radius,(error,contraction)
    # 解坐标重新包紧以验极小漂移符号。
    root = [z+arb(0,(error/(1-contraction)).upper()) for z in center]
    return root,jac


def verify(targets=HOLDOUT):
    edges,kernel,moments,change,share = setup()
    _,lib,bgs = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,(SOURCE,*targets),tight_neutral=True)
    model,_,_,_,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    total = sum(abs(w) for a,b,w in edges)
    order = kernel["ORDER"]
    wc = abs(edges[16][2])
    source_rows = bgs[0]
    denom = 2*sum(a for g,a,z in source_rows)
    source_mass = sum(g*z for g,a,z in source_rows)
    slopes_lower,shifts_s,shifts_d,relative = [],[],[],[]
    for number,rows in zip(targets,bgs[1:],strict=True):
        assert max(z.upper() for g,a,z in rows)<min(arb(1),arb(".0069")*min(a.lower() for a,b,w in edges))
        series = kernel["series"](rows,fields)
        coeff = [arb_mat([[int(i==n) for i in range(5)]])*series for n in (0,1)]
        tails = [2*size/(4*arb.pi()*factorial(order+1)) for size in
                 [sum(abs(2*a) for g,a,z in rows),sum(abs(g) for g,a,z in rows)]]
        args = (coeff,moments,change,[t*total for t in tails],[2*t*wc for t in tails])
        (s,d),jac = repair_box(args,share)
        schur = jac[1,1]-jac[1,0]*jac[0,1]/jac[0,0]
        assert jac[0,0]>0 and not schur.contains(0)
        lead = (sum(g*z for g,a,z in rows)-source_mass)/denom
        # 仅取源端确定的导数与质量平方和，预测更换重谱后所需的份额方向。
        if number==targets[0]:
            src = kernel["series"](source_rows,fields)
            source_coeff = [arb_mat([[int(i==n) for i in range(5)]])*src for n in (0,1)]
            assert max(z.upper() for g,a,z in source_rows)<min(arb(1),arb(".0069")*min(a.lower() for a,b,w in edges))
            source_tails = [2*z/(4*arb.pi()*factorial(order+1)) for z in
                            [sum(abs(2*a) for g,a,x in source_rows),sum(abs(g) for g,a,x in source_rows)]]
            _,J = bounds(arb(1),arb(0),source_coeff,moments,change,
                         [x*total for x in source_tails],[2*x*wc for x in source_tails])
            source_slope = J[1,1]-J[1,0]*J[0,1]/J[0,0]
        predicted = -lead/source_slope
        assert predicted*d>0
        rel = abs((d-predicted)/predicted)
        shifts_s.append(abs(s-1).upper());shifts_d.append(abs(d).upper())
        slopes_lower.append(abs(schur).lower());relative.append(rel.upper())
        print(number,"标尺比改变量、份额改变量、质量和预报的相对误差",
              (s-1).str(9),d.str(9),rel.str(8),flush=True)
    assert min(slopes_lower)>arb(".01")
    assert max(shifts_s)<arb(".0001") and max(shifts_d)<arb(".0001") and max(relative)<arb(".01")
    if targets==FORWARD:
        assert max(shifts_s)<arb(".00002379") and max(shifts_d)<arb(".00004678")
        assert max(relative)<arb(".0000437")
    if targets==HOLDOUT:
        assert max(shifts_s)<arb(".0000189") and max(shifts_d)<arb(".0000373") and max(relative)<arb(".00013")
    print("最大标尺比改变量、最大份额改变量、真空对份额斜率绝对值下界、最大预报相对误差",
          max(shifts_s),max(shifts_d),min(slopes_lower),max(relative))


if __name__=="__main__":
    verify()
    verify(FORWARD)
    load("fifth-shell-rescaled-gravity")["verify"](forward=FORWARD)
