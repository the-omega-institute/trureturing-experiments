"""自由改变十九个荷六截止时，用四行连续对偶排除第一枚 Gram 探针对。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb, ctx


# E1 表在每个二进倍程取 32 个有理节点；细分编号与点积全部用整数。
SUBDIVISIONS = 32
FINE = 256
E_SCALE = 2**20
B_SCALE = 2**27
MIN_POWER = -200
MAX_POWER = 6


def tabulate_e1():
    samples = []
    for power in range(MIN_POWER,MAX_POWER):
        base = arb(2)**power
        for j in range(SUBDIVISIONS):
            y = base*(SUBDIVISIONS+j)/SUBDIVISIONS
            value = y.expint(1)*E_SCALE
            nearest = round(float(value.mid()))
            assert abs(value-nearest)<1
            samples.append(nearest)
    # E1(64) <= exp(-64)/64.
    samples.append(0)
    assert (-arb(64)).exp()/64 < arb(1)/E_SCALE
    assert all(a>=b>=0 for a,b in zip(samples,samples[1:]))
    assert FINE*max(samples) < 2**63
    return np.array(samples,dtype=np.int64)


def e1_fixed(y,table):
    # 复制到最近的左细格点，在粗格点间以整数作线性插值。
    assert np.all(np.isfinite(y)) and np.all(y>=2.**MIN_POWER)
    m,k = np.frexp(np.minimum(y,2.**MAX_POWER))
    j = (m*(2*SUBDIVISIONS*FINE)).astype(np.int64)-SUBDIVISIONS*FINE
    index = (k.astype(np.int64)-1-MIN_POWER)*SUBDIVISIONS+j//FINE
    frac = j%FINE
    right = np.minimum(index+1,len(table)-1)
    return ((FINE-frac)*table[index]+frac*table[right])//FINE


def near_float(v):
    mid = float(v.mid())
    assert mid>0 and abs(arb(mid)/v-1)<arb("1e-13")
    return mid


def gauge_combination(species,rule,p0,p1,dual):
    nodes,quad_error = rule
    p,m = dual
    charge = [-(p*(arb(qq)-arb(8)*cc/3)+m*(mm-cc))/(4*arb.pi())
              for _,_,(qq,mm,cc) in species]
    charge_norm = sum(map(abs,charge))
    coefficients,low,high = [],[],[]
    for (kind,z,_),c in zip(species,charge,strict=True):
        for h,omega in nodes:
            coefficient = c*omega[kind]
            integer = round(float((coefficient*B_SCALE).mid()))
            assert abs(coefficient*B_SCALE-integer)<1
            coefficients.append(integer)
            low.append(near_float(z+p0*h))
            high.append(near_float(z+p1*h))
    log = (p1/p0).log()
    table = tabulate_e1()
    # int64 点积任一部分和都在此包围内；单个 E1 值也小于表的首值。
    assert sum(map(abs,coefficients))*int(table[0]) < 2**63
    size = arb(sum(map(abs,coefficients)))/B_SCALE
    coeff = np.array(coefficients,dtype=np.int64)
    low,high = np.array(low),np.array(high)
    # y^2 E1''(y) <= 1，|d E1/d log y| <= 1。前者给粗格线性插值余项，
    # 后者控制细格截断和两个经 Arb 校验的浮点输入及一次二进除法。
    e_error = 1/(8*arb(SUBDIVISIONS)**2)+1/arb(SUBDIVISIONS*FINE)+arb(2)/E_SCALE+arb("1e-11")
    error = charge_norm*quad_error*log + len(coefficients)*log/B_SCALE + 2*size*e_error
    def at(x):
        middle = near_float(x)
        e = e1_fixed(low/middle,table)-e1_fixed(high/middle,table)
        return arb(int(coeff@e))/(E_SCALE*B_SCALE)+arb(0,error.upper())
    return at,charge_norm,error


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number = 2_306_888_100_089
    _,source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec = 192
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    species,_ = neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    p0,p1 = data["p0"],data["p1"]
    lam0,lamd,lamp,lamm = map(arb,("5.81545","175.9075","1.223486","-7.127971"))
    a1,a2,a3 = data["initial"]
    target = lam0-lamp*(5*a1/3+a2-8*a3/3)-lamm*(a2-a3)
    gauge_at,charge_norm,error = gauge_combination(species,gauge["angular_rule"](11,25),p0,p1,(lamp,lamm))
    low_weak = lib("native-rh-ray-pvm")["weak_series"](56)
    x = arb("1.2768")*p0
    assert x<min(s.lower() for s in source)
    step = arb(17)/16
    # 对数坐标下 E1(a/x)-E1(b/x) 的二阶导数为两个 v exp(-v) 的差。
    # 插值误差在整个格子内统一；D 则对子段直接取带余项的弱差包围。
    curvature_error = charge_norm*step.log()**2/(8*arb(1).exp())
    high_weak = neutral["weak_column"]
    left = gauge_at(x)
    supremum = arb(0)
    cells = 0
    while x < arb(".58"):
        y = x*step
        right = gauge_at(y)
        gauge_band = left.union(right)+arb(0,curvature_error.upper())
        # 低端的弱差陡峭，配合较细的同一对数段线性插值。
        pieces = 32 if x<256*p0 else 1
        for j in range(pieces):
            def gauge_end(k):
                return ((pieces-k)*left+k*right)/pieces
            def part(k):
                return x*step**(arb(k)/pieces)
            sub = part(j).union(part(j+1))
            g = gauge_band if pieces==1 else gauge_end(j).union(gauge_end(j+1))+arb(0,curvature_error.upper())
            if x<arb("1e5")*p0:
                d = low_weak(sub/p0)
            else:
                d = high_weak(sub,p0,p1)
            supremum = max(supremum,abs(lam0+g+lamd*d).upper())
        x,left = y,right
        cells += 1
    assert x>max(s.upper() for s in source)
    lower = (target/supremum).lower()
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)[1]
    cap = sum(c*c for c in gram).sqrt()/abs(gram[1])
    print("截止格数、规范角积分端点误差、格间插值误差",cells,error.str(12),curvature_error.str(12))
    print("对偶目标、连续一致上界、总变差下界、第一探针上界",
          *[v.str(14) for v in (target,supremum,lower,cap)],flush=True)
    assert cells==1288 and supremum<arb("1.044454")
    assert lower>arb("26.27861") and cap<arb("24.97937")
    return lower


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
