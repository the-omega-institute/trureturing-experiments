"""第四壳十九个实际截止下，认证消去标量门户的中性组合在全能标半轴的符号。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb, ctx


INPUT_SCALE = 2**256
E1_SCALE = 2**20
WEIGHT_SCALE = 2**27


def bounds(x,scale):
    return int((x*scale).lower().floor().unique_fmpz()),int((x*scale).upper().ceil().unique_fmpz())


def ceiling(x):
    return -((-x.numerator)//x.denominator)


def angular_grid():
    # 在半角区间按二倍尺度分段，每一倍程再等分八段；其余半区由对称性给出。
    bits,parts = 20,8
    den = 2**bits*parts
    u = [j*2**i for i in range(bits-1) for j in range(parts,2*parts)]+[den//2]
    assert u[0]*2**bits==den and all(a<b for a,b in zip(u,u[1:]))
    h = np.array([i*(den-i) for i in u],dtype=object)
    def primitive(i):
        t = Fraction(i,den)
        return (3*t*t-2*t**3, 3*t-6*t*t+4*t**3,
                (5*t+6*t*t-4*t**3)/7, (9*t+6*t*t-4*t**3)/11)
    values = [primitive(i) for i in u]
    areas = [[2*(b[k]-a[k])*WEIGHT_SCALE for a,b in zip(values,values[1:])] for k in range(4)]
    wlo = np.array([[v.numerator//v.denominator for v in row] for row in areas],dtype=np.int64)
    whi = np.array([[ceiling(v) for v in row] for row in areas],dtype=np.int64)
    return h,den**2,wlo,whi,arb(2)**-bits


def evaluate(root,cutoffs,weights):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_085+4
    data=lib("all-field-signed-crossing")["joint_system"](root,"1.176397e14",order=4,last=8,number=number)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=192
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=lib("full-neutral-free-cutoff")["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    p0=data["p0"]
    a1,a2,a3=data["initial"]
    initial=a2-5*a1/3+2*a3/3
    assert len(cutoffs)==len(weights)==19 and all(0<x<1 for x in cutoffs)
    # 电磁、混合和色系数组合为 2 beta_QT3 - beta_QQ + 2 beta_c/3。
    # 标量行在模型中代数相消；所有含零却非精确零的区间也照常保留。
    rows=[(kind,mass,2*arb(qt)-qq+arb(2)*colour/3) for kind,mass,(qq,qt,colour) in species]
    assert all(b.contains(0) for kind,_,b in rows if kind==1)
    rows=[(kind,mass,b) for kind,mass,b in rows if not b.is_zero()]
    kinds=np.array([kind for kind,_,_ in rows])
    mass_lo,mass_hi = [np.array(col,dtype=object)[:,None] for col in zip(*(bounds(z,INPUT_SCALE) for _,z,_ in rows))]
    assert min(mass_lo[:,0])>=0
    h,hden,wlo,whi,delta=angular_grid()
    wlo,whi=wlo[kinds],whi[kinds]
    # 两端遗漏角区域用 E1(y)<=log(1+1/y)、h>=u/2 与 0<=omega<=3。
    tail=(6*delta*((2/(p0*delta)).log1p()+1)).upper()

    # E1 在正轴递减。在每倍程 32 个有理格点上预先向外取整函数值。
    count=32
    knots=[Fraction(2)**e*Fraction(j,count) for e in range(-180,8) for j in range(count,2*count)]
    e1=[bounds((arb(v.numerator)/v.denominator).expint(1),E1_SCALE) for v in knots]
    elo=np.array([a for a,_ in e1]+[0],dtype=np.int64)
    ehi=np.array([b for _,b in e1],dtype=np.int64)
    assert int(max(ehi))*int(max(whi.sum(axis=1)))<2**63
    nodes=[]
    for x,c in zip(cutoffs,weights,strict=True):
        xl,xr=bounds(x,INPUT_SCALE)
        floors=np.array([xl*v.numerator//v.denominator for v in knots],dtype=object)
        ceils=np.array([ceiling(xr*v) for v in knots],dtype=object)
        coeff=[c*b/(4*arb.pi()) for _,_,b in rows]
        nodes.append((floors,ceils,coeff))

    def endpoints(p):
        pl,pr=bounds(p,INPUT_SCALE)
        small=mass_lo+(pl*h//hden)[None,:]
        big=mass_hi+((pr*h+hden-1)//hden)[None,:]
        plus,minus=arb(0),arb(0)
        for floors,ceils,coeff in nodes:
            # big<=floor(x_lo*k) 给出上格点；ceil(x_hi*k)<=small 给出下格点。
            i=np.searchsorted(floors,big[:,1:],side="left")
            j=np.searchsorted(ceils,small[:,:-1],side="right")-1
            assert j.min()>=0
            lower=(elo[i]*wlo).sum(axis=1)
            upper=(ehi[j]*whi).sum(axis=1)
            for a,b,c in zip(lower,upper,coeff,strict=True):
                band=(arb(int(a))/(E1_SCALE*WEIGHT_SCALE)).union(
                    arb(int(b))/(E1_SCALE*WEIGHT_SCALE)+tail)
                # 每个角端点函数随 p 递减，正负系数分别累加。
                plus+=c.max(0)*band
                minus+=c.min(0)*band
        return plus,minus

    left=endpoints(p0)
    baseline=initial-left[0]-left[1]
    edges=[arb(2)**e for e in range(-113,11)]
    assert p0<edges[0]
    ceilings=[]
    for right in edges:
        pair=endpoints(right)
        ceilings.append((baseline+left[0]+pair[1]).upper())
        left=pair
    # 无穷端点的每个正核均为零；最后一段仍由同一单调性夹住。
    ceilings.append((baseline+left[0]).upper())
    worst=max(ceilings)
    print("初始组合、闭段含无穷尾数、全半轴上界、最差段序号、无穷尾上界",
          initial.str(15),len(ceilings),worst.str(13),ceilings.index(worst)+1,ceilings[-1].str(13),flush=True)
    assert worst < arb("-54.30904") and ceilings[-1] < arb("-1498.565")
    return worst


def verify(root):
    xs,ws=run_path(str(root/"ch39-fourth-full-ckm-gravity.py"))["verify"](root,return_cutoffs=True)
    return evaluate(root,xs,ws)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
