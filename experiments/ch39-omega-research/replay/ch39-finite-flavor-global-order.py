"""截断六面坐标约束下，全局真空所需最大绕数的最小位数。"""
import json
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, acb, arb_mat, ctx


def competitor_rule(basis, costs, center):
    # 由一组活动圈数解中心处的支持平面。裁去负权后重新归一，并用全部 35 个圈验证。
    active, points = basis
    matrix, rhs = [], []
    for q in active:
        en,gx,gy,*_ = center[q-1]
        matrix.append([*[costs[i][q-1] for i in points], gx, gy, arb(1)])
        rhs.append([en])
    matrix.append([*[arb(1) for _ in points],arb(0),arb(0),arb(0)])
    rhs.append([arb(1)])
    try:
        solution=arb_mat(matrix).solve(arb_mat(rhs)).entries()
    except ZeroDivisionError:
        return None
    weights=[max(arb(0),v.mid()) for v in solution[:-3]]
    mass=sum(weights)
    if not mass>0:
        return None
    return [w/mass for w in weights],solution[-3].mid(),solution[-2].mid()


def covers_box(basis, costs, center, rt, rd):
    result=competitor_rule(basis,costs,center)
    if result is None:
        return False
    weights,u,v=result
    for q,(en,gx,gy,hxx,hxy,hyy,ca,cb,cd) in enumerate(center,1):
        ref=sum(w*costs[i][q-1] for w,i in zip(weights,basis[1],strict=True))
        diff=en-ref-u*gx-v*gy
        if not diff>0:
            return False
        # F=E−uGx−vGy；以 (t,d) 参数化 (x,y)=(t,t−d)。
        dx=gx-u*hxx-v*hxy
        dy=gy-u*hxy-v*hyy
        first=abs(dx+dy)*rt+abs(dy)*rd
        second=q*q/2*(ca*(1+q*abs(u))*rt**2+cb*(1+q*abs(v))*(rt+rd)**2
                      +cd*(1+q*abs(u-v))*rd**2)
        if not diff-first-second>0:
            return False
    return True


def check_lower(root,lib,previous,tri,angles):
    delta=arb(0).union(arb("4.48e-9"))
    assert arb(23)/25*arb(2)**-20/14**2 < delta.upper()
    raw=[lib["coeff"](tri,q,delta) for q in range(1,36)]
    # 同时处理正负两支。共轭改变虚部符号，直接取其共同的包围。
    coeffs=[[acb(z.real,z.imag.union(-z.imag)) for z in row] for row in raw]
    # JSON 只供选择竞争相位及支撑。后面的严格逐格比较不依赖采样点是否满足三窗。
    data=json.loads(Path(__file__).with_name("ch39-eight-bit-competitors.json").read_text(),parse_float=str)
    points=[tuple(arb(t) for t in row) for row in data["points"]]
    bases=data["bases"]
    samples=[(float(x),float(y),i) for x,y,i in data["samples"]]
    costs=[[lib["energy"](q,c,x,y) for q,c in enumerate(coeffs,1)] for x,y in points]
    cut1,cut2=arb("18.94088212"),arb("18.94414547")
    lower,upper=run_path(str(root/"ch39-curvature-window-density-cap.py"))["eigen_band"]()
    assert cut1<lower<upper<cut2
    lo,hi=[previous["spectral_cut"](root,c) for c in (cut1,cut2)]
    phase,amplitude=angles
    dbox=((arb(253)/256/amplitude).acos()-abs(phase)).union((arb(252)/256/amplitude).acos()+abs(phase))
    pending=[(-arb.pi()/16,arb.pi()/16,dbox.lower(),dbox.upper(),-1)]
    used=set()
    leaves=[0,0]
    while pending:
        tl,tr,dl,dr,last=pending.pop()
        x,d=tl.union(tr),dl.union(dr)
        y=x-d
        if abs(y)>arb.pi()/16 or all(lo(s*x,s*y,s*d)[0]>=16 or hi(s*x,s*y,s*d)[1]<16 for s in (-1,1)):
            leaves[0]+=1
            continue
        rt,rd=(tr-tl)/2,(dr-dl)/2
        if max(rt,rd)<arb(".00065"):
            tm,dm=(tl+tr)/2,(dl+dr)/2
            ym=tm-dm
            center=[]
            for q,c in enumerate(coeffs,1):
                a,b,k=[z*acb(0,q*t).exp() for z,t in zip(c,(tm,ym,dm),strict=True)]
                en=sum(z.real for z in c)-(a+b+k).real
                center.append((en,q*(a+k).imag,q*(b-k).imag,q*q*(a+k).real,-q*q*k.real,
                               q*q*(b+k).real,*[abs(z).upper() for z in c]))
            xf,yf=float(tm),float(ym)
            near=sorted(samples,key=lambda p:(p[0]-xf)**2+(p[1]-yf)**2)
            choices=list(dict.fromkeys(([last] if last>=0 else [])+[p[2] for p in near]))
            proved=False
            for idx in choices[:(len(choices) if max(rt,rd)<arb(".00002") else 4)]:
                if covers_box(bases[idx],costs,center,rt,rd):
                    used.add(idx)
                    leaves[1]+=1
                    proved=True
                    break
            if proved:
                continue
            last=near[0][2]
        assert max(rt,rd)>arb("1e-8"), (x,d,leaves,used)
        if rt>rd:
            mid=(tl+tr)/2
            pending.extend(((tl,mid,dl,dr,last),(mid,tr,dl,dr,last)))
        else:
            mid=(dl+dr)/2
            pending.extend(((tl,tr,dl,mid,last),(tl,tr,mid,dr,last)))
    print("八位全码的谱和竞争排除叶格、所用基数",leaves,len(used),flush=True)


def seven_bit_example(lib, previous, shell, physical, tri):
    # 这份作用只需一枚真实高度来固定刚度；两只联络角在环面上自由变分。
    height=acb.zeta_zero(20000).imag
    offset=shell["soft"]["point"]["study"]["offset"]
    delta=arb(23)/25*arb(2)**-20/(height**2+offset)
    ticks=(415914529,2030364,1595884644,5604668226,2381502237)
    assert sum(ticks)==10**10
    modes=[(q,[arb(t)*c/10**10 for c in lib["coeff"](tri,q,delta)])
           for q,t in zip((1,30,32,33,34),ticks,strict=True)]
    x,y=arb(".0045"),arb("-.1376")
    box=arb("1e-6")
    def value(a,b):
        return sum(lib["energy"](q,c,a,b) for q,c in modes)
    trial=value(x,y)
    pending=[(1,0,0)]
    leaves=[0,0,0]
    while pending:
        n,i,j=pending.pop()
        a,b=[arb.pi()*(arb(2*k+1,1)/n-1) for k in (i,j)]
        if value(a,b)>trial:
            leaves[0]+=1
            continue
        if any(abs(a-sign*x)<box and abs(b-sign*y)<box for sign in (-1,1)):
            leaves[1]+=1
            continue
        gx,gy,*_=lib["derivatives"](modes,a,b)
        if not gx.contains(0) or not gy.contains(0):
            leaves[2]+=1
            continue
        assert n<2**33,(a,b)
        pending.extend((2*n,2*i+u,2*j+v) for u in (0,1) for v in (0,1))
    print("七位见证的全环面能量、方框、非驻点排除叶格",leaves,flush=True)
    flavor=physical["joint"]["flavor"]
    for sign in (-1,1):
        a,b=sign*x,sign*y
        _,_,xx,xy,yy=lib["derivatives"](modes,a+arb(0,box),b+arb(0,box))
        assert xx>60 and (xx-60)*(yy-60)>xy**2
        gx,gy,*_=lib["derivatives"](modes,a,b)
        radius=(gx.abs_upper()**2+gy.abs_upper()**2).sqrt()/60
        assert radius<arb("5.3e-11")
        clocks={724:acb(0,a+arb(0,radius)).exp(),804:acb(0,b+arb(0,radius)).exp()}
        c=flavor["singlet_at"](clocks)[1]
        assert (128*c).floor().unique_fmpz()==126 and (256*c).floor().unique_fmpz()==253
        print("七位全局极小点至中心的上界、六面坐标",radius.str(14),c.str(16),flush=True)
        previous["full_windows"](physical,a+arb(0,radius),b+arb(0,radius))


def verify(root,winding,ref):
    lib,previous=run_path(str(winding)),run_path(str(ref))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    physical=run_path(str(root/"ch39-first-three-mode-observables.py"))
    flavor=physical["joint"]["flavor"]
    ctx.prec=160
    tri=lib["terms"](shell)
    center=flavor["singlet"](11661)[1]
    assert (center*256).floor().unique_fmpz()==252 and (center*128).floor().unique_fmpz()==126
    a,b=[flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    check_lower(root,lib,previous,tri,(acb(a,b).arg(),abs(acb(a,b))))

    seven_bit_example(lib,previous,shell,physical,tri)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    winding=Path(sys.argv[2]) if len(sys.argv)>2 else root/"ch39-mincut-global-winding.py"
    ref=Path(sys.argv[3]) if len(sys.argv)>3 else root/"ch39-mincut-reference-ckm-action.py"
    verify(root,winding,ref)
