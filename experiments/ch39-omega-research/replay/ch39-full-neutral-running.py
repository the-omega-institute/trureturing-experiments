"""八个自由截止的整条中性动能与色路径，并核查十二个截止的中间能标。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def weak_bound(data, sites, weights, axial, boson):
    v, mp, yt, g1, g2 = map(arb, ("246.22","1.22089e19",".8616",".467453",".63811"))
    r = (v/mp)**2/2
    w,z,t = [k*r for k in (g2**2/2,(g2**2+3*g1**2/5)/2,yt**2)]
    p0,pi = data["p0"],arb.pi()
    # 从 M_Z 到任一端点，各正的无截止矩不超过 M_Z 端点的读数。
    def entropy(mass):
        a=mass/p0
        return a*(1+(1/(4*a)).log1p())
    bound = ((2*axial["J"](w/p0)+boson["j_neutral"](z/p0))/(4*pi)
             +(entropy(w)+entropy(z))/(8*pi)
             +3*axial["J"](t/p0)/(8*pi)+axial["L"](t/p0)/(24*pi))
    coeff = (2*w+z)/(24*pi)+(w+z)/(48*pi)+t/(16*pi)+t/(24*pi)
    # 精确归一使其余权和等于 1-c[0]。极小节点单用原核的两端绝对界。
    return (abs(1-weights[0])*bound+
            coeff*sum(abs(c)/x for c,x in zip(weights[1:],sites[1:],strict=True))+
            4*abs(weights[0])*sites[0]/p0)


def direct(data, sites, weights, species, energy, gauge, rule):
    p0=data["p0"]
    p=(arb(energy)/arb("1.22089e19"))**2
    log=(p/p0).log()
    a1,a2,a3=data["initial"]
    out=[5*a1/3+a2,a2,a3]
    for x,c in zip(sites,weights,strict=True):
        for kind,z,b in species:
            exposure=gauge["angular_exposure"](kind,z,x,p0,p,log,rule)/(4*arb.pi())
            for j in range(3):
                out[j]-=c*b[j]*exposure
    return out


def path_floors(data, sites, weights, projected, initial, gauge):
    # 先按有符号系数对同一单调核分离，再只细分下界不足的动量段。
    p0,p1 = data["p0"],data["p1"]
    rule = gauge["angular_rule"](7,18)
    total_log = (p1/p0).log()
    entries = [(x,kind,z,[-c*b for b in beta]) for x,c in zip(sites,weights,strict=True)
               for kind,z,beta in projected if any(b!=0 for b in beta)]
    def endpoint(frac):
        log = total_log*frac
        p = p1 if frac==1 else p0*log.exp()
        pos,neg=[arb(0)]*3,[arb(0)]*3
        for x,kind,z,coeff in entries:
            term=gauge["angular_exposure"](kind,z,x,p0,p,log,rule)/(4*arb.pi())
            for j,b in enumerate(coeff):
                if b>0:
                    pos[j]+=b*term
                elif b<0:
                    neg[j]+=b*term
        return pos,neg
    floors=list(initial)
    start=([arb(0)]*3,[arb(0)]*3)
    stack=[(arb(0),arb(1),start,endpoint(arb(1)))]
    leaves=0
    while stack:
        lo,hi,left,right=stack.pop()
        bounds=[(a+p+n).lower() for a,p,n in zip(initial,left[0],right[1],strict=True)]
        if bounds[0]>15 and bounds[1]>40 and bounds[2]>2:
            floors=[min(f,b) for f,b in zip(floors,bounds,strict=True)]
            leaves+=1
            continue
        assert hi-lo>arb("1e-6"), (lo,hi,bounds)
        mid=(lo+hi)/2
        middle=endpoint(mid)
        stack.extend(((mid,hi,middle,right),(lo,mid,left,middle)))
    print("自适应动量闭段数",leaves,flush=True)
    return floors


def verify(root):
    build=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["verify"]
    sites,c8,c12,species,data=build(root)
    gauge=run_path(str(root/"ch39-all-field-signed-crossing.py"))
    rule=gauge["angular_rule"](12,33)
    vals=direct(data,sites,c12,species,"1e12",gauge,rule)
    print("十二位置，在 1e12 GeV 的 QQ、QT3、色",*[v.str(12) for v in vals],flush=True)
    assert arb("-79000") < vals[0] < arb("-77000") and vals[2]>0
    axial=run_path(str(root/"ch39-neutral-axial-top.py"))
    boson=run_path(str(root/"ch39-neutral-bosonic.py"))
    err=weak_bound(data,sites[:8],c8,axial,boson)
    print("八位置全程弱差绝对上界",err.upper().str(14),flush=True)
    a1,a2,a3=data["initial"]
    projected=[]
    ambiguous=arb(0)
    for kind,z,beta in species:
        # 代数上抵消的 beta 球可能含零；为这几行保留统一的双侧误差。
        values=[]
        for coefficient in (beta[1],beta[0]-beta[1],beta[2]):
            coefficient=arb(coefficient)
            if coefficient.contains(0):
                ambiguous+=abs(coefficient).upper()
                coefficient=arb(0)
            values.append(coefficient)
        projected.append((kind,z,values))
    correction=sum(abs(c) for c in c8)*ambiguous*(data["p1"]/data["p0"]).log()/(4*arb.pi())
    floors=[v-correction for v in path_floors(data,sites[:8],c8,projected,(a2,5*a1/3,a3),gauge)]
    print("QT3、QQ-QT3、色的沿程下界",*[v.str(13) for v in floors],flush=True)
    assert min(floors)>0
    b,y,colour=floors
    threshold=b*y/(b+y)
    assert err<threshold and colour>0

if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
