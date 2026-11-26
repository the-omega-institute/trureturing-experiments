"""荷六首枢轴的正保序 Hill 截止同时匹配两项引力矩和弱差。"""
from math import factorial
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, ctx


def verify(root):
    degree = 6000
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    _,layers=run_path(str(root/"ch39-three-index-free-atoms.py"))["spectra"](root)
    rows=layers[0]
    heat=run_path(str(root/"ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    moments=run_path(str(root/"ch39-signed-cutoff-gravity.py"))["moments"]
    weak=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    for q,levels,pivots,_ in native["shell_weights"](root):
        if q==6:
            weights=[z.real for z in pivots[0]]
            break
    ctx.prec=256
    pi=arb.pi()
    mp,mz=arb("1.22089e19"),arb("91.1876")
    p0,p1=(mz/mp)**2,(arb("1.176397e14")/mp)**2
    relative=[(v/levels[17])**degree for v in levels]
    relative[17]=arb(1)
    r18=levels[17]/levels[18]
    weak_low=native["weak_series"]()
    def nodes(h,par):
        k,t=par
        y=[v/(v+t) for v in relative]
        return [(h+k*v,v,-k*v/(u+t)) for u,v in zip(relative,y,strict=True)]
    def values(h,par):
        node=nodes(h,par)
        cols=[moments(rows,heat,x) for x,_,_ in node]
        return arb_mat([[sum(c*v[j] for c,v in zip(weights,cols,strict=True))-b]
                        for j,b in enumerate((1,arb("1e-122")))])
    def jacobian(h,par):
        out=[[arb(0)]*2 for _ in range(2)]
        for c,(x,dk,dt) in zip(weights,nodes(h,par),strict=True):
            waves=[(-z/x).exp() for g,a,z in rows]
            k=sum(a*v for (_,a,_),v in zip(rows,waves,strict=True))/(2*pi)
            w=-x*sum(g*v for (g,_,_),v in zip(rows,waves,strict=True))/(4*pi)
            for i,f in enumerate((k,w)):
                out[i][0]+=c*f*dk
                out[i][1]+=c*f*dt
        return arb_mat(out)
    low=p0*arb(1).union(2)
    center=[arb(".635"),arb("1.868")]
    for _ in range(5):
        step=jacobian(low.mid(),center).solve(values(low.mid(),center))
        center=[(x-v.mid()).mid() for x,v in zip(center,step.entries(),strict=True)]
    radius=arb("1e-10")
    box=[v+arb(0,radius.upper()) for v in center]
    assert 0<box[0]<1-low and box[1]>0
    pre=jacobian(low.mid(),center).inv()
    residual=pre*values(low,center)
    defect=arb_mat([[1,0],[0,1]])-pre*jacobian(low,box)
    row_norm=[sum(abs(z).upper() for z in row) for row in defect.tolist()]
    assert all(abs(z)+t*radius<radius for z,t in zip(residual.entries(),row_norm,strict=True))
    def D(t):
        terms=[]
        for x,_,_ in nodes(t*p0,box):
            if x < arb("1e-24"):
                term=weak_low(x/p0)
            else:
                term=weak(x,p0,p1)
            terms.append(term)
        return sum(c*d for c,d in zip(weights,terms,strict=True)),sum(abs(c)*d.rad() for c,d in zip(weights,terms,strict=True))
    left,right=arb(1),arb(2)
    dl,e1=D(left);dr,e2=D(right)
    assert dl>0>dr,(degree,dl,dr)
    for _ in range(22):
        m=(left+right)/2
        val,err=D(m)
        if val>0:
            left=m
        elif val<0:
            right=m
        else:
            break
    dl,e1=D(left);dr,e2=D(right)
    assert dl>0>dr
    h=p0*left.union(right)
    x=nodes(h,box)
    assert arb("5.66600e-35")<h<arb("5.66604e-35")
    assert arb(".635595244")<box[0]<arb(".635595246")
    s=box[1]*r18**degree
    assert arb("3.291030650e-40")<s<arb("3.291030654e-40")
    qmin=mp*x[0][0].sqrt()
    assert arb("91.90004")<qmin<arb("91.90006")
    assert all(a<b for a,b in zip(levels[:-1],levels[1:],strict=True))
    assert max(row_norm)<arb("2.2e-9") and max(abs(v).upper() for v in residual.entries())<arb("1.8e-26")
    trunc=10*arb(4)**33/(33*factorial(33))
    error=arb(0)
    for c,(position,_,_) in zip(weights,nodes(low,box),strict=True):
        if position < arb("1e-24"):
            error+=abs(c)*(trunc+2*position/p1+arb("1.1e-21"))
        else:
            error+=abs(c)*(arb("2.400e-21")+arb("7.468e-36")/position)
    print("全参数盒十九项弱差核加权截断界",error.upper().str(12),flush=True)
    assert error<arb("8e-18") and max(e1,e2)<arb("8e-18")
    print("Hill 指数、h、κ、t=s/r18^d、s",degree,h.str(12),box[0].str(14),box[1].str(14),
          (box[1]*r18**degree).str(14),flush=True)
    print("最小 Q/GeV、十六至十九层 x、D 两端与最大加权弱核半宽",
          (mp*x[0][0].sqrt()).str(12),*[v[0].str(12) for v in x[15:]],dl.str(9),dr.str(9),
          max(e1,e2).str(6),flush=True)
    print("引力预条件中心残差和两行雅可比缺陷",*[v.str(7) for v in residual.entries()],
          *[v.str(7) for v in row_norm],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
