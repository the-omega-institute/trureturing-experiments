"""三循环字符块各自保迹时，第三枢轴在第四壳实现荷六五式。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


ANCHORS = (
    "8.92128431297471e-10", "6.21384726050197e-05", "8.89403672056501e-24",
    "1.57784599778674e-16", "1.99998055005156e-17", "0.080655233577267",
    "5.75149231390457e-18", "2.64152436436895e-26", "0.0559646514495884",
    "4.11126946866653e-34", "3.50979292048591e-05", "9.52513147234784e-35",
    "0.0232178217067686", "4.71512328802887e-06", "0.0343042184799948",
    "0.0187986513372345", "0.00166426052251678", "0.0176946404722149",
    "2.00298919733206e-23",
)
SELECTED = (17,18,15,16,5,12,0)
PARAMETER = 9


def angular_means(t):
    # I=∫exp(-t h)，H=∫h exp(-t h)，h=u(1-u)。
    # 低端用 1-t/4 <= exp(-t h) <=1，高端用 h>=u/2（半角）。
    if t<arb("1e-12"):
        return [arb(1).union(1-t/4)]*4
    if t>arb("1e8"):
        return [arb(0).union(12/t)]*4
    def point(z):
        integral=(arb.pi()/z).sqrt()*(-z/4).exp()*(z.sqrt()/2).erfi()
        first=((z+2)*integral-2)/(4*z)
        return [6*first,3*integral-12*first,(5*integral+12*first)/7,(9*integral+12*first)/11]
    # 四种角权都非负，故各均值对 t 递减。
    return [a.union(b) for a,b in zip(point(t.lower()),point(t.upper()),strict=True)]


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_089
    _,_,incoming,readout,(low,high,highest),_,thermal=lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    moments,heat,fields=thermal
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,number) if row[0]==6)
    weights=[z.real for z in pivots[2]]
    assert sum(weights).contains(1)
    sites=lib("three-ray-chiral-boundary")["geo"]["sites"]
    cycle=[sites.index((b,c,a)) for a,b,c in sites]
    groups=lib("charge-six-projector-algebra-transport")["groups"](readout.conjugate().transpose(),cycle)
    overlap=readout*incoming*acb_mat([[(high-low).sqrt(),0],[0,(highest-low).sqrt()]])
    totals,maxima=[],[]
    for group in groups:
        block=acb_mat([[overlap[i,j] for j in range(2)] for i in group])
        gram=block.conjugate().transpose()*block
        tr=gram.trace().real
        gap=((gram[0,0].real-gram[1,1].real)**2+4*abs(gram[0,1])**2).sqrt()
        assert tr>gap>0
        totals.append(len(group)*low+tr)
        maxima.append(low+(tr+gap)/2)
    gauge=lib("all-field-signed-crossing")
    neutral=lib("full-neutral-free-cutoff")
    data=gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=number)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=256
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    rows=[(kind,z,(qq-arb(8)*colour/3,mix-colour,colour)) for kind,z,(qq,mix,colour) in species]
    p0,p1=data["p0"],data["p1"]
    a1,a2,a3=data["initial"]
    pi=arb.pi()
    rule=gauge["angular_rule"](19,54)
    log=(p1/p0).log()
    anchors=list(map(arb,ANCHORS))
    center=[anchors[i] for i in SELECTED]
    def values(x):
        columns=[arb(0)]*3
        for kind,z,coef in rows:
            val=gauge["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*pi)
            for j,c in enumerate(coef):
                columns[j]-=c*val
        return [*columns[:2],*moments(fields,heat,x),columns[2]]
    def jacobian(coords):
        cols=[]
        for v,x0,i in zip(coords,center,SELECTED,strict=True):
            x=x0*(1+v)
            diff=[u-d for u,d in zip(angular_means(p0/x),angular_means(p1/x),strict=True)]
            gauge_cols=[arb(0),arb(0)]
            for kind,z,coef in rows:
                val=(-z/x).exp()*diff[kind]/(4*pi)
                for j in range(2):
                    gauge_cols[j]-=coef[j]*val
            ex=[(-z/x).exp() for _,_,z in fields]
            k=x*sum(a*e for (_,a,_),e in zip(fields,ex,strict=True))/(2*pi)
            w=-x*x*sum(g*e for (g,_,_),e in zip(fields,ex,strict=True))/(4*pi)
            trace=x.sqrt()/(2*(1+v))
            cols.append([weights[i]*val/(1+v) for val in (*gauge_cols,k,w)]+
                        [trace if i in group else arb(0) for group in groups])
        return arb_mat(cols).transpose()

    cols=[values(x) for x in anchors]
    f=[t+sum(w*c[j] for w,c in zip(weights,cols,strict=True))
       for j,t in enumerate((5*a1/3+a2-8*a3/3,a2-a3,-1,-arb("1e-122")))]
    f += [sum(anchors[i].sqrt() for i in group)-total for group,total in zip(groups,totals,strict=True)]
    span=arb("1e-9")
    param=anchors[PARAMETER]*(1+arb(0,span.upper()))
    delta=-(1-span).log()
    bounds=[sum(abs(b[j]) for _,_,b in rows)/(4*pi) for j in range(3)]
    derivs=[*bounds[:2],param*sum(abs(a) for _,a,_ in fields)/(2*pi),
            param**2*sum(abs(g) for g,_,_ in fields)/(4*pi)]
    errors=[abs(weights[PARAMETER])*d*delta for d in derivs]+[
        param.sqrt()*delta/2 if PARAMETER in group else arb(0) for group in groups]
    f=arb_mat([[v+arb(0,e.upper())] for v,e in zip(f,errors,strict=True)])
    radius=arb("1e-5")
    box=[arb(0,radius.upper())]*7
    pre=jacobian([arb(0)]*7).mid().inv()
    assert not pre.det().contains(0)
    offset=max(abs(v).upper() for v in (pre*f).entries())/radius
    identity=arb_mat([[int(i==j) for j in range(7)] for i in range(7)])
    defect=identity-pre*jacobian(box)
    slope=max(sum(abs(v).upper() for v in row) for row in defect.tolist())
    print("七维中心位移占比、导数余量行和",offset.str(14),slope.str(14),flush=True)
    assert offset+slope<arb(".019")
    output=anchors[:]
    output[PARAMETER]=param
    for i in SELECTED:
        output[i]*=1+box[0]
    amps=[x.sqrt() for x in output]
    for group,maximum in zip(groups,maxima,strict=True):
        assert all(low<amps[i]<maximum for i in group)
    weak=native["weak_series"](56)
    def weak_at(x):
        return weak(x/p0) if x<10**7*p0 else neutral["weak_column"](x,p0,p1)
    rest=sum(w*weak_at(x) for i,(w,x) in enumerate(zip(weights,output,strict=True)) if i!=PARAMETER)
    ends=[rest+weights[PARAMETER]*weak_at(anchors[PARAMETER]*(1+j*span)) for j in (-1,1)]
    print("弱差左右端点的下界与上界",*[(v.lower().str(13),v.upper().str(13)) for v in ends],flush=True)
    assert ends[0]>0>ends[1]
    colour=a3+sum(w*c[4] for w,c in zip(weights,cols,strict=True))
    colour+=arb(0,(bounds[2]*(sum(abs(weights[i]) for i in SELECTED)*(-(1-radius).log())+
                              abs(weights[PARAMETER])*delta)).upper())
    assert 34<colour<36
    print("色逆耦合上下界、三块迹与最大本征值",colour.lower().str(13),colour.upper().str(13),[v.str(13) for v in totals],
          [v.str(13) for v in maxima],flush=True)
    return amps,colour


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
