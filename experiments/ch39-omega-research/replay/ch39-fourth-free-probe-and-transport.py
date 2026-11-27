"""第四壳同时重选荷六输运和归一交叉权，用六个加权截止配平引力与中性规范。"""
from pathlib import Path
from runpy import run_path

from flint import arb,arb_mat,ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    number = 2_306_888_100_089
    _,ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    source = [x.sqrt() for x in ordered]
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    gravity = lib("three-index-free-atoms")
    moment,heat,(fields,) = gravity["broken_fields"](root,[number])
    data = gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec = 256
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    species,_ = neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    p0,p1 = data["p0"],data["p1"]
    a1,a2,a3 = data["initial"]
    starts = (5*a1/3+a2-8*a3/3,a2-a3)
    xs = [arb("1.31")*p0,*map(arb,("1e-28","1e-12","1e-6",".18",".45"))]
    amps = [x.sqrt() for x in xs]
    fill = (sum(source)-sum(amps))/13
    output = amps+[fill]*13
    gaps = lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,output)
    assert min(gaps)>0 and source[0]<fill<source[-1]

    rule = gauge["angular_rule"](15,42)
    low = lib("native-rh-ray-pvm")["weak_series"](56)
    coeffs = [(qq-arb(8)*colour/3,mixed-colour,colour) for _,_,(qq,mixed,colour) in species]
    cols,colours = [],[]
    for x in xs:
        exposure = [gauge["angular_exposure"](kind,z,x,p0,p1,(p1/p0).log(),rule)/(4*arb.pi())
                    for kind,z,_ in species]
        def shift(j):
            return -sum(b[j]*v for b,v in zip(coeffs,exposure,strict=True))
        d = low(x/p0) if x<256*p0 else neutral["weak_column"](x,p0,p1)
        k,w = moment(fields,heat,x)
        cols.append([1,k,w,d,shift(0),shift(1)])
        colours.append(shift(2))

    matrix = arb_mat(cols).transpose()
    target = arb_mat([[z] for z in (1,1,arb("1e-122"),0,-starts[0],-starts[1])])
    weights = gravity["solve"](matrix,target)
    tv = sum(abs(c) for c in weights)
    colour = a3+sum(c*g for c,g in zip(weights,colours,strict=True))
    print("最小主序余量、六权总变差、端点色逆耦合",*[v.str(14) for v in (min(gaps),tv,colour)],flush=True)
    print("六权",*[v.str(9) for v in weights],flush=True)
    assert arb("34.50")<tv<arb("34.51") and arb("37.61")<colour<arb("37.65")
    assert all(w*s>0 for w,s in zip(weights,(1,-1,-1,1,-1,1),strict=True))
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    caps = [sum(v*v for v in gram[j]).sqrt()/abs(gram[j][j]) for j in (1,2)]
    assert caps[0]<tv<caps[1]
    # 在未加权的第7、8层分别补足两个探针的范数。归一内积及全部投影权有显式实现。
    rho_probe = 1/caps[1]
    norm_slack = 1-rho_probe*tv
    assert norm_slack>0
    return weights


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
