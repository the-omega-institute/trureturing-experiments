"""固定第四壳十九个截止及能标，求同时配平六行的有符号实权最小总变差。"""
from math import comb, factorial
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def response_matrix(root,cutoffs):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    n = 2_306_888_100_089
    moments,heat,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[n])
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=n)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,n)
    ctx.prec = 256
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    species,_ = neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    assert len(species)==119 and len(cutoffs)==19
    p0,p1=data["p0"],data["p1"]
    log=(p1/p0).log()
    init1,init2,init3=data["initial"]
    start=(5*init1/3+init2-8*init3/3,init2-init3)
    rows=[(kind,z,(arb(qq)-arb(8)*colour/3,arb(mixed)-colour)) for kind,z,(qq,mixed,colour) in species]
    pi4=4*arb.pi()
    # 四种角权在 h=u(1-u) 上分别是 6h、3−12h、(5+12h)/7、(9+12h)/11。
    angular=((0,6,1),(3,-12,1),(5,12,7),(9,12,11))
    order=40
    beta=[arb(factorial(j))**2/factorial(2*j+1) for j in range(order+2)]
    h_moments=[[(a*beta[j]+b*beta[j+1])/c for j in range(order+1)] for a,b,c in angular]
    coeff=[[arb(0)]*(order+1) for _ in start]
    # E1(a/x)−E1(b/x) 的对数项独立于 x。先合并各物种的整个反幂级数。
    for kind,z,charges in rows:
        powers=[z**j for j in range(order+1)]
        for j in range(1,order+1):
            difference=sum(comb(j,k)*powers[j-k]*(p1**k-p0**k)*h_moments[kind][k]
                           for k in range(1,j+1))
            term=(-1)**j*difference/(pi4*j*factorial(j))
            for r,c in enumerate(charges):
                coeff[r][j]-=c*term
    largest=max((z+p1/4).upper() for _,z,_ in rows)
    charge_size=[sum(abs(c[r]).upper() for _,_,c in rows) for r in range(2)]
    def corrections(x):
        t=largest/x
        assert t<2
        # 交错 E1 尾项按两个端点各自估计；四种角权非负且积分为一。
        tail=2*t**(order+1)/((order+1)*factorial(order+1)*pi4)
        out=[]
        for col,size in zip(coeff,charge_size,strict=True):
            v=arb(0)
            for a in reversed(col[1:]):
                v=(v+a)/x
            out.append(v+arb(0,(size*tail).upper()))
        return out
    def direct(x):
        nodes,quad_error=rule
        result=[arb(0),arb(0)]
        e_limit=(-arb(100)).exp()
        def e1(y):
            lower,upper=y.lower(),y.upper()
            if lower>100:
                return arb(0).union(e_limit/100)
            return lower.expint(1).union(upper.expint(1))
        for kind,z,charges in rows:
            if (z/x).lower()>100:
                exposure=arb(0).union(log*e_limit)
            else:
                exposure=sum(weights[kind]*(e1((z+p0*h)/x)-e1((z+p1*h)/x)) for h,weights in nodes)
                exposure+=arb(0,(quad_error*log).upper())
            for r,c in enumerate(charges):
                result[r]-=c*exposure/pi4
        return result
    rule=gauge["angular_rule"](16,43)
    # x→∞ 的角积分直接取对数。原 Gauss 误差对 x 一致，故也覆盖这一极限。
    inf=[arb(0),arb(0)]
    nodes,quad_error=rule
    for kind,z,charges in rows:
        constant=sum(weights[kind]*((z+p1*h)/(z+p0*h)).log() for h,weights in nodes)
        constant=(constant+arb(0,(quad_error*log).upper()))/pi4
        for r,c in enumerate(charges):
            inf[r]-=c*constant
    g16=direct(cutoffs[15])
    weak_low=lib("native-rh-ray-pvm")["weak_series"](56)
    ref=arb("-.074256577194663365545")
    d16=weak_low(cutoffs[15]/p0)
    ell_target=-ref/(d16-ref)
    target=arb_mat([[x] for x in (1,1,arb("1e-122"),ell_target,
                   *[-s-g-(v-g)*ell_target for s,g,v in zip(start,inf,g16,strict=True)])])
    cols=[]
    for i,x in enumerate(cutoffs):
        # E_k(v) ≤ exp(−v)/(k−1) 对 k>1 成立。
        def bounded_heat(k,v):
            if v>160:
                return arb(0).union((-arb(160)).exp()/(k-1))
            return heat(k,v)
        K,W=moments(field,bounded_heat,x)
        if i==15:
            ell=arb(1)
            dG=[arb(0),arb(0)]
        else:
            ell=(neutral["weak_column"](x,p0,p1)-ref)/(d16-ref)
            dG=[p-(v-g)*ell for p,v,g in zip(corrections(x),g16,inf,strict=True)]
        cols.append([1,K,W,ell,*dG])
    return arb_mat(cols).transpose(),target


def certify_minimum(matrix, target, original):
    support = (4, 6, 13, 15, 17, 18)
    signs = (1, -1, 1, 1, -1, 1)
    basis = arb_mat([[matrix[i,j] for j in support] for i in range(6)])
    assert not basis.det().contains(0)
    weights = basis.solve(target, algorithm="precond").entries()
    assert all(v*s>0 for v,s in zip(weights,signs,strict=True))
    sign_vector = arb_mat([[s] for s in signs])
    dual = basis.transpose().solve(sign_vector, algorithm="precond")
    products = (matrix.transpose()*dual).entries()
    # 对每个精确矩阵，支撑上依定义有 B^T λ = sign(w)，只需界住支撑外各列。
    inactive = max(abs(p).upper() for j,p in enumerate(products) if j not in support)
    assert inactive < arb(".997246")
    cost = sum(abs(v) for v in weights)
    assert cost.overlaps((target.transpose()*dual)[0,0])
    assert arb("3070445.7") < cost < arb("3070447.0")
    residual = (matrix*arb_mat([[v] for v in original])-target).entries()
    assert all(v.contains(0) for v in residual[:4])
    assert arb("502.659") < residual[4] < arb("502.661")
    assert arb("90.492") < residual[5] < arb("90.493")
    for j,v in zip(support,weights,strict=True):
        print("最优权，第",j+1,"层",v.lower().str(15),v.upper().str(15))
    print("最小总变差端点",cost.lower().str(15),cost.upper().str(15))
    print("非支撑列的最大对偶绝对值上界",inactive.str(12))
    full = dict(zip(support,weights,strict=True))
    lower = cost.lower()-sum(abs(v) for v in original).upper()
    upper = sum(abs(full.get(i,0)-v) for i,v in enumerate(original)).upper()
    assert lower > arb("3070432.3") and upper < arb("3070437.2")
    print("离原权的最小距离，下界与可行上界",lower.str(15),upper.str(15))
    return cost


def compare_sources(root, original, cost):
    pivot = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ctx.prec = 256
    gram = pivot["gram_cholesky"](root)
    _,_,direct,_ = next(row for row in pivot["shell_weights"](root,2_306_888_100_089) if row[0]==6)
    actual = [sum(abs(v.real) for v in weights) for weights in direct]
    assert actual[0].overlaps(sum(abs(v) for v in original))
    caps = [sum(v*v for v in gram[j]).sqrt()/abs(gram[j][j]) for j in (1,2,3)]
    assert all(a<c for a,c in zip(actual,caps,strict=True))
    assert cost > 1176*max(c.upper() for c in caps)
    print("第四壳实际三枢轴总变差",*[v.str(12) for v in actual])
    print("三枢轴与投影无关的总变差上界",*[v.upper().str(12) for v in caps])
    cross = run_path(str(root/"ch39-rh-cross-spectral-cutoff.py"))
    groups = cross["verify"]()
    pairs = [row for group in groups.values() for row in group]
    assert len(pairs)==388
    # 原枚举给负变差上界，归一实权的总变差等于一加负变差的两倍。
    cap = 1+2*max(high for _,high,_,_ in pairs)
    assert cap < Fraction(7690)
    assert cost > 399*7690
    print("388 对的总变差统一上界",cap)


def verify(root):
    joint = run_path(str(root/"ch39-fourth-full-ckm-gravity.py"))
    cutoffs,original = joint["verify"](root,return_cutoffs=True)
    matrix,target = response_matrix(root,cutoffs)
    cost = certify_minimum(matrix,target,original)
    compare_sources(root,original,cost)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
