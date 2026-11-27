"""在三条最短弦上旋转第五壳完整联合算子，检验第七壳原生第三枢轴双矩。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


SEVENTH = 2_306_888_100_092


def adj(matrix):
    return matrix.conjugate().transpose()


def geometry(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    shell = lib("three-ray-chiral-boundary")
    _, accepted = lib("behrend-rational-chord-pairs")["inspect"](shell["geo"]["sites"], 867)
    assert accepted == [(6,11,12),(7,12,12),(8,18,12)]
    geo = shell["geo"]
    pole = (17,17,17)
    for i,j,_ in accepted:
        a,b = geo["sites"][i],geo["sites"][j]
        assert geo["det"](pole,a,b)==0
        assert sum((x-y)**2 for x,y in zip(pole,a,strict=True))<sum((x-y)**2 for x,y in zip(pole,b,strict=True))
    ctx.prec = 320
    matrix = shell["laplacian"](6, acb.zeta_zero(SEVENTH).imag)
    assert all((-matrix[i,j]).contains(1) for i,j,_ in accepted)
    eigenvalues, basis = lib("adjacent-pvm-transport")["spectrum"](shell,6,SEVENTH)
    _, nodes, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,SEVENTH) if row[0]==6)
    assert all(a.overlaps(b) for a,b in zip(eigenvalues,nodes,strict=True))
    weights = [v.real for v in pivots[2]]
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root,[SEVENTH])
    ctx.prec = 192
    assert 0 < levels[0] < levels[1] < levels[2]
    station = fifth*columns
    baseline = (adj(basis)*fifth)*columns
    cosine, sine = [], []
    for a,b,_ in accepted:
        cosine.append(acb_mat([[basis[a,i].conjugate()*station[a,j]+basis[b,i].conjugate()*station[b,j]
                                for j in range(2)] for i in range(19)]))
        sine.append(acb_mat([[basis[b,i].conjugate()*station[a,j]-basis[a,i].conjugate()*station[b,j]
                              for j in range(2)] for i in range(19)]))
    return baseline,cosine,sine,basis,station,weights,heat,field


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    baseline,cosine,sine,_,_,weights,heat,field = geometry(root,fifth,columns,levels)
    low,h,k = levels
    increases = (h-low,k-low)
    fixed = baseline-cosine[0]-cosine[1]-cosine[2]
    def amplitudes(projection):
        return [low+sum(w*abs(projection[i,j])**2 for j,w in enumerate(increases)) for i in range(19)]

    def angles(theta):
        parts = [a.sin_cos() for a in theta]
        z = fixed
        directions = []
        for (s,c),u,v in zip(parts,cosine,sine,strict=True):
            z = z+c*u+s*v
            directions.append((-s)*u+c*v)
        return z,[directions[0],directions[2]]

    def bands(amps):
        return [(max(a.lower(),low.lower())**2,min(a.upper(),k.upper())**2) for a in amps]

    def physical_heat(bounds):
        totals = [arb(0),arb(0)]
        for (left,right),weight in zip(bounds,weights,strict=True):
            assert 0<left and left<=right
            for g,alpha,mass in field:
                m2 = (left*heat(2,mass/left)).union(right*heat(2,mass/right))
                m4 = (left**2*heat(3,mass/left)).union(right**2*heat(3,mass/right))
                totals[0] += weight*alpha*m2/(2*arb.pi())
                totals[1] -= weight*g*m4/(4*arb.pi())
        return totals

    def eqs(center, want_values):
        projection, deriv = angles([center[0],5*arb.pi()/8,center[1]])
        amps = amplitudes(projection)
        bounds = bands(amps)
        values = physical_heat(bounds) if want_values else [arb(0),arb(0)]
        values[0] -= 1
        values[1] -= arb("1e-122")
        jac = [[arb(0),arb(0)],[arb(0),arb(0)]]
        for i,(a,(left,right)) in enumerate(zip(amps,bounds,strict=True)):
            assert 0<left and left<=right
            dk, dw = arb(0),arb(0)
            for g,alpha,mass in field:
                ex = (-mass/left).exp().union((-mass/right).exp())
                dk += alpha*ex
                dw += g*ex
            slopes = (a*dk/arb.pi(),-a**3*dw/(2*arb.pi()))
            for col,d in enumerate(deriv):
                da = 2*sum(w*(projection[i,j].conjugate()*d[i,j]).real for j,w in enumerate(increases))
                for row,slope in enumerate(slopes):
                    jac[row][col] += weights[i]*slope*da
        return arb_mat([[v] for v in values]),arb_mat(jac)

    center = [arb("3.2176023525"),arb("-2.5142248392")]
    residual,j = eqs(center,True)
    pre = j.mid().inv()
    radius = arb(".0001")
    _,wide = eqs([a+arb(0,radius) for a in center],False)
    unit = arb_mat([[1,0],[0,1]])
    defect = max(sum(abs(v).upper() for v in row) for row in (unit-pre*wide).tolist())
    offset = max(abs(v).upper() for v in (pre*residual).entries())
    ratio = offset/radius+defect
    print("独立站位角的自映比与盒内导数行列式",ratio.upper().str(15),wide.det().str(14),flush=True)
    assert ratio<arb(".169") and wide.det()>arb("2.18")

    # 同角时逐线重叠是常数、正弦、余弦的线性组合。以统一的质量修正包围所有热矩。
    assert k<arb(".75553")
    normalized = lib("fifth-seventh-cutoff-separation")["moment_errors"](field,heat,weights,[".75553"]*19)
    slopes = (sum(a for _,a,_ in field)/(2*arb.pi()),-sum(g for g,_,_ in field)/(8*arb.pi()))
    csum = cosine[0]+cosine[1]+cosine[2]
    ssum = sine[0]+sine[1]+sine[2]
    u,v = arb(168)/457,arb(425)/457
    assert 168**2+425**2==457**2
    fixed_bands = bands(amplitudes(fixed+u*csum+v*ssum))
    net, vacuum = physical_heat(fixed_bands)
    colour = lib("fifth-joint-coherent-holdout")["colour_probe"](root,SEVENTH,fixed_bands,weights)
    print("近端向远端的固定有理弦距门：牛顿、真空、56 万 GeV 色项",
          net.str(16),vacuum.str(16),colour.str(16),flush=True)
    assert arb("4.17077")<net<arb("4.17083")
    assert arb("-1.84768")<vacuum<arb("-1.84766")
    assert arb("18.33768")<colour<arb("18.33772")
    # 逐种在振幅上端控制有质量与无质量热矩的差。损失随平方截止增加。
    loss_cache = {}
    def loss(tick):
        if tick not in loss_cache:
            x = (arb(tick)/1000)**2
            if not tick:
                loss_cache[tick] = (arb(0),arb(0))
            else:
                loss_cache[tick] = (x*sum(abs(alpha)*(1-heat(2,z/x)) for _,alpha,z in field)/(2*arb.pi()),
                                    x*x*sum(abs(g)*(arb(1)/2-heat(3,z/x)) for g,_,z in field)/(4*arb.pi()))
        return loss_cache[tick]
    pending = [(arb(0),2*arb.pi())]
    covered = 0
    least = arb(100)
    while pending:
        l,r = pending.pop()
        s,c = l.union(r).sin_cos()
        amps = amplitudes(fixed+c*csum+s*ssum)
        gaps = [abs(slope*(sum(w*a**p for w,a in zip(weights,amps,strict=True))-target)).lower()
                for p,slope,target in zip((2,4),slopes,normalized,strict=True)]
        gap = max(gaps)
        if gap<=arb(".01") and r-l<arb(".02"):
            caps = [int((min(k.upper(),a.upper())*1000).ceil().unique_fmpz()) for a in amps]
            errors = [sum(abs(w)*loss(t)[j] for w,t in zip(weights,caps,strict=True)) for j in range(2)]
            target = (1,arb("1e-122"))
            gap = max((abs(slopes[j]*sum(w*a**p for w,a in zip(weights,amps,strict=True))-target[j])
                       -errors[j]).lower() for j,p in enumerate((2,4)))
        if gap>arb(".01"):
            covered += 1
            least = min(least,gap)
        else:
            assert r-l>arb("1e-9"),(l,r,gaps)
            m = (l+r)/2
            pending.extend(((l,m),(m,r)))
    print("同角完整周期的闭段数、两残差绝对值较大者的统一下界",covered,least.str(15),flush=True)
    assert least>arb(".01022")
    return ratio,least


def verify(root):
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
