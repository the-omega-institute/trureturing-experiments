"""右壳第十八层像在保持逐级荷六振幅时的相位自由度及奇夸克所需射线位移。"""
from itertools import combinations
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def tilted_phases(base, derivative, waves, ids, picks, phase13):
    # 第16层只进入有界正弦与余弦，覆盖任意相位。
    def system(xs, wide):
        phase = [arb(0)]*7
        phase[ids.index(12)] = arb(phase13) + (arb(0,arb(".000001").upper()) if wide else 0)
        for i,x in zip(picks,xs,strict=True):
            phase[i] = x
        trig = [x.sin_cos() for x in phase]
        if wide:
            trig[ids.index(15)] = (arb(0,1),arb(0,1))
        residual = arb_mat([[sum(b[k]*(c-1)+d[k]*s for b,d,(s,c) in zip(base,derivative,trig,strict=True))]
                            for k in range(5)])
        jac = arb_mat([[-base[i][k]*trig[i][0]+derivative[i][k]*trig[i][1] for i in picks] for k in range(5)])
        return residual,jac,trig

    point = [arb(0)]*5
    for _ in range(10):
        val,jac,_ = system(point,False)
        point = [(x-y).mid() for x,y in zip(point,jac.mid().solve(val.mid()).entries(),strict=True)]
    radius = arb(".0025")
    box = [x+arb(0,radius.upper()) for x in point]
    val,j0,_ = system(point,True)
    pre = j0.mid().inv()
    _,jac,trig = system(box,True)
    defect = arb_mat([[int(i==j) for j in range(5)] for i in range(5)])-pre*jac
    factor = max((abs(v)/radius+sum(abs(x) for x in row)).upper()
                 for v,row in zip((pre*val).entries(),defect.tolist(),strict=True))
    assert factor<1, factor
    phases = [acb(c,s) for s,c in trig]
    shifted = sum((v*(z-1) for v,z in zip(waves,phases,strict=True)),acb_mat(19,1))
    weights = [(v.conjugate().transpose()*v)[0,0].real for v in waves]
    gap = sum(weights[i]*weights[j]*abs(phases[i]-phases[j])**2 for i in range(7) for j in range(i))
    assert gap>arb(".062")**2
    return shifted,factor,gap.sqrt().lower()


def check_strange_floor(source, target, labels, fixed, images):
    adj = lambda z:z.conjugate().transpose()
    origin = adj(source)*fixed
    rhs = adj(images)*target
    shifts = {0:".00017",5:".0013",6:".0128",9:".0208",11:".00041",15:".0023",18:".082"}
    floors = [arb(0)]*19
    for j,value in shifts.items():
        shift = arb(value)
        scaled = acb_mat([[origin[i,k]/(labels[i]+shift) for k in range(4)] for i in range(19)])
        metric = adj(origin)*scaled
        data = acb_mat([[rhs[i,j]] for i in range(4)])
        trial = metric.mid().solve(data.mid()).mid()
        weighted = origin*trial
        loss = sum(abs(weighted[i,0])**2/(labels[i]+shift) for i in range(19))
        bound = 2*(adj(trial)*data)[0,0].real-loss-shift
        floors[j] = max(arb(0),bound.lower())
    rt3 = arb(3).sqrt()
    rows = acb_mat([[target[14,j] for j in range(19)],
                    [sum(target[i,j] for i in (8,7,6))/rt3 for j in range(19)]])
    block = acb_mat([[rows[i,j]*floors[j] for j in range(19)] for i in range(2)])*adj(rows)
    for i in range(2):
        block[i,i] -= arb(".00160")
    assert block[0,0].real>0 and block.det().real>0


def verify(root):
    joint = run_path(str(root / "ch39-one-generator-up-and-gravity.py"))
    rank = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    shells = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    orbits = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    _, (a,_), gev, (low,high,_) = joint["verify"](root,return_parameters=True)
    ctx.prec = 384
    adj, col = rank["adj"],joint["col"]
    n = rank["N"] + 1
    source = spectrum(shells,6,n)[1]
    target = spectrum(shells,6,n+1)[1]
    E = acb_mat([[1/arb(len(o)).sqrt() if i in o else 0 for o in orbits] for i in range(19)])
    seed = acb_mat([[acb(x,y) for x,y in zip(re,im,strict=True)]
                    for re,im in zip(rank["REAL"],rank["IMAG"],strict=True)])
    kp = E*seed*(1/arb(10000))
    P = kp*(adj(kp)*kp).inv()*adj(kp)
    S = col(source,17)
    p = P*S
    s = S-p
    _,offset,sm,bg,tangent,length = joint["mode_parameters"](P,E,source,target)
    e = ((gev/arb("1.22089e19")-low)/(high-low)).sqrt()
    lin = offset+sm*e+bg*a
    rad = (length-(adj(lin)*lin)[0,0].real).sqrt()
    mid = offset+sm*e.mid()+bg*a.mid()
    w = mid+tangent*(length-(adj(mid)*mid)[0,0].real).sqrt()
    for var,direction in ((e,sm),(a,bg)):
        w += (direction-tangent*((adj(lin)*direction)[0,0].real/rad))*(var-var.mid())
    Y = w-p*acb(0,1)
    ids = [j for j in range(19) if (adj(col(target,j))*E*adj(E)*col(target,j))[0,0].real>arb(".9")]
    assert len(ids)==7 and 15 in ids
    waves = [col(target,j)*(adj(col(target,j))*Y)[0,0] for j in ids]
    waves[ids.index(15)] = joint["align"](col(target,15),S)*e
    derivative = []
    for wave in waves:
        v = wave*acb(0,1)
        c = (adj(kp)*v).entries()
        derivative.append([c[0].real,c[0].imag,c[1].real,c[1].imag,(adj(s)*v)[0,0].imag])
    candidates = []
    for picks in combinations([i for i,j in enumerate(ids) if j!=15],5):
        det = arb_mat([derivative[i] for i in picks]).det()
        candidates.append((abs(det).lower(),picks))
    gap,picks = max(candidates,key=lambda pair:pair[0])
    assert gap>0
    free = [j for i,j in enumerate(ids) if i not in picks]
    assert all(abs((adj(col(target,j))*Y)[0,0])>0 for j in free if j!=15)
    assert abs((adj(col(target,15))*Y)[0,0]).overlaps(e)
    assert e>0 and (adj(p)*p)[0,0].real>arb(".04507")
    print("五个被解相位层、行列式绝对值下界、两个自由相位层",
          [ids[i]+1 for i in picks],gap.str(12),[j+1 for j in free])
    # X 的实际第十八层像靠近参考射线时，从其在 P 和 span(s,w) 上的作用估计所需距离。
    pnorm = (adj(p)*p)[0,0].real.sqrt()
    snorm = (adj(s)*s)[0,0].real.sqrt()
    inner = abs((adj(s)*w)[0,0])
    least = snorm**2-inner
    assert least>arb(".077")
    phase_amplification = 1+snorm/pnorm
    core_error = max((1/pnorm).upper(),(arb(2)/least).sqrt().upper()*phase_amplification.upper())
    prior_radius = arb(".0010818")
    chord_floor = prior_radius/(2*core_error)
    assert chord_floor>arb(".00001906")
    print("补空间 Gram 最小本征值下界、四维核心误差系数上界、射线弦距离下界",
          least.lower().str(12),core_error.str(12),chord_floor.lower().str(12))
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    _,labels,_ = train(root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    down_source,down_target = [spectrum(shells,2,i)[1] for i in (n,n+1)]
    base = []
    for v in waves:
        x,y = (adj(kp)*v).entries()
        base.append([x.real,x.imag,y.real,y.imag,(adj(s)*v)[0,0].imag])
    first = [col(kp,0),col(kp,1)]
    sign_line = p*adj(p)*(1/(adj(p)*p)[0,0].real)
    for phase in ("-.5",".5"):
        change,factor,separation = tilted_phases(base,derivative,waves,ids,picks,phase)
        w_new = w+change
        fixed = joint["side_by_side"]([*first,s,w_new])
        for action in (P,2*sign_line-P):
            images = joint["side_by_side"]([*(action*v*acb(0,1) for v in first),w_new,s])
            check_strange_floor(down_source,down_target,labels[n,2],fixed,images)
        print("第13层相位中心、统一压缩比、像射线正弦距离下界",
              phase,factor.str(10),separation.str(12))
    reference,_ = run_path(str(root/"ch39-fixed-reflection-quark-spectral-neighbourhood.py"))["reflection_spectrum"](root)
    window = [arb(x).union(arb(y)) for x,y in rank["WINDOWS"]]
    char_mass = sum(prod(v-m for m in window)/prod(v-z for k,z in enumerate(reference) if k!=i)
                    for i,v in enumerate(reference[:2]))
    assert arb(".00160")*char_mass>arb(".0002793")
    print("新相位两族在上型三窗内的奇质量下界",(arb(".00160")*char_mass).lower().str(12))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
