"""枚举参考上型反射与一条荷六高源像射线所容许的全部四维核心作用。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify(root):
    rank = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    joint = run_path(str(root / "ch39-one-generator-up-and-gravity.py"))
    prior = run_path(str(root / "ch39-fixed-core-light-down-mass-floors.py"))
    train = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    shells = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    orbits = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    adj, col = rank["adj"], joint["col"]
    columns = joint["side_by_side"]
    _, (a,_), gev, (low,high,_) = joint["verify"](root,return_parameters=True)
    n = rank["N"] + 1
    _, label, _ = train(root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    six_source = spectrum(shells,6,n)[1]
    six_target = spectrum(shells,6,n+1)[1]
    source = spectrum(shells,2,n)[1]
    target = spectrum(shells,2,n+1)[1]
    sites = (shells["soft"]["pole"],*shells["boundary_sites"]())
    E = acb_mat([[1/arb(len(o)).sqrt() if i in o else 0 for o in orbits] for i in range(19)])
    seed = acb_mat([[acb(x,y) for x,y in zip(re,im,strict=True)]
                    for re,im in zip(rank["REAL"],rank["IMAG"],strict=True)])
    kp = E*seed*(arb(1)/10000)
    P = kp*(adj(kp)*kp).inv()*adj(kp)
    S = col(six_source,17)
    p = P*S
    pnorm2 = (adj(p)*p)[0,0].real
    assert pnorm2>arb(".04507")
    line = p*adj(p)*(1/pnorm2)
    s = S-p
    _, offset, small, big, tangent, length = joint["mode_parameters"](P,E,six_source,six_target)
    e = ((gev/arb("1.22089e19")-low)/(high-low)).sqrt()
    linear = offset + small*e + big*a
    radius = (length-(adj(linear)*linear)[0,0].real).sqrt()
    mid = offset + small*e.mid() + big*a.mid()
    w = mid + tangent*(length-(adj(mid)*mid)[0,0].real).sqrt()
    for var, direction in ((e,small),(a,big)):
        w += (direction-tangent*((adj(linear)*direction)[0,0].real/radius))*(var-var.mid())
    gram_det = (adj(s)*s)[0,0].real*(adj(w)*w)[0,0].real-abs((adj(s)*w)[0,0])**2
    assert gram_det>arb(".14")
    first = [col(kp,0),col(kp,1)]
    fixed = columns([*first,s,w])
    floors = []
    for action in (P,2*line-P):
        phase = action*acb(0,1)
        images = columns([*(phase*v for v in first),w,s])
        level = prior["destination_trivial_floor"](source,target,sites,label[n,2],fixed,images,adj)
        floors.append(level)
        print("秩二子空间上的分支、下型平凡块严格下界",
              "同号" if action is P else "异号",level.str(8))

    # 两分支保持同一个四维平凡核心。重排界只依赖其正交补上的源压缩谱。
    char_sum, _ = prior["character_sum_bound"](root)
    assert all(0<x<arb(".293784") for x in label[n,2]) and char_sum>arb(".00058130")
    strange = char_sum/2
    distance = (strange.sqrt()-arb(".000271").sqrt())/arb(".293784").sqrt()
    assert min(floors)>strange and distance>arb(".0010818")
    print("源投影模方、源像补空间的 Gram 行列式",pnorm2.str(12),gram_det.str(8))
    print("全部上型反射平方根分支的奇质量下界、发生元距离下界",
          strange.lower().str(13),distance.lower().str(13))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
