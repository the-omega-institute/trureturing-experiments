"""保持荷六对合，以另一立方根让同一四站三平面同时达到三个上型窗与底质量窗。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


SEED = (
    ((-6943, 29910), (-131167, 22029), (-65610, 149895), (60809, 491595), (-83622, -59931), (4731, 204557), (109380, 293535), (-49257, 27148)),
    ((198322, 194495), (197137, 107878), (95291, 301072), (-114491, 111595), (-254219, -55086), (-98442, -108263), (-42346, 282099), (-153405, -164015)),
    ((-207985, 58792), (-25074, 75997), (-222175, -6918), (139787, 184699), (-237473, -128922), (-17895, -156812), (-302957, -141790), (-60584, -202698)),
    ((-57457, -28616), (179353, 256113), (314456, 90664), (258825, 60008), (-21613, -118634), (81779, 270052), (-131316, -83036), (52428, 48191)),
    ((-82885, 44393), (226816, -129848), (73127, -208867), (-213948, 104315), (-103456, -142984), (99975, 170579), (-337088, 313699), (19809, -50393)),
    ((112690, -225229), (113640, -100503), (203285, -127857), (-50662, 368358), (-181902, -77501), (-16123, -62781), (194448, -155542), (-73495, -102971)),
    ((229815, 306842), (-88351, -138095), (-27354, -112161), (-30545, -62202), (-122765, -220343), (-31534, 296711), (88296, -180295), (-132502, 101415)),
    ((236769, -224684), (-28115, -168623), (-134838, 113475), (157463, -37712), (-186127, -286692), (13750, -70644), (-141261, -19062), (149156, 251892)),
    ((-227023, 35057), (-84671, -205525), (114311, 93741), (84669, -212777), (-211906, -198645), (-327031, 73817), (113150, 99124), (138643, -67269)),
    ((13320, 306817), (119449, 253526), (-300657, -98149), (-37165, -1936), (-52390, 55635), (-198156, 49248), (192477, -77939), (180959, 246328)),
    ((-67744, 19786), (76411, 228755), (104071, 121608), (151560, 120530), (-184018, 1992), (-132111, -47429), (-155639, -41694), (-66909, 235888)),
    ((-340878, 62935), (52967, 26727), (13032, 90423), (-15665, 78705), (-3577, -204184), (407958, 273136), (108114, 36026), (67200, 143046)),
    ((-159461, 68990), (215122, 108547), (52336, 68375), (-10953, -67110), (-215927, -188359), (153214, -131214), (300549, -164253), (29823, -282366)),
    ((117, -107294), (71151, 117141), (-265438, -212801), (120115, 226656), (-175929, -75193), (-132820, 56339), (-26439, -44852), (11178, -138688)),
    ((120849, -94374), (-337148, 20011), (16560, -16878), (-108661, -240255), (-478289, -95722), (-4567, 169117), (6620, -26716), (65064, 75888)),
    ((111056, -135337), (412828, 68431), (-286095, 112324), (-166807, -140204), (-108896, -19854), (38222, -3766), (133924, -41675), (-18147, 249820)),
    ((-320729, 1098), (105841, -93095), (-111101, -5523), (-221729, 41203), (-146018, 66826), (-125879, -181968), (194606, -9099), (84430, 420506)),
    ((76241, 15015), (170441, 168392), (208923, -365600), (43929, 146360), (-113130, -29596), (-222331, 3841), (-88304, 68617), (296526, 116225)),
    ((-244926, -91150), (204903, -72793), (13590, 45361), (-139589, 1717), (-96047, -12371), (-69591, 291575), (30885, -230766), (-319906, -37291)),
)


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    other = run_path(str(root/"ch39-rephased-involution-up-centers.py"))
    data = other["verify"](root,return_matrices=True)
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    prior = run_path(str(root/"ch39-rephased-four-phase-down-floors.py"))
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    n = rank["N"]+1
    _, labels, _ = train(root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    us,ut,ds,dt,gs = [spec(shell,q,j)[1] for q,j in ((-4,n),(-4,n+1),(2,n),(2,n+1),(6,n))]
    identity = rank["eye"](19)
    sites = shell["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    rot = acb_mat([[int(i==cycle[j]) for j in range(19)] for i in range(19)])
    w = acb(-1,arb(3).sqrt())/2
    s19 = acb_mat([[gs[i,18]] for i in range(19)])
    chars = [(identity+rot*z+rot*rot*z*z)*(arb(1)/3) for z in (w,w.conjugate())]
    high, = [c for c in chars if (adj(s19)*c*s19)[0,0].real > arb(".99")]
    f = data["first_center"]
    g = data["source_differences"][1]
    gc = g.mid()
    free = identity-high-prior["line"](f)
    seed = acb_mat([[acb(a,b) for a,b in row] for row in SEED])
    first = free*acb_mat([[seed[i,j] for j in range(4)] for i in range(19)])
    p1 = other["span"](first)
    second = (free-p1)*acb_mat([[seed[i,j] for j in range(4,8)] for i in range(19)])
    p2 = other["span"](second)
    assert p1.trace().contains(4) and p2.trace().contains(4)
    u = identity+(w-1)*p1+(w.conjugate()-1)*p2
    h0 = identity-2*(prior["line"](f)+prior["line"](gc))
    x = h0*u
    g_radius = prior["frobenius"](g-gc)
    gnorm = (adj(gc)*gc)[0,0].real.sqrt()
    fnorm = (adj(f)*f)[0,0].real.sqrt()
    sine = data["first_motion"]/(fnorm-data["first_motion"])
    ex = data["up_error"] + 2*(sine + g_radius/(gnorm-g_radius))
    def four(target, amps):
        c = acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
        return acb_mat([[c[i,j]*amps[j] for j in range(19)] for i in range(4)])*adj(c)
    up = four(ut,rank["read"](us,ut,labels[n,-4],u))
    down = four(dt,rank["read"](ds,dt,labels[n,2],x))
    values, vectors = up.eig(right=True,algorithm="rump")
    seq = sorted(range(4),key=lambda j:values[j].real.mid())
    lam = [values[j].real for j in seq]
    basis = []
    for j in seq:
        v = acb_mat([[vectors[i,j]] for i in range(4)])
        basis.append(v*(1/(adj(v)*v)[0,0].real.sqrt()))
    eu = data["up_error"]*max(v.upper() for v in labels[n,-4]).sqrt()
    ed = ex*max(v.upper() for v in labels[n,2]).sqrt()
    interval = (arb(".8587"),arb(".8645"))
    def actual_box(central, shift):
        return ((central.sqrt()-shift)**2).union((central.sqrt()+shift)**2)
    for val,window in zip((arb(".00000615"),lam[2],interval[0].union(interval[1])),rank["WINDOWS"],strict=True):
        ball = actual_box(val,eu)
        assert arb(window[0]) < ball < arb(window[1])
    dlevels = sorted((v.real for v in down.eig(algorithm="rump")),key=lambda a:a.mid())
    strange = actual_box(dlevels[1],ed)
    assert strange > arb(".002")
    print("第二个根在任意三平面上的奇质量下界",strange.lower().str(12))
    def bottom(mt):
        # 选取包含第三条上型本征线的三平面。另三条线压缩为两个指定本征值。
        selected = (0,1,3)
        mu = (arb(".00000615"), mt)
        probs = [prod(lam[i]-m for m in mu)/prod(lam[i]-lam[j] for j in selected if j!=i) for i in selected]
        assert min(probs)>0 and sum(probs).contains(1)
        columns = []
        for m in mu:
            v = sum((basis[i]*(p.sqrt()/(lam[i]-m)) for i,p in zip(selected,probs,strict=True)),acb_mat(4,1))
            columns.append(v*(1/(adj(v)*v)[0,0].real.sqrt()))
        frame = acb_mat([[v[i,0] for v in (*columns,basis[2])] for i in range(4)])
        assert (adj(frame)*frame).contains(rank["eye"](3))
        centre = max((v.real for v in (adj(frame)*down*frame).eig(algorithm="rump")),key=lambda v:v.mid())
        return actual_box(centre,ed)
    a,b = [bottom(t) for t in interval]
    assert b < arb(".01392") < arb(".01410") < a
    print("另一立方根下的中心三平面粲质量和两个顶质量参数",lam[2].str(15),[s.str(7) for s in interval])
    print("顶参数低、高端的实际底质量外包",[(v.lower().str(10),v.upper().str(10)) for v in (a,b)])
    return a,b


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
