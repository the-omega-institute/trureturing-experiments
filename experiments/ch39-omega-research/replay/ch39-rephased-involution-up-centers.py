"""荷六第一条像在原生目的层重相后，秩二对合的四相位根实现三个上型中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ANGLES = ('-0.381640198', '0.133278742', '-0.213256179', '-0.499017664', '-0.106573142', '0.549802520', '0.245137920')
# 八个高斯整数种子；前后四列各指定一个非平凡三次单位根。
SEED = (
    ((39896, 57472), (-67646, -60840), (-172077, 127105), (4920, 409838), (-146028, -23113), (70024, 187460), (-165180, 225653), (-82149, 21144)),
    ((236547, 251997), (218314, 14734), (204158, 313271), (-190223, 107466), (-222396, 22476), (-41268, 64328), (-80975, 136731), (-207108, -40061)),
    ((-286899, 55582), (48447, 89665), (-261329, 134653), (214988, -15043), (-301621, -37663), (-71816, -153325), (-305580, 59209), (-80531, -174707)),
    ((-159143, -50659), (171041, 182128), (386896, 155660), (285209, 13399), (-95410, -70641), (107048, 178333), (-44173, 197991), (151342, 69533)),
    ((-153268, 175779), (262114, -146066), (132140, -342702), (-198262, -31119), (-161880, -90713), (80631, 150833), (-128633, 107027), (76503, -156871)),
    ((144924, -299528), (68488, -94910), (126102, -179161), (118478, 397134), (-209928, -76204), (131952, -70953), (-126142, -96544), (-27213, 2235)),
    ((153140, 78040), (-59536, -253636), (31342, 47191), (113429, -11321), (-70372, -189360), (-118282, 159432), (70648, -73276), (-48610, 83124)),
    ((63472, -3862), (-37923, -205772), (-30571, 3733), (233809, -96130), (-74254, -241987), (-29425, 43548), (-25233, -61309), (88983, 145189)),
    ((37376, 114744), (-7279, -248422), (-37250, 79080), (100172, -157978), (-117890, -212311), (-174212, 24538), (33071, 15743), (73936, -5003)),
    ((64691, 208633), (110791, 228815), (-265406, -236251), (100080, 12948), (-73791, 136223), (-193180, 287926), (252737, -144866), (47062, 274130)),
    ((-92073, 2689), (132774, 96732), (47194, -44489), (88038, 40279), (-152850, 31711), (-46691, 146630), (-252566, -208665), (-186906, 328240)),
    ((-454687, 104583), (-78661, -151580), (-10318, 58442), (-105677, 186584), (-71665, -119082), (288243, 286225), (238495, -18582), (12466, 22592)),
    ((-2518, 168554), (94168, 220767), (40943, 76960), (180438, 61552), (-209174, -165740), (269644, -38987), (320039, -172870), (49620, -202893)),
    ((109677, -44413), (94603, 93186), (-37353, -104121), (81478, 214557), (-224843, -11487), (-294298, 86359), (-4656, -58988), (26209, -278461)),
    ((167122, -57220), (-361437, -104390), (32809, -1369), (-8385, -285687), (-389888, -217250), (-45020, 142609), (6910, -15943), (28786, 86532)),
    ((136275, -141390), (331622, -3604), (-393676, 138262), (-123686, -84204), (-175629, -34878), (186490, -99634), (193144, 34410), (113045, 307512)),
    ((-202647, -28373), (19926, -78639), (-5204, 62090), (-95939, 158159), (-272888, 138867), (-338981, -234658), (378198, -3619), (108963, 233678)),
    ((49594, 153858), (116759, 61902), (56564, -147525), (-66171, 121737), (-151436, -52236), (-118298, -97939), (-258235, 80869), (445071, 209452)),
    ((-173202, -255021), (330216, -115080), (31349, 23308), (-70900, -123715), (-180826, -23325), (-2698, 139726), (145650, -166345), (-164232, -122327)),
)


def adj(a):
    return a.conjugate().transpose()


def span(a):
    return a*(adj(a)*a).inv()*adj(a)


def verify(root, *, return_matrices=False):
    grav = run_path(str(root/"ch39-gravity-two-householder-roots.py"))
    rk = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    label, target = run_path(str(root/"ch39-up-unrestricted-character-floor.py"))["shell_data"](root)
    spectral = run_path(str(root/"ch39-adjacent-pvm-transport.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    source = spectral["spectrum"](shell,-4,rk["N"]+1)[1]
    gs = spectral["spectrum"](shell,6,rk["N"]+1)[1]
    # 用高荷六字符的正交补消去第二条差向量。循环字符投影是精确的常数矩阵。
    geo = shell["geo"]["sites"]
    perm = [geo.index((y,z,x)) for x,y,z in geo]
    R = acb_mat([[int(i==perm[j]) for j in range(19)] for i in range(19)])
    identity = rk["eye"](19)
    omega = acb(-1,arb(3).sqrt())/2
    s19 = acb_mat([[gs[i,18]] for i in range(19)])
    characters = [(identity+R*a+R*R*a*a)*(arb(1)/3) for a in (omega,omega.conjugate())]
    high, = [p for p in characters if (adj(s19)*p*s19)[0,0].real > arb(".99")]
    branches, f0, motion = grav["gravity"](root,angles=tuple(map(arb,ANGLES)),with_motion=True)
    ds = branches[-1,1]
    projector = grav["line"](f0)
    free = identity-high-projector
    seed = acb_mat([[acb(x,y) for x,y in row] for row in SEED])
    first = free*acb_mat([[seed[i,j] for j in range(4)] for i in range(19)])
    p1 = span(first)
    second = (free-p1)*acb_mat([[seed[i,j] for j in range(4,8)] for i in range(19)])
    p2 = span(second)
    up = identity + (omega-1)*p1 + (omega.conjugate()-1)*p2
    assert (p1*p2).contains(acb_mat(19,19))
    assert p1.trace().contains(4) and p2.trace().contains(4)
    unit = f0*(1/grav["norm"](f0))
    true = ds[0]*(1/grav["norm"](ds[0]))
    align = (adj(true)*unit)[0,0]
    assert abs(align) > 0
    true = true*(align/abs(align))
    k = true*adj(unit)-unit*adj(true)
    rotation = identity+k+k*k*(1/(1+(adj(unit)*true)[0,0].real))
    actual = rotation*up*adj(rotation)
    h = identity-2*(grav["line"](ds[0])+grav["line"](ds[1]))
    x = h*actual
    assert (x*x*x-h).contains(acb_mat(19,19)) and (x*x*actual).contains(identity)
    # 最小酉旋转把中心差直线送至真实差直线，同时在高字符上保持恒等。
    sine = motion/(grav["norm"](f0)-motion)
    error = 2*sine*(2/(1+(1-sine**2).sqrt())).sqrt()
    cu = acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
    amplitudes = rk["read"](source,target,label,up)
    raw = acb_mat([[cu[i,j]*amplitudes[j] for j in range(19)] for i in range(4)])*adj(cu)
    levels = sorted((z.real for z in raw.eig(algorithm="rump")), key=lambda x:x.mid())
    factor = max(w.upper() for w in label).sqrt()
    boxes = [(max(arb(0),a.sqrt()-factor*error)**2).union((a.sqrt()+factor*error)**2)
             for a in levels]
    print("输运偏差上界、四站中心算子谱",error.str(10),[x.str(13) for x in levels])
    print("实际解的升序四站谱区间", [(x.lower().str(11),x.upper().str(11)) for x in boxes])
    trivial = (identity+R+R*R)*(arb(1)/3)
    orbits = run_path(str(root/"ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    e = acb_mat([[1/arb(len(o)).sqrt() if i in o else 0 for o in orbits] for i in range(19)])
    escape = (identity-trivial)*up*e
    eig = (adj(escape)*escape).eig(algorithm="rump")
    transfer = (max(v.real.lower() for v in eig).sqrt()-error)**2
    assert transfer > arb(".8267")
    print("上型输运离开七维循环不变空间的最大概率下界", transfer.lower().str(12))
    for j,(lo,hi) in enumerate(rk["WINDOWS"]):
        assert boxes[j] < arb(lo) < arb(hi) < boxes[j+1]
    if return_matrices:
        return {"levels":boxes, "gravity":h, "source_differences":ds, "first_center":f0,
                "first_motion":motion, "up_center":up, "up_error":error, "up":actual, "root":x}
    return boxes


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
