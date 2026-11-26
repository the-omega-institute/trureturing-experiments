"""Exact CKM fit with conditional phase stability, and a lower-energy joint rearrangement."""
from pathlib import Path
from itertools import permutations
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

joint = run_path(str(Path(__file__).with_name("ch39-shell-gw-joint-reference.py")))
kin, m = joint["kin"],joint["m"]
adj,hs,count,rows = kin["adj"],kin["hs"],kin["count"],kin["rows"]
K = kin["scalar_kinetic"]()
t, lam = arb(".29"),arb(".05")
q = [acb(0,(-1)**i*t/2).exp() for i in range(6)]
arc = kin["arc"]
p = m["c"]["p"]
kp = acb_mat([[K[9*i,9*j] for j in range(6)] for i in range(6)])
pts = m["c"]["m"]["m"]["pts"]
short = {(i,j) for i in range(6) for j in range(6)
         if sum(a*b for a,b in zip(pts[i],pts[j]))==25}


def energy_ref(w):
    v = acb_mat([[z] for z in w])
    return hs(v,kp*v).real + 676*sum(((w[i].conjugate()*kp[i,j]*w[j]).real+arb(25)/26)**2
                                      for i,j in short if i<j)


def tangent(z, q=q, lam=lam):
    v = acb_mat([[z[9*i+j,0] if i==site else 0 for site in range(6)] for i in range(6) for j in range(9)])
    links = adj(v)*K*v
    for i in range(6):
        for j in range(6):
            u = q[i].conjugate()*kp[i,j]*q[j]
            factor = 1+676*(u.real+arb(25)/26) if (i,j) in short else 1
            links[i,j] += lam*u*factor
    for i in range(6):
        links[i,i] -= sum(links[i,j] for j in range(6)).real
    return links


def exact_fit(trial, t=t, lam=lam, target=None):
    q = [acb(0,(-1)**i*t/2).exp() for i in range(6)]
    if target is None:
        target = arc["central_target"]()
    P = acb_mat(rows(m["F"])+rows(arc["conjugate_output"](t)))
    inv = adj(P)*(P*adj(P)).inv()
    z = trial+inv*(target-P*trial)
    form = tangent(z,q,lam)
    res = arb_mat([[sum(form[i,j] for j in range(6)).imag] for i in range(5)])
    kz = K*z
    # For h row i, the real derivative in direction v is Im(h_i v).
    h = acb_mat([[sum(z[s,0].conjugate()*K[s,k] for s in range(9*i,9*i+9))
                  -(kz[k,0].conjugate() if k//9==i else 0) for k in range(54)] for i in range(5)])
    N = adj(h)-inv*(P*adj(h))
    gram = h*N
    C = acb(0,1)*N*acb_mat([[gram[i,j].real for j in range(5)] for i in range(5)]).inv()
    norm = hs(C,C).real.sqrt()
    radius = 2*sum(res[i,0]**2 for i in range(5)).sqrt()
    # Each quadratic phase equation has Hessian norm <= 2||K|| <= 16.
    assert 16*arb(5).sqrt()*norm**2*radius < arb(".001")
    margin = arb(0,(norm*radius).upper())
    return acb_mat([[z[i,0]+acb(margin,margin)] for i in range(54)])


def shuffled(z,perm,w):
    return acb_mat([[w[i]*z[9*j+a,0]/q[j]] for i,j in enumerate(perm) for a in range(9)])


# Floating trial values only seed the exact affine projection and contraction.
seed = [
    (0.041844826733501296, 0.058645226237655293),
    (-0.064021302032267177, -0.0011557839471081017),
    (0.073927507765517333, 0.0038802930916949784),
    (0.48186044294875652, 0.069679640268519458),
    (0.02035542452791814, -0.47274783422964861),
    (0.58489283353293053, -0.34837484746431835),
    (0.29685194934463988, 0.59890764746092295),
    (3.4032077692102112, -0.8357043111875101),
    (0.80899906202582594, 3.4137801031744974),
    (-0.051004074296511305, -0.076254875532483612),
    (-0.068848095094373229, -0.042463604430645292),
    (0.0021571838255308932, -0.043786816643538604),
    (0.34966018844254171, -0.31316429502954457),
    (-0.34092018280505376, -0.33584668296578352),
    (0.66153694880262548, -0.62465461508430487),
    (0.59394740456876727, 0.66326472412674797),
    (3.2876878417186224, -1.6573022721413295),
    (1.6427830503404697, 3.3036259053549504),
    (-0.029240592553516808, 0.077197460665275719),
    (-0.026772194017035088, 0.035342815248291325),
    (0.082358444492782082, -0.023705296992561239),
    (0.57428636132026223, 0.041818013463726732),
    (-0.0069957574149027346, -0.55675028950847316),
    (0.82415148185804232, -0.30782242965299567),
    (0.26424750406831987, 0.79770745067749405),
    (3.1840529241098836, -0.63349330185070774),
    (0.60955390386194019, 3.1811067942844597),
    (-0.081159577406497915, -0.022804449143929491),
    (0.038101877056951616, -0.042223785869325126),
    (-0.068412228480884657, 0.011854513756519847),
    (0.63290216355035811, -0.32701436119744859),
    (-0.33355958208600123, -0.59689668773053417),
    (0.4494258321249282, -0.54336939333413603),
    (0.49611400945198791, 0.46358963128789077),
    (2.5662443215220585, -1.9802810047851154),
    (1.9719323974483214, 2.5781502833036583),
    (0.029211056446812329, 0.091485212435499369),
    (0.010350522496914507, 0.063209981257269079),
    (-0.017603226178853834, 0.053149450266619755),
    (0.73723793539534466, -0.063902689660850998),
    (-0.093431663335908507, -0.72624057866009528),
    (0.59074858271021258, -0.058627996836795045),
    (0.019570374907095133, 0.59749614457087241),
    (3.4049575667576502, -1.4129840460868797),
    (1.3883074154038029, 3.4116180372347071),
    (-0.02502303070805894, -0.058882527129099548),
    (0.012102180115087981, -0.027388326713629323),
    (-0.039401402926007863, 0.058756751632228155),
    (0.52289027366279983, -0.15789360690705626),
    (-0.20267476843660465, -0.50120224324769735),
    (0.34451684468824106, -0.49007587129384367),
    (0.45693127592787169, 0.36319788944860487),
    (2.8466116777644848, -2.0332765877523218),
    (2.0177710181669437, 2.858213172227503),
]
other_phases = [0.11785554241659095, -0.13559209672863962, 0.24713877327526804, -0.28396516363986746, 0.1430579168367648, -0.08849497216011154]


if __name__=="__main__":
    init = acb_mat([[acb(x,y)] for x,y in seed])
    z = exact_fit(init)
    norm = hs(z,z).real.sqrt()
    en = hs(z,K*z).real
    tan = tangent(z)
    anchored = acb_mat([[tan[i,j] for j in range(5)] for i in range(5)])
    assert count(rows(anchored),arb(1))==0
    assert norm < 4*arb.pi() and en/norm**2 < arb("1.347")
    # Fix the six shapes and q, and compare every nontrivial permutation.
    min_gap = arb(".44")
    ident = tuple(range(6))
    for perm in permutations(range(6)):
        if perm != ident:
            moved = shuffled(z,perm,q)
            assert hs(moved,K*moved).real > en+min_gap
    w = [acb(0,v).exp() for v in other_phases]
    moved = shuffled(z,(0,1,5,2,4,3),w)
    diff = hs(moved,K*moved).real-en+lam*(energy_ref(w)-energy_ref(q))
    assert diff < arb("-.68")
    print("strength",norm.str(9),"kinetic quotient",(en/norm**2).str(9),
          "pure phase global minimum; fixed-phase permutation gaps > .44; joint move",diff.str(9))
