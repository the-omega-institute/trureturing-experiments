"""Flavor readout on the first Wilson-conditioned replacement zero."""
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx

folder=Path(__file__).resolve().parent
m=run_path(str(folder/"ch39-shell-gw-higgs-low-band.py"))
yuk=m["yuk"]; selected=yuk["model"]
adj,hs,eye=m["adj"],m["hs"],yuk["identity"]
pts=m["pts"]; count=m["model"]["model"]["count"]
ctx.prec=512
height=acb.zeta_zero(195).imag
angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
p,q=[height*(z-angles[0]) for z in angles[1:]]


def zero_basis(sign):
    _,kernel=selected["action_kernel"](sign,p,q)
    gap=arb(".005")
    assert count(kernel,-gap)==count(kernel,gap)==18-3*sign
    assert max(sum(abs(z) for z in row).upper() for row in kernel)<8
    x=acb_mat(kernel)/8;b=gap/8;unit=eye(36)
    while 1-b>arb("1e-48"):
        x=x*(3*unit-x*x)/2
        b=b*(3-b*b)/2
    err=arb(0,(1-b).upper()/2)
    start=0 if sign==1 else 18
    proj=acb_mat([[(unit[i,j]-sign*x[i,j])/2+acb(err,err)
                   for j in range(start,start+18)] for i in range(start,start+18)])
    rest=acb_mat([[proj[i,j] for j in range(3,18)] for i in range(3,18)])
    cross=acb_mat([[proj[i,j] for j in range(3)] for i in range(3,18)])
    tail=rest.solve(-cross)
    cols=[]
    for j in range(3):
        v=acb_mat([[int(i==j) if i<3 else tail[i-3,j]] for i in range(18)])
        for u in cols:
            v-=hs(u,v)*u
        cols.append(v/hs(v,v).real.sqrt())
    return acb_mat([[col[i,0] for col in cols] for i in range(18)])


def scalar_models():
    sing,_,phases,G=m["transported_spectra"]()
    family=eye(48)*4
    for i,j in product(range(6),repeat=2):
        if j==i or j==(i^1):continue
        U=acb_mat(m["model"]["link"](pts[j],pts[i],p,q))
        for a,h in enumerate(G):
            transported=U*h*adj(U)
            for b,z in enumerate(G):
                family[8*i+b,8*j+a]=sing[i,j]*hs(z,transported)
    l,r=zero_basis(1),zero_basis(-1)
    unit=eye(3)/arb(3).sqrt();cols=[]
    for i,phase in enumerate(phases):
        L=acb_mat([[l[3*i+j,k] for k in range(3)] for j in range(3)])
        R=acb_mat([[r[3*i+j,k] for k in range(3)] for j in range(3)])
        cols.extend([(phase*L.transpose()*v.transpose()*R).entries() for v in [unit]+G])
    return sing,family,acb_mat(cols).transpose()


def array(v):
    return np.array([[complex(v[i,j]) for j in range(v.ncols())] for i in range(v.nrows())])


bands=run_path(str(folder/"ch39-shell-gw-asymmetric-bands.py"))
ctx.prec=512
cuts={2:".02719",3:".04",4:"1.068",5:"1.070",6:"1.08",9:"1.116"}


def sources(sing,fam,k):
    sf=bands["band_frame"](sing,1,arb(cuts[2]))
    assert count(bands["m"]["rows"](sing),arb(cuts[k]))==1
    count_family=k-1
    if k==9:
        af=bands["band_frame"](fam,count_family,arb(cuts[k]))
    else:
        af=[bands["band_frame"](fam,j-1,arb(cuts[j]))[-1] for j in range(2,k+1)]
    cols=[]
    for v,is_singlet in [(sf[0],True)]+[(f,False) for f in af]:
        col=[]
        for i in range(6):
            for j in range(9):
                col.append([v[i,0] if is_singlet and j==0 else
                            v[8*i+j-1,0] if not is_singlet and j else acb(0)])
        cols.append(acb_mat(col))
    return cols


def image_basis(F,cols):
    imgs=[acb_mat(3,3,(F*v).entries()) for v in cols]
    return bands["orb"]["orthonormalize"](imgs)


def complement(modes):
    M=np.stack([array(b).flatten() for b in modes]).T
    comp=[]
    for col in np.linalg.svd(M,full_matrices=True)[0][:,len(modes):].T:
        v=acb_mat([[acb(float(z.real),float(z.imag)) for z in row] for row in col.reshape(3,3)])
        for b in modes:
            v-=hs(b,v)*b
        comp.append(v)
    return bands["orb"]["orthonormalize"](comp)


def four_exclusion(bs, frames=None):
    four=bands["four"]
    if frames is None:
        L=four["unitary_trial"](np.linalg.eig(array(bs[1])@np.linalg.inv(array(bs[0])))[1])
        R=four["unitary_trial"](array(bs[0]).conj().T@L)
        frames=four["exact_frame"](L),four["exact_frame"](R)
    L,R=frames
    B=[adj(L)*b*R for b in bs]
    forward=[(0,1),(1,2),(2,0)]
    if abs(B[3][0,1])<arb(".2"):
        forward=[(j,i) for i,j in forward]
    for i,j in forward:
        assert arb(".574")<abs(B[3][i,j])<arb(".576") and abs(B[3][j,i])<arb(".051")
        for a,b in ((i,j),(j,i)):
            assert sum(abs(z[a,b])**2 for z in B[:3]).sqrt()<arb(".00130")
    d=arb(".00131")
    assert (1+2*arb(".01914")**2).sqrt()*arb(".00130")<d
    radii=[]
    for r in map(arb,(".00362",".01914")):
        T=2*(r+d)/arb(".523")
        for _ in range(3):
            off=arb(".627")*T+2*d
            T=((r+(arb(".051")*T+d)**2)/(1-off)+d)/arb(".574")
        off=arb(".627")*T+2*d
        col=((arb(".576")*T+d)**2+(arb(".051")*T+d)**2).sqrt()
        assert r+2*off<1
        radii.append(col/(1-r-2*off))
    assert radii[1]<arb(".02428") and sum(radii)<arb(".032")
    return sum(radii)


# Unit points on the physical orbit; each is certified by exact projection.
trials={
    (3,6): [
        [
            (0.021788648872271034-0.005373019197458424j), (0.007972342469834833-0.0066531810225726245j), (-0.04210816847922959-0.998807020665242j),
            (0.006726153026457017+0.010609999306504405j), (-0.8534068857295298-0.5209959607853614j), (0.009421702436889409-0.0036489091702297046j),
            (0.8392073834368331-0.5432030373721267j), (0.006325033577822337-0.010586176024669968j), (0.008263103064144389+0.02100180301689052j),
        ],
        [
            (-0.015622585435418293+0.02479445185658445j), (-0.08586687288025277-0.04665728456475975j), (-0.7408671644625303+0.6638576598499868j),
            (0.001028991369043199-0.006201466327965763j), (0.7950392435064352+0.5985687406144731j), (-0.08169163286428424+0.053990963651619546j),
            (-0.8934046838165581-0.44825185321226174j), (-0.0061932856132482905+0.006299618721759797j), (0.021307804363164908+0.019138508742954394j),
        ],
        [
            (0.028097785094262306+0.023365165805281606j), (0.043988686930430426-0.07653948066908896j), (0.9710314592314835-0.21901870240577737j),
            (-0.006784760215697101+0.013039177081187597j), (-0.7675961478180906-0.634716870982451j), (0.004411147967658265+0.08772192547221642j),
            (-0.9277964864749971+0.3710016308241781j), (0.0007738003931697578+0.011694105313313066j), (0.010673285240178493-0.03605828235655791j),
        ],
    ],
    (6,3): [
        [
            (-0.018930312887731428+0.0065654884606872405j), (0.05365079144762015+0.00818439275393575j), (0.9983074704211286+0.005944779971446942j),
            (0.13738808384900386-0.1931859724576381j), (0.579224564113417-0.7784257418999261j), (-0.021131289252556834+0.04369704511280927j),
            (-0.8997017965501849+0.3659766513884636j), (0.22011408321247375-0.08459708824183232j), (-0.03064599067327327+0.007191259350700412j),
        ],
        [
            (0.0009116157787885966-0.0014011226597127671j), (0.31414281279826084-0.9493227719879624j), (0.001963355074677624+0.009691171297067277j),
            (-0.0034189914633074876-0.048001831196456686j), (0.008548956650489948+0.005066691937125035j), (0.9987393861957854+0.010247786757886717j),
            (-0.12476009961462702-0.9910177889701424j), (0.00020008144631155218+0.0013310494773790492j), (-0.04808781653783022+0.0020999200467147956j),
        ],
        [
            (-0.0006326831517622725+0.028523845680014422j), (0.018799110906503845-0.020327372251601226j), (-0.2771988984312189+0.9599896624519693j),
            (0.008866824751282913+0.0024775919399205953j), (-0.4846742125196287+0.8742115105102649j), (-0.001361154074420813+0.027541438350890235j),
            (-0.9804543065262497+0.19444943741443216j), (-0.0010319863052263062+0.008801646469629543j), (0.028567588812848627+0.0020317471417453815j),
        ],
    ],
    (4,5): [
        [
            (0.01362909810315912+0.014325500068983048j), (-0.0022066711173400724+0.015434463054664397j), (0.295940177763721-0.9548745189136855j),
            (0.18578695187350128+0.26956769751557247j), (0.8197078825570486-0.4698679079728839j), (0.010920244591390892+0.0005912349757677979j),
            (0.895325027029203+0.30136122623445327j), (-0.12991938997154528+0.3003015091424643j), (0.000698375611040646+0.022671023031866157j),
        ],
        [
            (0.00015010053994761208-0.028868077493295796j), (-0.0992571443867921+0.0948183255022613j), (0.21555659282614345-0.966364046981252j),
            (-0.010951448865056673+0.0030865994198481478j), (-0.26275842251065457-0.955006074274657j), (-0.10349354018541208-0.08989457886818852j),
            (0.8302793465490167+0.5564830248422574j), (0.0070348429968051926-0.005424474108888864j), (-0.019811248331863967-0.02216867746932838j),
        ],
        [
            (0.007612107921747882+0.01575580005008833j), (-0.14644162686101356+0.6406705486699127j), (0.3644762415368711-0.6595051007552689j),
            (0.009854578028583498-0.0028974864526200715j), (-0.594833132525708-0.46287033138489153j), (-0.3867520537439196-0.5312644742143873j),
            (-0.8750587644924501-0.4835912128304738j), (0.003468341423462395+0.0029443794517712905j), (-0.00037469385691003577-0.019770139099742377j),
        ],
    ],
}


if __name__=="__main__":
    sing,fam,F=scalar_models()
    cols=sources(sing,fam,6)
    B=image_basis(F,cols)
    radius=four_exclusion(B)
    comp={k:complement(B[:k]) for k in (3,4,5,6)}
    for pair,seed in trials.items():
        bands["certify"](*pair,seed,comp)
    for lo,hi,j in (("1.0698814590","1.0698814594",4),("1.1158923785","1.1158923790",8)):
        r=bands["m"]["rows"](fam)
        assert count(r,arb(lo))==j-1 and count(r,arb(hi))==j
    full=sources(sing,fam,9)
    X=acb_mat([[v[i,0] for v in full] for i in range(54)])
    Y=F*X
    assert count(bands["m"]["rows"](adj(Y)*Y-arb(".0285")**2*adj(X)*X),arb(0))==0
    print("zero 195 four-mode heavy-line upper sum",radius.str(9),"nine-mode gain > .0285")
