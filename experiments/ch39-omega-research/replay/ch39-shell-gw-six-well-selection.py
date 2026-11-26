"""Certify a homogeneous six-well potential selecting an exact common-field CKM fit."""
from pathlib import Path
from itertools import product
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-soft-reference.py")))
K,kp,short = m["K"],m["kp"],m["short"]
hs,adj,count,rows = (m[x] for x in ("hs","adj","count","rows"))
t,lam = arb(".38"),arb(".08")
q = [acb(0,(-1)**i*t/2).exp() for i in range(6)]
depths = list(map(arb,(".1689","-39.6764","-8.491475","-8.718125","-39.51565","0")))
seed = [
    (0.023488771192787958, 0.03830501610207418),
    (-0.058576646429702145, 0.0046574612194317383),
    (0.060626743915047972, 0.00053353928226335383),
    (0.55621836949098724, 0.16425438898111697),
    (0.11454961739219686, -0.53925934446157764),
    (0.63515290390688806, -0.50781446681046782),
    (0.45075659229399084, 0.64971694951046777),
    (2.1726866385350756, -0.81627962044376279),
    (0.79163730699925539, 2.1758144315946639),
    (-0.042849580082847338, -0.046817393808901479),
    (-0.078786041813694405, -0.04186622534537944),
    (0.050291970684543462, -0.023992616315360756),
    (0.34990530550109328, -0.36557107398344907),
    (-0.38762267461715783, -0.34443064337471935),
    (0.41994166033279934, -0.58633740307073923),
    (0.57271100889796356, 0.40560989100021494),
    (4.000365061147189, -1.7778639016819808),
    (1.7550615514231074, 4.0158913636191311),
    (-0.022583056550046197, 0.060253249560610338),
    (-0.039566020297912993, 0.040249130311011717),
    (0.067336929300841536, -0.016977711862567687),
    (0.6550020100954298, -0.17079322438981137),
    (-0.21836924067566563, -0.62283631957150942),
    (0.72806612319930153, -0.093779614381185594),
    (0.063780017710238029, 0.70632457823201178),
    (2.8883078907848758, -1.0249811129762227),
    (0.99898688567758254, 2.8888828026277826),
    (-0.064811866123357972, -0.018032378990736453),
    (0.017345813606439196, -0.031334706522548397),
    (-0.06892284140389042, 0.014690430623735205),
    (0.45544485804770513, -0.39057636966832032),
    (-0.38757082260768966, -0.43165960794951796),
    (0.69244460754858905, -0.45836862469917888),
    (0.40218091028728375, 0.69558317066972053),
    (2.6567474766075061, -1.5401281542261609),
    (1.5289685832028836, 2.6704060756892471),
    (0.012364691052090008, 0.066784152690428478),
    (0.021124392271100276, 0.024482060094087515),
    (-0.031468521328137468, 0.081534984456082094),
    (0.58309231189799682, 0.048838982190008783),
    (0.04327924089255486, -0.57085533834230573),
    (0.62998252650711706, -0.019217010409300135),
    (-0.013118529190606189, 0.63763386939033972),
    (3.9487783847640086, -1.8778088760390832),
    (1.8537737225028372, 3.9615732361668838),
    (-0.017598655264737199, -0.033539356835268777),
    (0.0046404763151324429, -0.020801657606592531),
    (-0.039526238770247635, 0.045053623287204289),
    (0.65888712550531281, -0.079393418856617395),
    (-0.12999387732530751, -0.63864781715844532),
    (0.33447938898799318, -0.6305696702640704),
    (0.5916198777806384, 0.3492159489391341),
    (2.0433152876808922, -1.1040081511752042),
    (1.0953589462101017, 2.0577744058199401),
]


def weights(z):
    table = {}
    shapes = [acb_mat([[z[9*i+b,0]/q[i]] for b in range(9)]) for i in range(6)]
    differences = [v-shapes[-1] for v in shapes[:-1]]
    affine = acb_mat([[hs(x,y).real for y in differences] for x in differences])
    assert count(rows(affine),arb(".035"))==0
    for i in range(6):
        for j in range(i,6):
            block = acb_mat([[K[9*i+a,9*j+b] for b in range(9)] for a in range(9)])
            for a,left in enumerate(shapes):
                for b,right in enumerate(shapes):
                    table[i,j,a,b] = hs(left,block*right)+lam*kp[i,j]
    return table


def phase_proposals(H):
    w = np.tile(np.array(list(map(complex,q))), (len(H),1))
    base = np.array([[complex(kp[i,j]) for j in range(6)] for i in range(6)])
    mask = np.array([[(i,j) in short for j in range(6)] for i in range(6)])
    ii = np.arange(6)
    for _ in range(32):
        rot = w.conj()[:,:,None]*w[:,None,:]
        ref = base*rot
        signed = float(lam)*676*mask*(ref.real+25/26)
        tangent = H*rot+signed*ref
        grad = 2*tangent.sum(axis=2).imag
        lap = 2*tangent.real
        lap[:,ii,ii] -= lap.sum(axis=2)
        positive = float(lam)*1352*mask*ref.imag**2
        lap -= positive
        lap[:,ii,ii] += positive.sum(axis=2)
        ev,V = np.linalg.eigh(lap[:,:5,:5])
        step = (V @ ((V.transpose(0,2,1)@grad[:,:5,None])/np.maximum(ev,10)[:,:,None]))[:,:,0]
        step /= np.maximum(1,np.linalg.norm(step,axis=1)/.4)[:,None]
        w[:,:5] *= np.exp(-.8j*step)
    ref = base*w.conj()[:,:,None]*w[:,None,:]
    tangent = H*w.conj()[:,:,None]*w[:,None,:]+float(lam)*676*mask*(ref.real+25/26)*ref
    return ref.real,tangent.sum(axis=2).real


def check_batch(batch,table,base_energy):
    edges = [(i,j) for i in range(6) for j in range(i+1,6)]
    data = []
    for seq,depth in batch:
        T = acb_mat(6,6)
        for i in range(6):
            T[i,i] = table[i,i,seq[i],seq[i]].real
        for i,j in edges:
            value = table[i,j,seq[i],seq[j]]
            T[i,j],T[j,i] = value,value.conjugate()
        data.append(T)
    numeric = np.array([[[complex(z) for z in row] for row in rows(T)] for T in data])
    refs,duals = phase_proposals(numeric)
    for (seq,depth),T,r,d in zip(batch,data,refs,duals):
        const = depth
        for i,j in short:
            if i<j:
                trial = arb(float(r[i,j]))
                extra = lam*676*(trial+arb(25)/26)*kp[i,j]
                T[i,j] += extra
                T[j,i] += extra.conjugate()
                const += lam*676*((arb(25)/26)**2-trial**2)
        # Real trial slopes and dual diagonals can be arbitrary. Only interval
        # positivity and the resulting lower energy bound certify the conclusion.
        ds = list(map(lambda x:arb(float(x)),d))
        shift = (const+sum(ds)-base_energy-arb(".19"))/6
        for i in range(6):
            T[i,i] += shift-ds[i]
        assert count(rows(T),arb(0))==0,seq


def certify(z=None, norm_limit="12.11"):
    if z is None:
        z = m["exact_fit"](acb_mat([[acb(x,y)] for x,y in seed]),t,lam)
    norm = hs(z,z).real.sqrt()
    assert norm < arb(norm_limit) < 4*arb.pi()
    tangent = m["tangent"](z,q,lam)
    anchored = acb_mat([[tangent[i,j] for j in range(5)] for i in range(5)])
    assert count(rows(anchored),arb(1))==0
    table = weights(z)
    base = hs(z,K*z).real+lam*m["energy_ref"](q)+sum(depths)
    diag = {(i,a):table[i,i,a,a].real+depths[a] for i in range(6) for a in range(6)}
    edges = [(i,j) for i in range(6) for j in range(i+1,6)]
    lower = {key: -2*abs(z) for key,z in table.items() if key[0]<key[1]}
    batch = []
    direct = 0
    refined = 0
    for seq in product(range(6),repeat=6):
        if seq==tuple(range(6)):
            continue
        floor = sum(diag[i,a] for i,a in enumerate(seq))+sum(lower[i,j,seq[i],seq[j]] for i,j in edges)
        if floor > base+arb(".19"):
            direct += 1
        else:
            batch.append((seq,sum(depths[a] for a in seq)))
            if len(batch)==128:
                check_batch(batch,table,base)
                refined += len(batch)
                batch.clear()
    if batch:
        check_batch(batch,table,base)
        refined += len(batch)
    assert direct+refined==6**6-1
    print("exact selected quark fit: strength",norm.str(9),"baseline",base.str(10),
          "assignment gaps > .19,",direct,"triangle bounds,",refined,"tangent certificates")


if __name__=="__main__":
    certify()
