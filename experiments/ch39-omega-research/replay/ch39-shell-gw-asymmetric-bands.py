"""Certify exact CKM orbits with different up- and down-sector truncations."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, arb_mat

folder = Path(__file__).resolve().parent
orb = run_path(str(folder/"ch39-shell-gw-five-exact-orbit.py"))
four, five = orb["four"], orb["five"]
m = four["band"]
hs, adj = orb["hs"], orb["adj"]
eye = m["yuk"]["identity"]
cuts = {2:".015882",3:".02",4:"1.085",5:"1.087",6:"1.09"}


def band_frame(kernel, number, cut):
    rows = m["rows"](kernel)
    assert m["count"](rows, cut) == number
    vals, vecs = np.linalg.eigh(four["float_matrix"](kernel))
    result = []
    for i in range(number):
        e = arb(float(vals[i]))
        v = acb_mat([[acb(float(x.real),float(x.imag))] for x in vecs[:,i]])
        shift = acb(0,arb("1e-11"))
        system = kernel-(e+shift)*eye(kernel.nrows())
        for _ in range(4):
            v = system.solve(v)*(-shift)
        e = hs(v,kernel*v).real/hs(v,v).real
        assert e < cut
        r = kernel*v-e*v
        delta = arb(0,(hs(r,r).real.sqrt()/(cut-e)).upper())
        assert delta.abs_upper() < arb("1e-20")
        result.append(acb_mat([[v[j,0]+acb(delta,delta)] for j in range(kernel.nrows())]))
    return result


def images(k):
    sing, family, phases, gen = five["m"]["transported_spectra"]()
    assert m["count"](m["rows"](sing), arb(cuts[k])) == 1
    # Each added column lies in its own spectral cut, so every prefix truly
    # spans the corresponding smaller band even when roundoff mixes eigenvectors.
    sources = [(v,[eye(3)/arb(3).sqrt()]) for v in band_frame(sing,1,arb(".015882"))]
    for n in range(2,k+1):
        sources.append((band_frame(family,n-1,arb(cuts[n]))[-1],gen))
    left,right = (m["yuk"]["zero_basis"](sgn) for sgn in (1,-1))
    out = []
    for v, basis in sources:
        y = acb_mat(3,3)
        for i,phase in enumerate(phases):
            l = acb_mat([[left[3*i+j,c] for c in range(3)] for j in range(3)])
            r = acb_mat([[right[3*i+j,c] for c in range(3)] for j in range(3)])
            h = sum((v[i*len(basis)+j,0]*g for j,g in enumerate(basis)),acb_mat(3,3))
            y += phase*l.transpose()*h.transpose()*r
        out.append(y)
    return out


def complement(k):
    modes = orb["orthonormalize"](images(k))
    basis = np.stack([four["float_matrix"](b).flatten() for b in modes]).T
    numeric = np.linalg.svd(basis,full_matrices=True)[0][:,k:]
    comp = []
    for col in numeric.T:
        v = acb_mat([[acb(float(z.real),float(z.imag)) for z in row] for row in col.reshape(3,3)])
        for b in modes:
            v -= hs(b,v)*b
        comp.append(v)
    return orb["orthonormalize"](comp)


def certify(ku,kd,units,comps):
    ck,du,dd = orb["reference"]()
    U,A,B = [four["exact_frame"](np.array(v).reshape(3,3)) for v in units]
    x,y = U*du*A, U*ck*dd*B
    def proj(a,b):
        z = [hs(c,a) for c in comps[ku]]+[hs(c,b) for c in comps[kd]]
        return [v.real for v in z]+[v.imag for v in z]
    r = arb_mat([[v] for v in proj(x,y)])
    g = orb["generators"]()
    zero = acb_mat(3,3)
    cols = [proj(t*x,t*y) for t in g]+[proj(x*t,zero) for t in g]+[proj(zero,y*t) for t in g]
    J = arb_mat([list(c) for c in zip(*cols)])
    C = J.transpose()*(J*J.transpose()).inv()
    eta,n = orb["frob"](r),orb["frob"](C)
    factor = 12*eta*n*n
    assert factor < arb(".5")
    print((ku,kd), "residual",eta.upper().str(7),"right inverse",n.upper().str(7),
          "contraction below",factor.upper().str(7))


# Numerical unitary trials are orthonormalized and verified by interval arithmetic.
trials = {
    (3,6): [
        [
            (0.008510283013095149-0.0021169621673882j), (0.0016211033996354253-0.003738552213220222j), (0.05429070730474049-0.9984783462377533j),
            (0.0004554518773579454+0.0057958362320681455j), (-0.9947344795487776-0.10224127373416494j), (0.0038961032500106217-0.001028769878367415j),
            (0.8462948148513639-0.5326109093474113j), (-0.0021644560051523017-0.005361671694296487j), (0.0023007733892785633+0.00848408206456224j),
        ],
        [
            (0.00612252744006263+0.009513598026869546j), (-0.013831294036895678-0.006409649076149233j), (0.31481395643004073+0.94896353494282j),
            (0.003429452903590915-0.0002046714431940188j), (-0.20951772774666486+0.9776795092819855j), (-0.013089181714802401+0.00787205020462843j),
            (-0.9851184326635619-0.17146975447346494j), (-0.0016320881364661633+0.003169681735629608j), (0.010261241727095104+0.004668324272854846j),
        ],
        [
            (4.1158346909902885e-05-0.019348126519748054j), (-0.29774329025248014-0.5741463781230685j), (0.5584643387417502+0.5190839048339781j),
            (0.011342418837209519+0.013737813135138572j), (-0.7026515515140498-0.2966218342962509j), (-0.6455020640130485-0.03613851067353851j),
            (-0.4695872949265501+0.8824942168506889j), (-4.621505744011109e-07+0.0010743235977747785j), (0.008647779295433025+0.024815144780949618j),
        ],
    ],
    (6,3): [
        [
            (-0.006045962589735042-0.00868926610560211j), (-0.03622727592695484-0.015824675011815607j), (-0.868964894498443-0.4931785875926644j),
            (0.22598348003370808-0.04127176765968167j), (0.9573658337146007-0.1703920881307316j), (-0.04034567496998938+0.0041870270721641005j),
            (-0.883169763582791-0.4088119612121418j), (0.20829189383630514+0.09730800501579719j), (0.0025144750036867228-0.00482344279291651j),
        ],
        [
            (0.0011616404620027408+0.0027641072473013964j), (0.4342940137554446-0.8493907991227242j), (0.2685663754084635-0.13336825815899242j),
            (0.020668295275175667+0.027619014861277186j), (-0.29970184057064514+0.0010844265019035388j), (0.7681473897869469+0.5647452788025111j),
            (-0.9966120804889201-0.07460149405081336j), (-0.008446708215459673+0.005501261663664725j), (0.03212829166389728-0.008070949393794949j),
        ],
        [
            (0.0006916274669312554+0.011237850350082332j), (0.0038626805417055657+0.0015436371512910563j), (-0.19211447513215604+0.9812991173434485j),
            (0.0033566690837520214-0.0013077054231634508j), (0.06470619321386689+0.9978893520391567j), (0.003988669346230848-0.001030942615782178j),
            (-0.9999286846451964+0.001697254805203814j), (-0.0010620447584093944+0.003393904970854858j), (0.01091512429815302+0.002821078054797648j),
        ],
    ],
    (4,5): [
        [
            (0.006903666041721197-0.003533334767363023j), (0.005074536213330186-0.0006834879703520463j), (-0.8755100835520052-0.4831104743904691j),
            (0.22502155353819006-0.008994944686770763j), (-0.7160741221988157-0.6606890531743121j), (0.00010395177817462024-0.003493871806555602j),
            (0.6401749891375831-0.7344387174139494j), (0.22463988737352897-0.01544439898647535j), (0.008469772216372013+0.0015518200313917373j),
        ],
        [
            (0.0019646362467850512-0.01019360836869681j), (-0.031131007081626635-0.03799316640447448j), (0.4654037404992122-0.8836735644947346j),
            (-0.003159080411787679+0.0023620103244799554j), (-0.322849142314483-0.9451686399276117j), (-0.036606094646405955+0.032697222467745896j),
            (-0.9895172760480933+0.14398691632707783j), (-0.00039393549273470437+0.0034493366987267625j), (0.010440994665428627+0.001503612496791352j),
        ],
        [
            (-0.00019066501330952378-0.023046375353381058j), (-0.7102877611270249-0.024801147707715823j), (0.03744060958752283+0.7020991585589436j),
            (0.0189621417959469+0.0011700051175518357j), (0.03300609443467751+0.7026895904705776j), (-0.7039884625838297+0.09579796443480408j),
            (-0.8513928275693722+0.5236775244937356j), (0.0002225918609976777-0.003753277615587859j), (-0.02330866920664525+0.018294064703063048j),
        ],
    ],
}


if __name__ == "__main__":
    comps = {k:complement(k) for k in (3,4,5,6)}
    for (ku,kd), seed in trials.items():
        certify(ku,kd,seed,comps)
