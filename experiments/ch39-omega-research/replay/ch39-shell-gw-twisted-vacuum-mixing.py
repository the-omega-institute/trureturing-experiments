"""Test the equal-coupling construction against global vacua of local singular-value potentials."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb

folder = Path(__file__).resolve().parent
native = run_path(str(folder/"ch39-shell-gw-native-arc-reference.py"))
arc = native["arc"]
m = arc["m"]
a = m["m"]
hs, adj = m["hs"],m["adj"]
count, rows = m["count"],m["rows"]


def scalar_basis():
    s,k,_,_ = a["five"]["m"]["transported_spectra"]()
    small = [(v,True) for v in a["band_frame"](s,1,arb(".02"))]
    small += [(v,False) for v in a["band_frame"](k,2,arb(".02"))]
    p = m["c"]["p"]
    columns = []
    for v,is_singlet in small:
        coeff = []
        for i in range(6):
            for j in range(9):
                if is_singlet:
                    x = v[i,0] if j==0 else 0
                else:
                    x = v[8*i+j-1,0] if j else 0
                coeff.append([p[i].conjugate()*x])
        columns.append(acb_mat(coeff))
    return a["orb"]["orthonormalize"](columns)


def sector_basis(readout,scalar,conjugate):
    output = []
    for col in scalar:
        val = readout*col
        if conjugate:
            val = val.conjugate()
        output.append(acb_mat(3,3,val.entries()))
    gram = acb_mat([[hs(x,y) for y in output] for x in output])
    assert count(rows(gram),arb(".38")**2)==0
    return a["orb"]["orthonormalize"](output)


def commutator_bound(up,down):
    u = [x*adj(y) for x in up for y in up]
    d = [x*adj(y) for x in down for y in down]
    cs = [x*y-y*x for x in u for y in d]
    return sum(hs(x,x).real for x in cs).sqrt()


if __name__=="__main__":
    t = native["relative_phase"](arb(26)**2)
    scalar = scalar_basis()
    u = sector_basis(m["F"],scalar,False)
    d = sector_basis(arc["conjugate_output"](t),scalar,True)
    bound = commutator_bound(u,d)
    assert count(rows(m["F"]*adj(m["F"])),arb(".515")**2)==9
    frac = arb(".004136")
    dist = arb(".515")*frac/(arb(".38")*(1-frac**2).sqrt()-arb(".515")*frac)
    total = bound+4*dist
    print("low-mode commutator",bound.str(9),"with arbitrary vacuum tail below",total.upper().str(9))
    assert total < arb(".02524")
    # Use the fixed t=0 image for every positive stiffness of the native phase potential.
    d0 = sector_basis(m["F"].conjugate(),scalar,True)
    X = acb_mat([[v[i,0] for v in scalar] for i in range(54)])
    parity = acb_mat([[(-1)**(i//9)*X[i,j] for j in range(3)] for i in range(54)])
    base = m["F"].conjugate()*X
    alternate = m["F"].conjugate()*parity
    assert count(rows(adj(base)*base),arb(".407")**2)==0
    assert count(rows(adj(alternate)*alternate),arb(".004")**2)==3
    tmax = (arb(25)/26).acos()-arb(".0779")
    assert tmax < arb(".201")
    error = arb(".201")*arb(".004")+arb(".515")*frac
    relative = error/(arb(".407")*tmax.cos()*(1-frac**2).sqrt()-error)
    uniform = commutator_bound(u,d0)+2*dist+2*relative
    assert uniform < arb(".030")
    print("all positive arc stiffnesses: commutator below",uniform.upper().str(9))
