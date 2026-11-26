"""Phase transversality of the one-loop crossing at one fixed shell mass scale."""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, acb_mat

folder = Path(__file__).resolve().parent
flow = run_path(str(folder/"ch39-shell-gw-first-heavy-flow.py"))
full = flow["full"]
shell = full["model"]["m"]
q, adj, hs = arb, shell["adj"], shell["hs"]


def family(sing, G, p, r):
    pts = shell["pts"]
    result = shell["yuk"]["identity"](48)*4
    for i in range(6):
        for j in range(6):
            if i == j or j == i^1:
                continue
            edge = sing[i,j]
            if 26-sum(a*b for a,b in zip(pts[j],pts[i])) != 1:
                for k in range(8):
                    result[8*i+k,8*j+k] = edge
                continue
            U = acb_mat(shell["model"]["link"](pts[j],pts[i],p,r))
            for a, src in enumerate(G):
                image = U*src*adj(U)
                for b, dest in enumerate(G):
                    result[8*i+b,8*j+a] = edge*hs(dest,image)
    return result


def derivative(calc, p, r, phase):
    def offset(t):
        return calc(p+t,r) if phase == 0 else calc(p,r+t)
    return (offset(q.pi()/4)-offset(-q.pi()/4)
            +(1-q(2).sqrt())*(offset(q.pi()/2)-offset(-q.pi()/2))/2)


def check(index):
    sing, original, _, G = shell["transported_spectra"]()
    h = acb.zeta_zero(index).imag
    a = [(1-q(d)/26).acos() for d in (3,5,7)]
    p,r = [h*(x-a[0]) for x in a[1:]]
    width = q("3e-8" if index==43 else "1e-8")
    calc = lambda x,y: family(sing,G,x,y)
    if index == 43:
        for entry in (original-calc(p,r)).entries():
            assert entry.contains(0)
    # Outward-rounded first singlet eigenvalue is independent of both phases.
    lo,hi = q(".01375939"),q(".01375942")
    assert full["count"](full["rows"](sing),lo)==0 and full["count"](full["rows"](sing),hi)==1
    eye = shell["yuk"]["identity"](48)
    x,y = p+arb(0,width),r+arb(0,width)
    K = calc(x,y)
    resolvent = (K+(q(0).union((q(125)/q("1e12"))**2)-lo.union(hi))*eye).inv()
    slopes = []
    for phase in (0,1):
        A = derivative(calc,x,y,phase)
        z = -shell["trace"](resolvent*A)/2
        assert z.imag.contains(0)
        slopes.append(z.real)
    print("window",index,"exposure slopes wrt two phases",
          [(float(s.lower()),float(s.upper())) for s in slopes])
    if index == 43:
        assert slopes[1] > q("32.63") and slopes[1] < q("32.76") and abs(slopes[0]) < q(".007")
    else:
        assert slopes[1] > q("12.816") and slopes[1] < q("12.818") and abs(slopes[0]) < q(".0035")
    # Only 53 shifted masses occur, and the smallest gap exceeds .001.
    assert full["count"](full["rows"](calc(p,r)),hi+q(".001")+4*width)==0
    mass_sensitivity = slopes[1]/(q(53)/(1+(q(125)/q("1e12"))**2/q(".001"))).union(q(53))
    print("d(log common mass)/d(second phase)",
          float(mass_sensitivity.lower()), float(mass_sensitivity.upper()))
    return slopes


if __name__ == "__main__":
    for n in (43,195):
        check(n)
