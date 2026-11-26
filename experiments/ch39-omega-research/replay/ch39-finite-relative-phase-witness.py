"""Construct and replay rational vertex rotors at one percent cover asymmetry."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx
from sympy import Rational

folder = Path(__file__).resolve().parent
local = run_path(str(folder/"ch39-finite-relative-yukawa.py"))
witness = run_path(str(folder/"ch39-finite-flat-flavor-witness.py"))


def ball(q):
    return arb(int(q.p))/int(q.q)


def complex_ball(q):
    re,im = q.as_real_imag()
    return acb(ball(re),ball(im))


def rational_phase(z):
    sign = 1 if z.real>0 else -1
    t = sign*z.imag/(1+sign*z.real)
    denominator = 1<<64
    numerator = int((denominator*t+arb(1)/2).floor().unique_fmpz())
    assert abs(numerator)<=denominator
    return sign*acb(*witness["cs"](arb(numerator)/denominator))


def main():
    ctx.prec = 384
    H,U,D = local["benchmark"]()
    pairs = list(combinations(range(4),2))
    axes = set(range(4))
    comp = lambda p:tuple(sorted(axes-set(p)))
    beta = Rational(1,100)
    def unmix(Y,p):
        q = comp(p)
        return -(Y[p]+Y[q])/2-(Y[p]-Y[q])/(2*beta)
    reservoirs = {}
    for p in pairs:
        for chi in (1,-1):
            m = complex_ball((unmix(U,p)+chi*unmix(D,p).conjugate()).expand())*4*arb(24).sqrt()/5
            norm2 = m.real**2+m.imag**2
            assert norm2>0 and norm2<1
            skew = acb(0,((1-norm2)/norm2).sqrt())
            reservoirs[p,chi] = [rational_phase(m*(1+skew)),rational_phase(m*(1-skew))]
    # Independently rebuild point-local matrices on the common normalized triplet.
    root_forms = {}
    for i,j in pairs:
        diff = H[i,:]-H[j,:]
        form = diff.T*diff
        root_forms[i,j] = acb_mat([[ball(form[a,b]) for b in range(3)] for a in range(3)])
    up,down = acb_mat(3,3),acb_mat(3,3)
    for point in permutations((0,1,3,4)):
        p = tuple(sorted((point.index(0),point.index(1))))
        chi = (-1)**sum(point[i]>point[j] for i,j in pairs)
        z = reservoirs[p,chi].pop()
        f = (root_forms[p]*ball(1+beta)+root_forms[comp(p)]*ball(1-beta))*(5/(32*arb(24).sqrt()))
        up += f*z
        down += f*(chi*z.conjugate())
    assert all(not remains for remains in reservoirs.values())
    for obtained, expected in [(up,H.T*U*H),(down,H.T*D*H)]:
        assert all(abs(obtained[i,j]-complex_ball(expected[i,j].expand())) < arb("1e-18")
                   for i in range(3) for j in range(3))
    target_up = [arb(v) for v in (".00000615",".00311",".8616")]
    target_down = [arb(v) for v in (".0000135",".000268",".01401")]
    P = witness["projectors"](up,target_up)
    S = witness["projectors"](down,target_down)
    ref_up,ref_down = [acb_mat([[complex_ball(Y[i,j].expand()) for j in range(3)] for i in range(3)])
                       for Y in (H.T*U*H,H.T*D*H)]
    P0 = witness["projectors"](ref_up,target_up)
    S0 = witness["projectors"](ref_down,target_down)
    tr = lambda A:sum(A[i,i] for i in range(3))
    err = max(abs(tr(p*s)-tr(p0*s0)).upper() for p,p0 in zip(P,P0) for s,s0 in zip(S,S0))
    jcp = tr(P[0]*S[0]*P[1]*S[1]).imag
    assert err < arb("1e-14") and jcp>arb("0.00003208") and jcp<arb("0.00003210")
    print("24 rational rotors, beta=1/100; mass matrix element errors < 1e-18, CKM squared-modulus error",
          err.str(6),"J",jcp.str(12))


if __name__=="__main__":
    main()
