"""A simple branch point gives the full symmetric group of the threshold mass law."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, acb_mat, arb, ctx

mono = run_path(str(Path(__file__).with_name("ch39-shell-gw-mass-monodromy.py")))
q, shell, full = mono["q"], mono["shell"], mono["full"]
ks = (0,1,2,-2,-1)


def symbol(index):
    single, _, _, G = shell["transported_spectra"]()
    phase = acb.zeta_zero(index).imag*((1-q(5)/26).acos()-(1-q(3)/26).acos())
    width = q("3e-8" if index == 43 else "1e-8")
    phase += arb(0,width)
    vals = np.linalg.eigvalsh(full["model"]["array"](single))
    singles = []
    for i,x in enumerate(vals):
        v = q(float(x))
        lo, hi = v-q("1e-13"), v+q("1e-13")
        assert full["count"](full["rows"](single),lo)==i and full["count"](full["rows"](single),hi)==i+1
        singles.append(lo.union(hi))
    eye = shell["yuk"]["identity"](48)
    a = 2*q.pi()/5
    samples = [mono["fixed"]["family"](single,G,phase,j*a) for j in range(5)]
    coeff = [sum((acb(0,-k*j*a).exp()*m for j,m in enumerate(samples)),acb_mat(48,48))/5
             for k in ks]
    coeff[0] -= singles[0]*eye
    return coeff,[v-singles[0] for v in singles[1:]],eye


def box(radius):
    return acb(arb(0,radius.upper()),arb(0,radius.upper()))


def setup(index, center_point=None, data=None):
    coeff, ds, eye = symbol(index) if data is None else data
    c = (acb("-4.20121883243","-1.15767341874") if index==43
         else acb("-4.20122876623","-1.15762007020"))
    if center_point is not None:
        c = center_point
    center = sum((m*c**k for m,k in zip(coeff,ks)),acb_mat(48,48))
    # Scale the nonsingular 47 directions to a near-identity block. Orthogonal
    # trial bases only precondition the calculation; their determinants are checked.
    U, singular, Vh = np.linalg.svd(full["model"]["array"](center))
    U = U.conj().T[::-1].copy()
    U[1:] /= singular[-2::-1,None]
    V = Vh.conj().T[:,::-1]
    def exact(M):
        return acb_mat([[acb(float(z.real),float(z.imag)) for z in row] for row in M])
    left,right = exact(U),exact(V)
    transform = left.det()*right.det()
    assert not transform.contains(0)
    ss = left*right
    mats = [left*m*right for m in coeff]
    anchor = left*center*right
    diag = q(1)
    for d in ds: diag *= d
    def values(s,z,deriv=False):
        # Center the Laurent series by an exact divided difference.
        scale = z-c
        pieces = {0:acb(0),1:acb(1),2:z+c,-1:-1/(z*c),-2:-(z+c)/(z*c)**2}
        K = anchor+scale*sum((m*pieces[k] for m,k in zip(mats,ks)),acb_mat(48,48))+s*ss
        rows = full["rows"](K)
        C = acb_mat([r[1:] for r in rows[1:]])
        inv = C.inv()
        x = acb_mat([rows[0][1:]])*inv
        y = inv*acb_mat([[r[0]] for r in rows[1:]])
        schur = K[0,0]-(x*acb_mat([[r[0]] for r in rows[1:]]))[0,0]
        prod = q(1)
        for d in ds: prod *= d+s
        factor = C.det()*prod/diag
        if not deriv:
            return factor*schur
        def rate(J,extra):
            L = acb_mat([[1,*[-v for v in x.entries()]]])
            R = acb_mat([[1]]+[[-v] for v in y.entries()])
            block = acb_mat([r[1:] for r in full["rows"](J)[1:]])
            return factor*((L*J*R)[0,0]+schur*(shell["trace"](inv*block)+extra))
        dz = sum((m*k*z**(k-1) for m,k in zip(mats,ks)),acb_mat(48,48))
        return rate(dz,0),rate(ss,sum(1/(d+s) for d in ds))
    return c,values,diag/transform


def certify(index, center_point=None, radius=None, data=None):
    c,H,normal = setup(index, center_point, data)
    # The determinant is normalized by a fixed nonzero scale, including the
    # phase-independent positive singlet product.
    zero = acb(0)
    test_radius = q(".00002") if radius is None else radius
    a,_ = H(zero,c,True)
    a0 = a.mid()
    deriv,_ = H(zero,c+box(test_radius),True)
    residual = abs(H(zero,c)/a0).abs_upper()
    lip = abs(1-deriv/a0).abs_upper()
    print("window",index,"root residual",residual.str(6),"derivative variation",lip.str(6),flush=True)
    assert residual+test_radius*lip < test_radius and lip < 1
    A,B = H(zero,c+box(test_radius),True)
    assert abs(B).abs_lower()>0 and abs(A).abs_lower()>0
    R = q(".0001")
    majorant = abs(H(box(R),c+box(R+test_radius))).abs_upper()
    # Rates A,B are derivatives of the normalized H at its exact local zero.
    bscale = (abs(B*normal).abs_upper()/q(2484).exp())**(q(1)/52)
    ratio = abs(B/A).abs_upper()
    # After v=(B/A)*b*t, s=b*u, the remainder divided by B*b
    # has total degree >=2. Bound it on |u|,|t|<=3.
    l = 3*(1+ratio)/R
    rem = majorant/abs(B).abs_lower()*bscale*l**2/(1-bscale*l)
    print("ratio",ratio.str(6),"normalization length",bscale.str(7),
          "scaled analytic remainder",rem.str(7),flush=True)
    assert bscale*l < q(".001") and rem < q("1e-12")
    # Uniform separation of the 52 simple critical points and their images.
    critical = (q(1)/53)**(q(1)/52)
    delta = q(".0001")
    drift = (1+delta/critical)**52-1
    assert 2*52*delta/critical-drift > q(".005")
    assert rem+delta*drift < q("1e-6")
    assert 2*critical*52/53*(q.pi()/52).sin() > q(".10")
    assert 2**53 > 3+rem


if __name__=="__main__":
    for n in (43,195):
        certify(n)
