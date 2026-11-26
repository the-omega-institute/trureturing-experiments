"""An overlap mass window with exactly one zero per charged octahedral line."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


data=run_path(str(Path(__file__).with_name("ch39-shell-gw-face-flux.py")))
count=data["linear"]["count_below"]
ctx.prec=240


def make(center_sign,with_berry):
    pts,faces=data["maps"]()
    flux=[]
    for face in faces:
        ds=[26-sum(a*b for a,b in zip(pts[i],pts[j])) for i,j in data["edges"](face)]
        flux.append(center_sign*(Fraction(1,4) if ds==[7]*3 else Fraction(5,12))-Fraction(1,8))
    k=data["twisted"](data["cochain"](faces,flux))
    if with_berry:
        for i in range(12):
            for j in range(12):
                x,y,z=pts[j%6]
                a,b,c=pts[i%6]
                overlap=acb(26+z*c+x*a+y*b+arb(26).sqrt()*(z+c),x*b-y*a)
                k[i][j]*=overlap/abs(overlap)
    return k


def verify(matrix,index):
    # The Ginsparg-Wilson operator uses sign(K+3 Gamma/4).
    k=[row.copy() for row in matrix]
    for j in range(12):
        k[j][j]+=(arb(3)/4 if j<6 else -arb(3)/4)
    gap=arb(251)/1000
    assert count(k,-gap)==count(k,gap)==6-index
    assert max(sum(abs(v) for v in row).upper() for row in k)<4
    x=acb_mat(k)/4
    unit=acb_mat([[int(i==j) for j in range(12)] for i in range(12)])
    b=gap/4
    while 1-b>arb("1e-12"):
        x=x*(3*unit-x*x)/2
        b=b*(3-b*b)/2
    # A zero on the negative chirality subspace is in ker(P+).
    # A zero on the positive chirality subspace is in ker(P-).
    off,perturb=(6,1) if index>=0 else (0,-1)
    comp=[[(unit[i,j]+perturb*x[i,j])/2 for j in range(off,off+6)] for i in range(off,off+6)]
    lower=arb(1)/500
    assert count(comp,lower)==0 and lower>(1-b)/2
    return gap,lower-(1-b)/2


if __name__=="__main__":
    for berry in (False,True):
        for sign in (1,-1,0):
            gap,comp=verify(make(sign,berry),sign)
            print("Berry",berry,"center charge",sign,"index",sign,"kernel gap >",gap.str(3),
                  "opposite-chirality exclusion margin >",comp.lower().str(7))
