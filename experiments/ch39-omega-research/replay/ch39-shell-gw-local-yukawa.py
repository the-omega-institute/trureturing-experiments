"""Interval-compressed local Yukawa vertices on the three zero modes at positive zero 43."""
from pathlib import Path
from runpy import run_path
from flint import acb,acb_mat,arb,ctx


folder=Path(__file__).resolve().parent
model=run_path(str(folder/"ch39-shell-gw-selected-zero.py"))
ctx.prec=480


def adjoint(a):
    return a.conjugate().transpose()


def identity(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])


def zero_basis(sign):
    h=acb.zeta_zero(43).imag
    angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
    _,matrix=model["action_kernel"](sign,h*(angles[1]-angles[0]),h*(angles[2]-angles[0]))
    gap=arb(1)/100
    count=model["model"]["model"]["count"]
    assert count(matrix,-gap)==count(matrix,gap)==18-3*sign
    unit=identity(36)
    assert max(sum(abs(z) for z in row).upper() for row in matrix)<8
    x=acb_mat(matrix)/8
    b=gap/8
    while 1-b>arb("1e-48"):
        x=x*(3*unit-x*x)/2
        b=b*(3-b*b)/2
    # Include the uniform sign error in each entry before solving the null equations.
    err=arb(0,(1-b).upper()/2)
    start=0 if sign==1 else 18
    proj=acb_mat([[(unit[i,j]-sign*x[i,j])/2+acb(err,err)
                   for j in range(start,start+18)] for i in range(start,start+18)])
    rest=acb_mat([[proj[i,j] for j in range(3,18)] for i in range(3,18)])
    cross=acb_mat([[proj[i,j] for j in range(3)] for i in range(3,18)])
    tail=rest.solve(-cross)
    z=acb_mat([[int(i==j) if i<3 else tail[i-3,j] for j in range(3)] for i in range(18)])
    gram=adjoint(z)*z
    tri=acb_mat(3,3)
    for i in range(3):
        for j in range(i+1):
            v=gram[i,j]-sum(tri[i,k]*tri[j,k].conjugate() for k in range(j))
            if i==j:
                assert v.real>0 and v.imag.contains(0)
                tri[i,i]=v.real.sqrt()
            else:
                tri[i,j]=v/tri[j,j]
    return z*adjoint(tri).inv()


def local_matrices():
    left,right=zero_basis(1),zero_basis(-1)
    out=[]
    for i in range(6):
        x,y,z=model["s"]["pts"][i]
        l=acb_mat([[left[3*i+a,b] for b in range(3)] for a in range(3)])
        r=acb_mat([[right[3*i+a,b] for b in range(3)] for a in range(3)])
        chi=acb((arb(3).sqrt()+1)*(arb(26).sqrt()+z)+x+y,y-x)
        out.append(l.transpose()*r*(chi.conjugate()/abs(chi))**2)
    return out


def is_positive(a):
    return all(acb_mat([[a[i,j] for j in range(k)] for i in range(k)]).det().real>0 for k in range(1,4))


if __name__=="__main__":
    ys=local_matrices()
    total=sum(ys[1:],ys[0])
    eye=identity(3)
    assert is_positive(eye-adjoint(total)*total)
    bound=arb(".03782")
    for y in ys:
        acc=(adjoint(total)*y+adjoint(y)*total)/2
        assert is_positive(acc-bound*eye)
        assert is_positive(arb(".34")**2*eye-adjoint(y)*y)
    y,z=ys[0],ys[1]
    a,b=y*adjoint(y),z*adjoint(z)
    cp=(a*b-b*a).det().imag
    assert arb("-9.70e-15")<cp<arb("-9.69e-15")
    print("uniform lightest/heaviest ratio above", (bound/arb(".34")).lower().str(14))
    print("opposite site pair, imaginary determinant of the mass-squared commutator",cp.str(13))
