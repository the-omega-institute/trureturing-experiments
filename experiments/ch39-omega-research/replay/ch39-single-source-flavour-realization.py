"""Real Hermite-free interpolation for one marked-site insertion on a robust flux-three clock."""
from pathlib import Path
import runpy
from flint import acb,arb,arb_mat,acb_mat

src=runpy.run_path(str(Path(__file__).with_name("ch39-single-arc-source.py")))
s=src["shell"]
p,q=(s[key]*src["ref"] for key in ("a","b"))
k,v,g=src["data"](p,q)
e,a,b=src["fourier_columns"](k,v)
strength=arb(".001")
one=acb_mat([[int(i==j) for j in range(3)] for i in range(3)])


def real_filter(target):
    rows=[[z**j for j in range(6)] for z in e]
    vandermonde=arb_mat([[(v.real if part==0 else v.imag) for v in row]
                         for row in rows for part in (0,1)])
    rhs=arb_mat([[part] for v in target for part in (v.real,v.imag)])
    weights=vandermonde.solve(rhs)
    values=[sum(weights[j,0]*z**j for j in range(6)) for z in e]
    assert all((x-y).contains(0) for x,y in zip(values,target))
    return values,max(abs(weights[j,0]).upper() for j in range(6))


def ckm():
    # Explicit benchmark angles in radians, as displayed in the 1 TeV input table.
    su,cu=arb(".22704").sin_cos()
    sb,cb=arb(".04275").sin_cos()
    sx,cx=arb(".00377").sin_cos()
    phase=acb(0,arb("1.139")).exp()
    return acb_mat([[cu*cx,su*cx,sx*phase.conjugate()],
                    [-su*cb-cu*sb*sx*phase,cu*cb-su*sb*sx*phase,sb*cx],
                    [su*sb-cu*cb*sx*phase,-cu*sb-su*cb*sx*phase,cb*cx]])


def adj(x):
    return x.conjugate().transpose()


def realize(masses,left):
    light,middle,heavy=masses
    cosine=(light*heavy+middle**2)/(middle*(light+heavy))
    sine=(1-cosine**2).sqrt()
    rotation=acb_mat([[cosine,0,-sine],[0,1,0],[sine,0,cosine]])
    diag=acb_mat([[masses[i]*int(i==j) for j in range(3)] for i in range(3)])
    target=left*diag*rotation*adj(left)
    rank_one=target-middle*one
    # Factor at the pivot farthest from zero. The 2x2 determinant vanishes analytically.
    _,i,j=max((float(abs(rank_one[i,j]).mid()),i,j) for i in range(3) for j in range(3))
    u=[rank_one[h,j] for h in range(3)]
    row=[rank_one[i,l]/rank_one[i,j] for l in range(3)]
    f,cost_f=real_filter([u[h]/a[h] for h in range(3)])
    h,cost_h=real_filter([-row[l]/(strength*b[l].conjugate()) for l in range(3)])
    rebuilt=acb_mat([[middle*int(i==j)-strength*f[i]*a[i]*b[j].conjugate()*h[j]
                      for j in range(3)] for i in range(3)])
    errors=rebuilt*adj(rebuilt)-left*diag*diag*adj(left)
    assert all(errors[i,j].contains(0) and abs(errors[i,j])<arb("1e-32")
               for i in range(3) for j in range(3))
    return (cost_f*cost_h).sqrt()


if __name__=="__main__":
    # A common real rescaling between the filters balances their maximum coefficients.
    up=realize(list(map(arb,("6.15e-6",".00311",".8616"))),one)
    down=realize(list(map(arb,("1.35e-5",".000268",".01401"))),ckm())
    print("balanced coefficient upper bounds in raw powers of T",up.str(12),down.str(12))
