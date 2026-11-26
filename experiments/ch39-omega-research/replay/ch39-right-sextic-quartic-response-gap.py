"""十个实际响应上的四次有理逼近下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_series, ctx


def verify(root):
    sextic=run_path(str(root/"ch39-common-sextic-full-ckm.py"))
    lib=run_path(str(root/"ch39-sparse-rational-full-ckm.py"))
    number,t,form,seed=sextic["CASES"][0]
    native,model,A,B,mass,nodes,coords=lib["setup"](root,number)
    ctx.prec=320
    spectra=sextic["projectors_and_knots"](nodes,coords)
    root_values=list(map(arb,seed))
    radius=arb("1e-9")
    pars=[v+arb(0,radius.upper()) for v in root_values]
    weak=native["target_basis"](A,arb(t)*B)*model["diag"]([1,1,acb(0,-1)])
    targets=[]
    for mat,ms in zip((A,B),mass,strict=True):
        M=model["adj"](weak)*mat*weak*ms[2]
        targets.append((M[0,0].real,M[1,1].real,M[2,2].real,M[0,2].real,M[0,2].imag,M[1,2].real))
    def equation(x):
        return lib["equations"](sextic["responses"](x,spectra,form),spectra[-1],targets,())
    def jac(x):
        cols=[[z[1] for z in equation([arb_series([v,int(i==j)]) for i,v in enumerate(x)])]
              for j in range(9)]
        return arb_mat(cols).transpose()
    ctx.cap=2
    pre=jac(root_values).inv()
    correction=pre*arb_mat([[v] for v in equation(root_values)])
    ident=arb_mat([[int(i==j) for i in range(9)] for j in range(9)])
    defect=ident-pre*jac(pars)
    gain=max(sum(abs(v).upper() for v in row) for row in defect.tolist())
    assert max(abs(z).upper() for z in correction.entries())+radius*gain<radius
    values=sextic["responses"](pars,spectra,form)[1:39]
    energy=nodes[0]+nodes[1]
    left=(16,28,15,18,10)
    right=(34,7,9,35,37)
    C=arb_mat([[1/(energy[i]-energy[j]) for j in right] for i in left])
    L=arb_mat([[C[k,l]*(values[i]-values[j])
                for l,j in enumerate(right)] for k,i in enumerate(left)])
    assert not L.det().contains(0)
    A=L.inv()
    absolute=lambda M:arb_mat([[abs(z).upper() for z in row] for row in M.tolist()])
    coefficient=absolute(A)*absolute(C)+absolute(A*C)
    bound=1/max(sum(row).upper() for row in coefficient.tolist())
    assert bound>arb("1.93131e-6"), bound
    print("十个实际节点上四次有理响应的逐层极大误差下界",bound.lower().str(12))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
