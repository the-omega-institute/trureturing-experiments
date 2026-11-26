"""独立六次有理谱见证，以两个三次多项式在全实轴上有界地重建完整夸克中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_series, ctx


CASES = (
    (2306888100086, ".0499143", "roots",
     ("1.14356957824", "-2.63488491886", "-.28671110465", "1.54005472517",
      "18.55192705919", "19.10001079832", "18.77148358431", ".002002102376", ".040178443742")),
)


def projectors_and_knots(nodes, coordinates):
    energies = nodes[0] + nodes[1]
    anchors = (9, 28, 16, 35)
    knots = [energies[i] for i in anchors]
    cards = []
    for i, energy in enumerate(energies):
        if i in anchors:
            cards.append([arb(int(j==anchors.index(i))) for j in range(4)])
        else:
            basis = []
            for j, knot in enumerate(knots):
                product = arb(1)
                for k, other in enumerate(knots):
                    if k != j:
                        product *= (energy-other)/(knot-other)
                basis.append(product)
            cards.append(basis)
    columns = [[(col,19*q+i,None) for i,col in enumerate(rows)] for q,rows in enumerate(coordinates)]
    return energies,knots,cards,columns


def responses(variables,spectra,form):
    energy,knots,cards,_ = spectra
    r=variables[:4]
    a,b,c=variables[4:7]
    if form=="roots":
        def poly(x):
            return (x-a)*(x-b)*(x-c)
    else:
        def poly(x):
            z=x-18
            return ((z+a)*z+b)*z+c
    values=[k*poly(e) for k,e in zip(r,knots,strict=True)]
    readings=[]
    for x,card in zip(energy,cards,strict=True):
        P=poly(x)
        Q=sum(k*v for k,v in zip(card,values,strict=True))
        S=P*P
        readings.append(arb("1e-7")+arb(".9999998")*S/(S+Q*Q+arb("1e-18")))
    return [arb("1e-7"),*readings,*variables[-2:]]


def verify(root,case):
    number,t,form,seed=case
    lib=run_path(str(root/"ch39-sparse-rational-full-ckm.py"))
    native,model,A,B,mass,nodes,coords=lib["setup"](root,number)
    ctx.prec=320
    ctx.cap=2
    basis=native["target_basis"](A,arb(t)*B)*model["diag"]([1,1,acb(0,-1)])
    target=[]
    for operator,ms in zip((A,B),mass,strict=True):
        M=model["adj"](basis)*operator*basis*ms[2]
        assert M[0,1].contains(0) and M[1,2].imag.contains(0)
        target.append((M[0,0].real,M[1,1].real,M[2,2].real,M[0,2].real,M[0,2].imag,M[1,2].real))
    spectra=projectors_and_knots(nodes,coords)
    def equation(box):
        return lib["equations"](responses(box,spectra,form),spectra[-1],target,())
    def jacobian(box):
        jets=[]
        for j in range(9):
            vector=[arb_series([z,int(i==j)]) for i,z in enumerate(box)]
            jets.append([h[1] for h in equation(vector)])
        return arb_mat(jets).transpose()
    x=list(map(arb,seed))
    for _ in range(4):
        f=arb_mat([[v] for v in equation(x)])
        dx=jacobian(x).solve(f)
        x=[(a-b).mid() for a,b in zip(x,dx.entries(),strict=True)]
    radius=arb("1e-24")
    box=[z+arb(0,radius.upper()) for z in x]
    pre=jacobian(x).inv()
    error=max(abs(z).upper() for z in (pre*arb_mat([[v] for v in equation(x)])).entries())
    defect=arb_mat([[int(i==j) for j in range(9)] for i in range(9)])-pre*jacobian(box)
    gain=max(sum(abs(z).upper() for z in row) for row in defect.tolist())
    assert error+gain*radius<radius,(error,gain)
    u,v=box[-2:]
    assert u>0 and v>0 and u*v<1
    lib["reproduce"](root,responses(box,spectra,form),spectra[-1],target,(),mass,model)
    print("编号、六次振幅的收缩常数、九个参数",number,gain.str(10),
          *[z.str(15) for z in box],flush=True)


if __name__=="__main__":
    for case in CASES:
        verify(Path(__file__).resolve().parent,case)
