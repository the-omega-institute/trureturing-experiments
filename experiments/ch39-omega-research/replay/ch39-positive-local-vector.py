"""Realize a mixed spectral witness with a single positive-kinetic massive vector."""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

mixed = run_path(str(Path(__file__).with_name("ch39-mixed-spectral-witness.py")))
poles, mat, q = (mixed[key] for key in ("poles","mat","q"))
site = arb_mat([[1],[q("2/5")],[0]])
family = mixed["family"]
basis = arb_mat([[family[i,0],site[i,0]] for i in range(3)])


def local_witness(logs, model, time):
    scale,_,rho,eta,_,masses = model
    n1 = poles["reference_readings"](logs)[3]
    e = n1**2/scale
    t = (4*masses[0])**2
    h = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
    k = 1/(poles["offset"]+poles["study"]["comp"]["model"]["target"](11430,time)/(15*scale))
    kinetic_inverse = arb_mat([[k/24,0,0],[0,k/4,0],[0,0,k*rho/eta]])+mat["hypercharge"]*(h-k/24)
    assert mat["positive"](kinetic_inverse)
    forms,packets = mixed["reference_problem"](logs,*model,time,basis=basis)
    difference = forms[0]-forms[1]
    a,b,c = difference[0,0],difference[0,1],difference[1,1]
    discriminant = b*b-a*c
    assert discriminant>0 and not a.contains(0)
    candidates = []
    for sign in (-1,1):
        slope = (-b+sign*discriminant.sqrt())/a
        v = arb_mat([[slope],[1]])
        norm2 = mixed["bilinear"](v,forms[0],v)
        assert norm2>0
        residue = basis*v*(1/norm2.sqrt())
        cost = mixed["bilinear"](residue,kinetic_inverse,residue)/(t+e)
        candidates.append((cost,slope,residue))
    cost,slope,residue = min(candidates,key=lambda row:float(row[0].mid()))
    return k, h, e, t, kinetic_inverse, cost, slope, residue, packets


def verify():
    time = poles["study"]["nominal"][2]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        model,_,_ = poles["retune"]["predictions"](logs,time)
        _, _, e, t, kinetic_inverse, cost, slope, residue, packets = local_witness(logs,model,time)
        assert q(".825") < cost < q(".836")
        bare = kinetic_inverse.inv()
        u = residue*(1/t.sqrt())
        matched = bare+u*u.transpose()*(e/(t+e))
        block = arb_mat([[matched[i,j] if i<3 and j<3 else
                          (u[min(i,j),0] if min(i,j)<3 else 1) for j in range(4)] for i in range(4)])
        assert mat["positive"](matched) and block.det()>0
        for factor,inverse,ref,normalization in packets:
            updated = inverse-residue*residue.transpose()*factor
            assert mat["positive"](updated)
            response = updated.inv()[2,2]/normalization
            assert response.overlaps(ref) and abs(response-ref) < q("1e-6")
        print("case",case,"family / first site overlap",slope.str(10),
              "kinetic fraction",cost.str(10),flush=True)


if __name__ == "__main__":
    verify()
