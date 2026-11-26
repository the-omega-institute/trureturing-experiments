"""Construct a gapped positive matrix spectral measure matching both Euclidean data."""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

shape = run_path(str(Path(__file__).with_name("ch39-heavy-family-shape.py")))
poles, mat, q = (shape[key] for key in ("poles","mat","q"))
family = arb_mat([[0],[0],[1]])
directions = arb_mat([[mat["y"][i,0], family[i,0]] for i in range(3)])


def bilinear(a, matrix, b):
    return (a.transpose()*matrix*b)[0,0]


def reference_problem(logs, scale, z, model_rho, eta, c, masses, time, *, basis=directions):
    n1 = poles["reference_readings"](logs)[3]
    e = n1**2/scale
    b,d = c/model_rho,eta/model_rho
    gy = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()
    k0 = 1/model_rho
    k1 = 1/(poles["offset"]+poles["study"]["comp"]["model"]["target"](11430,time)/(15*scale))
    bare = mat["massive_exchange"](e,z,b,d)[2,2]
    spectral_mass = 4*masses[0]
    forms, packets = [], []
    for s in (q(".5"),q(1)):
        P,A,x = shape["probe"](s,e,z,b,d,masses[0])
        H = mat["hypercharge"]*gy**2
        direct = ((P*k0+H).inv()-A).inv()
        inverse = (P*k1+H).inv()-A
        alternate = inverse.inv()
        desired = direct[2,2]*k1/k0
        gap = desired-alternate[2,2]
        assert gap>0 and mat["positive"](inverse)
        factor = e*(x-1)/((spectral_mass**2+e)*(spectral_mass**2+x*e))
        v = alternate*family
        form = (alternate+v*v.transpose()*(1/gap))*factor
        forms.append(basis.transpose()*form*basis)
        packets.append((factor,inverse,direct[2,2]/(x*k0*bare),x*k1*bare))
    return forms, packets


def verify():
    time = poles["study"]["nominal"][2]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        model,_,_ = poles["retune"]["predictions"](logs,time)
        forms,packets = reference_problem(logs,*model,time)
        delta = forms[0]-forms[1]
        a,b,c = delta[0,0],delta[0,1],delta[1,1]
        discriminant = b*b-a*c
        assert discriminant > 0
        assert not a.contains(0)
        slope = (-b-discriminant.sqrt())/a
        vector = arb_mat([[slope],[1]])
        norm2 = bilinear(vector,forms[0],vector)
        assert norm2 > 0
        residue = directions*vector*(1/norm2.sqrt())
        assert not slope.contains(0)
        for factor,inverse,ref,normalization in packets:
            updated = inverse-residue*residue.transpose()*factor
            assert mat["positive"](updated)
            answer = updated.inv()[2,2]/normalization
            assert answer.overlaps(ref) and abs(answer-ref) < q("1e-7")
        print("case",case,"Y/X overlap slope",slope.str(12),"residue",
              *(v.str(10) for v in residue.entries()),flush=True)


if __name__=="__main__":
    verify()
