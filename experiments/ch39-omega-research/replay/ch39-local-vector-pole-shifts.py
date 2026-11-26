"""Enclose the kinetic-mixing pole shifts and test all single-vector directions preserving both poles."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

local = run_path(str(Path(__file__).with_name("ch39-positive-local-vector.py")))
poles, mat, q = (local[key] for key in ("poles", "mat", "q"))
family = local["family"]
bilinear = local["mixed"]["bilinear"]


def pole_direction(energy, b, d):
    return [poles["colour"](energy,mat["unit"])/8,
            -poles["right"](energy,mat["unit"])/4, 2/(b-d*energy)]


def check_preserving_direction(model, k, h_ref, e, e2):
    _,z,rho,eta,c,masses = model
    b,d = c/rho,eta/rho
    a,r = (pole_direction(energy,b,d) for energy in (e,e2))
    cross = [a[(i+1)%3]*r[(i+2)%3]-a[(i+2)%3]*r[(i+1)%3] for i in range(3)]
    assert not cross[2].contains(0)
    w = arb_mat([[v/cross[2]] for v in cross])
    strengths = []
    for log_ratio in (q(".5"),q(1)):
        P,A,_ = local["mixed"]["shape"]["probe"](log_ratio,e,z,b,d,masses[0])
        ref = ((P*(1/rho)+mat["hypercharge"]*h_ref).inv()-A).inv()
        base = ((P*k).inv()-A).inv()
        y = mat["y"]
        v = y+base*A*y
        a_y = bilinear(y,A+A*base*A,y)
        assert 0 < a_y < 1
        # With h = g_Y^2 in [0,1], the updated exchange is
        # base + [h/(1-h*a_y)] v v'. Bound its three relevant projections.
        step = arb(0).union(1/(1-a_y))
        proj_w,proj_j = (v.transpose()*w)[0,0], v[2,0]
        gap = k*rho*ref[2,2]-(base[2,2]+step*proj_j**2)
        mixed_entry = bilinear(w,base,family)+step*proj_w*proj_j
        self_entry = bilinear(w,base,w)+step*proj_w**2
        assert gap>0 and self_entry>0
        strengths.append(1/(self_entry+mixed_entry**2/gap))
    ratio = strengths[0]/strengths[1]
    x1,x2 = q(".5").exp(),q(1).exp()
    maximum = (x1-1)*x2/((x2-1)*x1)
    assert ratio>q(".73500")>q(".62246")>maximum
    print("pole-preserving direction", *(v.str(10) for v in w.entries()),
          "required shift ratio lower bound",ratio.lower().str(10))


def verify():
    time = poles["study"]["nominal"][2]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        model,_,_ = poles["retune"]["predictions"](logs,time)
        scale,z,rho,eta,c,_ = model
        k,h,e,t,_,_,_,residue,_ = local["local_witness"](logs,model,time)
        b,d = c/rho,eta/rho
        e2 = mat["geometry"]["second_neutral_mass"](scale,z,b,d)**2/scale
        def propagation(energy):
            return mat["massive_exchange"](-energy,z,b,d)*k+mat["hypercharge"]*h
        def secular(energy):
            return (t+e)*(t-energy)/(energy+e)+bilinear(residue,propagation(energy),residue)

        for xi in (e,e2):
            direction = pole_direction(xi,b,d)
            assert not sum(a*b for a,b in zip(direction,residue.entries(),strict=True)).contains(0)
        intervals = ((q(".00001"),e*(1-q("1e-6")),e,-q(6),-q("5.9")),
                     (e*(1+q(".000001")),e2*(1-q("1e-8")),e2,-q(124000),-q(121000)))
        for lo,hi,old,min_ppm,max_ppm in intervals:
            assert secular(lo)>0>secular(hi)
            root = poles["roots"](lambda x:-secular(x),lo,hi)
            shift = (root.sqrt()/old.sqrt()-1)*10**6
            assert min_ppm < shift < max_ppm
            assert not bilinear(family,propagation(root),residue).contains(0)
            print("case",case,"mixed pole",root.str(13),"mass shift (ppm)",shift.str(10),flush=True)
        check_preserving_direction(model,k,h,e,e2)


if __name__=="__main__":
    verify()
