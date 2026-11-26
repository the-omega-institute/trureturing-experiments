"""A simultaneous coupling crossing of the nonlinear, step-matched two-loop gauge system."""
from pathlib import Path
from runpy import run_path
from math import ceil
import numpy as np
from flint import acb_mat, arb

w = Path(__file__).resolve().parent
low = run_path(str(w/"ch39-shell-gw-coupled-light-flow.py"))
inert = run_path(str(w/"ch39-shell-gw-two-inert-heavy.py"))
flow, coeff, q = low["flow"],low["coeff"],low["q"]
fit = flow["fit"]
adj, hs, rows, count = (fit[k] for k in ("adj","hs","rows","count"))
loop = flow["loop"]


def frob(M):
    return sum(abs(z).abs_upper()**2 for z in M.entries()).sqrt()


def continuous_light(index,data,reference,e=arb(0),basis_move=arb(0),flavor_error=arb(0)):
    F,D = fit["readouts"](index)
    P = acb_mat(rows(F)+rows(D))
    T = adj(P)*(P*adj(P)).inv()
    assert count(rows(P*adj(P)),q(1,400))==0
    N = fit["common"]["old"]["m"]["eye"](54)-T*P
    U,A,B = [fit["orbit"]["exact_frame"](np.array(z).reshape(3,3)) for z in fit["common"]["old"]["units"]]
    y = acb_mat([[z/inert["g"],0,0] for z in
                 (U*reference[0]*A).entries()+(U*reference[1]*B).conjugate().entries()])
    grams = inert["fingerprints"](index)
    z = inert["trial"](P,y,N,grams)
    X = z+T*(y-P*z)
    initial2 = q(0)
    gradients = []
    for j in range(6):
        site = acb_mat(rows(X)[9*j:9*j+9])
        gram = adj(site)*site
        assert count(rows(gram),q(".42")**2)==3
        initial2 += frob(gram-grams[j])**2
        for axis in inert["axes"]():
            local = (2*site*axis).table()
            embed = acb_mat([local[i%9] if i//9==j else [0,0,0] for i in range(54)])
            gradients.append((N*embed).entries())
    V = acb_mat(gradients)
    real = acb_mat([[z.real for z in row] for row in rows(V*adj(V))])
    sigma = q(".115" if index==43 else ".14")
    assert count(rows(real),sigma**2)==0
    change = ((sum(frob(m-r)**2 for m,r in zip(data[:2],reference[:2]))).sqrt()+flavor_error)/inert["g"]
    shift = 20*e/(1-20*e)
    move = 20*change/(1-20*e)+shift*frob(X)
    eta = initial2.sqrt()+q(".84")*move+move**2+2*basis_move+basis_move**2
    # Project the fixed right inverse into the new readout kernel. The map
    # and its fixed point then depend continuously on the common mass scale.
    C = (1+shift)/sigma
    slope = q(".84")*shift/sigma+2*move*C
    assert slope+4*eta*C**2 < 1
    assert slope+2*eta*C**2 < q(".5")
    return sum((hs(M,M).real.sqrt()+3*e)**2 for M in (F,D))



def dot(u,v):
    return sum(x*y for x,y in zip(u,v))


def combine(row,stack):
    return [dot(row,col) for col in zip(*stack)]


def shifted(m,s):
    return (125**2+(m**2-125**2)*(2*s).exp()).sqrt()


def run_face(Q,m,s,tau,params,point,light_trace):
    weights = [5,-3,-2]
    rho = q(164,440)  # (41/10+7+53/10)/44
    norm = [1-5*rho,3*rho,-1+2*rho]
    combos = [weights,norm]
    low_data,z = point
    high = sum(frob(Y)**2 for Y in low_data)
    first = shifted(m[0],s)
    initial = [4*arb.pi()/v-b*(first/1000).log()/(2*arb.pi()) for v,b in zip(flow["start"],flow["beta"])]
    mixed = [dot(c,z) for c in combos]
    Qnew = Q*tau.exp()
    last = shifted(m[-1],s)
    assert Qnew>last
    heavy_bound = q("7.9")**2*(q(81,977)**2)*light_trace
    for k in range(1,54):
        if k==3:
            high += heavy_bound
        if k==53:
            dt = (Qnew/last).log()
        else:
            low_end,hi_end = params
            def ratio(v):
                return (shifted(m[k],v)/shifted(m[k-1],v)).log()
            # Correlated common scale; the ratio increases with s when the masses
            # are ordered. Adjacent intervals with repeated masses may touch.
            dt = ratio(low_end).union(ratio(hi_end))
            dt = arb.max(dt,arb(0))
        steps = max(1,ceil(float(dt.upper())*10))
        h = dt/steps
        slopes = [b+k*d for b,d in zip(flow["beta"],coeff["db"])]
        B = [[a+k*e for a,e in zip(b,ext)] for b,ext in zip(coeff["base"],coeff["extra"])]
        Bmix = [combine(c,B) for c in combos]
        Dmix = [combine(c,coeff["yukawa"]) for c in combos]
        for _ in range(steps):
            assert 10*h*high < loop
            bound = high/(1-10*h*high/loop)
            assert bound < q("3.1")
            span = q(0).union(h)
            inverse = [v-b*span/(2*arb.pi())+d+arb(0,".05") for v,b,d in zip(initial,slopes,z)]
            assert min(inverse)>10
            gbox = [4*arb.pi()/v for v in inverse]
            def rate(b,d):
                # t>=0 and sum(t)<=bound. Linear functional extrema on the simplex.
                y = min([q(0),*d]).union(max([q(0),*d]))*bound
                return (y-dot(b,gbox))/(32*arb.pi()**3)
            dz = [rate(b,d) for b,d in zip(B,coeff["yukawa"])]
            assert all(abs(h*v)<q(".049") for v in dz)
            z = [v+h*a for v,a in zip(z,dz)]
            mixed = [v+h*rate(b,d) for v,b,d in zip(mixed,Bmix,Dmix)]
            initial = [v-h*b/(2*arb.pi()) for v,b in zip(initial,slopes)]
            high = bound
    # mass response written as log(1 + (1-m_light²/m²)*(exp(2s)-1))/2.
    mass_shift = sum((1+(1-125**2/mass**2)*((2*s).exp()-1)).log()/2 for mass in m)
    x = -44*tau/(2*arb.pi())+mixed[0]
    y = mass_shift/(20*arb.pi())+mixed[1]
    inv3 = 4*arb.pi()/flow["start"][2]+7*(Qnew/1000).log()/(2*arb.pi())+z[2]
    return x,y,inv3,high


def check(n):
    Q,m = flow["masses"](n)
    s0,s1,t0,t1 = map(q,(".311",".338","-.114","-.094"))
    s,t = s0.union(s1),t0.union(t1)
    def light(ss):
        return low["light"]((shifted(m[0],ss)/1000).log()/loop)
    point = light(s)
    middle = light((s0+s1)/2)
    light_trace = continuous_light(n,point[0],middle[0])
    faces = [(s,t0),(s,t1),(s0,t),(s1,t)]
    values = [run_face(Q,m,ss,tt,(s0,s1),point,light_trace) for ss,tt in faces]
    for j,(a,b,_,_) in enumerate(values):
        print("zero",n,"face",j,[(float(z.lower()),float(z.upper())) for z in (a,b)])
    assert values[0][0]>0 and values[1][0]<0
    assert values[2][1]<0 and values[3][1]>0
    _,_,inverse,bound = run_face(Q,m,s,t,(s0,s1),point,light_trace)
    print("simultaneous crossing GeV",(Q*t.exp()).str(8),"inverse alpha",inverse.str(8),
          "total Yukawa squared norm max",bound.upper().str(8))


if __name__=="__main__":
    for n in (43,195):
        check(n)
