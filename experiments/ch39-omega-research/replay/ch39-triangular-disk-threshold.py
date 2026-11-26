"""Fundamental transport eigenvalues and a three-face threshold at the six matched heights."""
from functools import cache
from itertools import combinations_with_replacement, product
from math import comb, factorial, lcm
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

glue = run_path(str(Path(__file__).with_name("ch39-glued-face-stiffness.py")))
trace, q, full = glue["trace"], glue["q"], glue["full"]
D, Q = 55,13
origin = (0,)*5
facts = [factorial(i) for i in range(110)]


@cache
def powers(shape):
    v = trace["vandermonde"](shape)
    return [trace["multiplicity"](shape,v,p,facts) for p in range(101)]


def inserted_tables():
    records = []
    for r in combinations_with_replacement(range(Q+2),4):
        shape = (*reversed(r),0)
        chi = trace["alternant"](shape)//trace["alternant"](origin)
        if chi:
            records.append((shape,chi,trace["vandermonde"](shape)//trace["vandermonde"](origin)))
    common = lcm(*(d for _,_,d in records))
    l,lj,r,rj = [],[],[],[]
    for shape,chi,dim in records:
        a = powers(shape)
        neighbors = [powers(row) for row in glue["fundamental_neighbors"](shape)]
        bj = [sum(z[p] for z in neighbors) for p in range(D+1)]
        factor = chi*(common//dim)
        l.append(a[:D+1])
        lj.append(bj)
        r.append([v*factor for v in a[:Q+1]])
        rj.append([v*factor for v in bj[:Q+1]])
    left,right = fmpz_mat(l).transpose(),fmpz_mat(r)
    return [left*right,fmpz_mat(lj).transpose()*right,left*fmpz_mat(rj)],common


def scalar_coefficients(x,shape):
    # Fourier character coefficient of exp(x*f), f=D-2.
    u = arb(1)
    amount = arb(0)
    for n,mult in enumerate(powers(shape)):
        amount += u*mult
        u *= x/(n+1)
    return (-2*x).exp()*amount


def bound_other(x,Z,delta):
    tail = (64*x).exp()*(64*x)**101/factorial(101)
    assert tail<q("1e-22")
    large = (scalar_coefficients(2*x,origin)+tail).sqrt()/100
    bound = large
    for dynkin in product(range(5),repeat=4):
        shape = tuple(sum(dynkin[i:]) for i in range(5))
        dim = trace["vandermonde"](shape)//trace["vandermonde"](origin)
        if 5<dim<100:
            bound = max(bound,scalar_coefficients(x,shape)/dim+tail)
    return (bound+delta)/Z


def verify():
    ctx.prec = 180
    data,common = inserted_tables()
    high = full["high"]
    power = (high["geo"]["BOOK_PLANCK"]/(1000*q(full["L0"]).exp()))**2/(8*arb.pi())
    t = q(1)/15
    for i in range(11485,11491):
        s = 1/(acb.zeta_zero(i).imag**2+q(1)/4)
        shift = 1-2*(1-t)*s
        alpha = s/shift
        X,Y = 32*power*alpha,32*power*alpha*t
        tail = (X+Y).exp()*((X+Y)**(D+1)/factorial(D+1)+Y**(Q+1)/factorial(Q+1))
        sums = [arb(0) for _ in data]
        u = arb(1)
        for n in range(D+1):
            for j in range(min(n,Q)+1):
                v = u*comb(n,j)*(-t)**j/common
                for h,M in enumerate(data):
                    sums[h] += v*int(M[n-j,j])
            u *= (power-n)*alpha/(n+1)
        den = sums[0]+arb(0,tail)
        Z = shift**power*den
        J,Jz = [(num+arb(0,10*tail))/den for num in sums[1:]]
        color,weak = (J+Jz)/12,(J-Jz)/8
        c = 30+2*t
        trunc = (power*s*c).exp()*(power*s*c)**61/factorial(61)
        x = power*s*(1+t)
        other = bound_other(x,Z,trunc)
        print(i,"color",color.str(10),"weak",weak.str(10),"other <",other.upper().str(9),"Z",Z.str(10),flush=True)
        assert other<min(color,weak) and Z>1 and trunc<q("1e-23")
        # Rest-mass bounds for the cube, including a possible negative part.
        density_max = (power*s*c).exp()/Z
        delta = trunc/Z
        second = (scalar_coefficients(2*x,origin)+q('1e-22'))/Z**2
        assert second<3 and second<density_max
        assert weak>color>other and weak<q('.078') and other<q('.0448')
        r = q("1.538")**3/((1+t)**3*(2*q(5).sqrt()+2*q(10).sqrt())**2)-delta*second
        correction = delta*q(".08580")*second
        H = (r/10-correction)/(1+r)
        bare = 4*power*s/(1+c*s)
        print("triangular disk lower bound",H.lower().str(9),"minimax", (35*H/(2*bare+H)).lower().str(9))
        assert H>q(".00250") and 35*H/(2*bare+H)>q(".0445")
        energy3 = q(".1")*10*(3*color**3+2*weak**3)
        energy2 = q(".1")*10*(3*color**2+2*weak**2)
        H_upper = (energy3+other*(q(".08580")*second-energy2))/(1+r)
        print("Three-face upper bound",H_upper.upper().str(9),"squared density norm <",second.upper().str(9))
        assert H_upper<q(".012")


if __name__ == "__main__":
    verify()
