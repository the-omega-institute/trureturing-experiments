"""Conjugacy-averaged boundary stiffness of two glued amplified SU(5) faces."""
from functools import cache
from itertools import combinations_with_replacement
from math import comb, factorial, lcm
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

trace = run_path(str(Path(__file__).with_name("ch39-amplified-su5-trace.py")))
full, q = trace["full"], trace["q"]
DEGREE, WEAK = 75, 18


def fundamental_neighbors(shape):
    # Invert the fundamental and conjugate-fundamental Pieri steps on SU(5).
    for sign in (-1,1):
        for i in range(5):
            row = list(shape)
            row[i] += sign
            if all(a>=b for a,b in zip(row,row[1:])):
                yield tuple(j-row[-1] for j in row)


def tables():
    facts = [factorial(i) for i in range(DEGREE+8)]
    origin = (0,)*5
    vd0, alt0 = trace["vandermonde"](origin),trace["alternant"](origin)
    @cache
    def powers(shape):
        vd = trace["vandermonde"](shape)
        return [trace["multiplicity"](shape,vd,j,facts) for j in range(DEGREE+2)]

    records = []
    for r in combinations_with_replacement(range(WEAK+2),4):
        shape = (*reversed(r),0)
        alt = trace["alternant"](shape)
        assert alt%alt0==0
        if alt:
            records.append((shape,alt//alt0,trace["vandermonde"](shape)//vd0))
    common = lcm(*(dim for _,_,dim in records))
    plain, inserted, zplain, zinserted = [],[],[],[]
    for shape,chi,dim in records:
        weights = powers(shape)
        near = [powers(other) for other in fundamental_neighbors(shape)]
        # 5*C = 18*(D-2)-6*(chi_5+chi_bar5), where C=-L D.
        cw = [18*(weights[p+1]-2*weights[p])-6*sum(row[p] for row in near)
              for p in range(DEGREE+1)]
        scalar = chi*(common//dim)
        plain.append(weights[:-1])
        inserted.append(cw)
        zplain.append([v*scalar for v in weights[:WEAK+1]])
        zinserted.append([v*scalar for v in cw[:WEAK+1]])
    left = fmpz_mat(plain).transpose()
    right = fmpz_mat(zplain)
    M = left*right
    # 75 * (C(U)-C(zU)/15) insertion.
    C = 15*(fmpz_mat(inserted).transpose()*right)-left*fmpz_mat(zinserted)
    assert M[0,0]==common and C[0,0]==0
    assert all(M[i,j]==M[j,i] for i in range(WEAK+1) for j in range(i))
    return M,C,common


def stiffness(index, power, data):
    M,C,common = data
    strength = 1/(acb.zeta_zero(index).imag**2+q(1)/4)
    t = q(1)/15
    offset = 1-2*(1-t)*strength
    alpha = strength/offset
    X,Y = 64*power*alpha,64*power*alpha*t
    err = (X+Y).exp()*((X+Y)**(DEGREE+1)/factorial(DEGREE+1)+Y**(WEAK+1)/factorial(WEAK+1))
    denominator,numerator = arb(0),arb(0)
    coeff = arb(1)
    for n in range(DEGREE+1):
        for j in range(min(n,WEAK)+1):
            w = coeff*comb(n,j)*(-t)**j
            denominator += w*int(M[n-j,j])/common
            numerator += w*(1-n/(2*power))*int(C[n-j,j])/(75*common)
        coeff *= (2*power-n)*alpha/(n+1)
    # |C(U)-tC(zU)| <= (1+t)*(2*5*12/5+2*10*18/5).
    numerator += arb(0,96*(1+t)*err)
    denominator += arb(0,err)
    assert denominator>1
    result = power**2*alpha/(24*(2*power-1))*numerator/denominator
    bare = 4*power*strength/(1+(30+2*t)*strength)
    error = 35*result/(2*bare+result)
    print(index,"isotropic two-face stiffness",result.str(10),"bare color",bare.str(9),
          "common-normalization minimax error at target 35",error.str(9),flush=True)
    return result, error


def verify():
    ctx.prec = 180
    data = tables()
    high = full["high"]
    power = (high["geo"]["BOOK_PLANCK"]/(1000*q(full["L0"]).exp()))**2/(8*arb.pi())
    values = [stiffness(index,power,data) for index in range(11485,11491)]
    centers = (85796954,85762708,85710456,85687985,85640540,85618137)
    for (value,error),center in zip(values,centers):
        assert abs(value-q(center)/10**9)<q("0.000000002") and error>q("1.4608")
    print("Stiffness range",values[-1][0].str(10),values[0][0].str(10))


if __name__ == "__main__":
    verify()
