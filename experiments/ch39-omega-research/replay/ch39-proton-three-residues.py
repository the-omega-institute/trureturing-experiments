"""三个重矢量极点对全部电流分布作证，并检验静态最优解在正动量处的符号。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx

static = run_path(str(Path(__file__).with_name("ch39-selected-phase-proton-exchange.py")))
heavy = static["heavy"]
L,G = heavy["laplacian"],heavy["green"]


def projection_floor(modes,coeffs,phase=None):
    half = len(modes)
    assert len(coeffs)==2*half
    def minimum(factors):
        vals = []
        for u in factors:
            for col in zip(*modes):
                vals.append(sum(a*(u*z).real+b*(u*z).imag
                                for z,a,b in zip(col,coeffs[:half],coeffs[half:])))
        return min(v.lower() for v in vals)
    lower = minimum((1,))
    if phase is not None:
        lower = (5*lower+21*minimum((phase,phase.conjugate())))/26
    return lower/sum(abs(v) for v in coeffs)


def certify_residues(phase):
    eigen,vecs = L.eig(right=True)
    order = sorted(range(len(eigen)),key=lambda i:float(eigen[i].real.mid()))
    lambdas = [eigen[i].real for i in order]
    assert all(a<b for a,b in zip(lambdas,lambdas[1:])) and lambdas[0]>0
    chosen = (0,3,6)
    modes = []
    for j in chosen:
        vect = [vecs[i,order[j]] for i in range(19)]
        norm = sum(abs(z)**2 for z in vect).sqrt()
        rotation = vect[0].conjugate()/abs(vect[0])
        modes.append([(rotation*z/norm).conjugate() for z in vect])

    assert projection_floor(modes[:1],[1,0]) > arb(".198")
    assert projection_floor(modes[1:],[-3944,6738,89318,0]) > arb(".0032")
    assert projection_floor(modes[:2],[95911,0,-3520,569],phase)>arb(".00028")
    assert projection_floor([modes[0],modes[2]],[85825,-9942,-3151,-1082],phase)>arb(".000109")
    floor = arb(".0032")*arb(".000109")
    momenta = [arb(2)**(i-9) for i in range(19)]
    def denominator(x):
        return prod(x+l for l in lambdas)
    lower = arb(1)
    for j in chosen:
        pole = -lambdas[j]
        d = prod(pole+l for i,l in enumerate(lambdas) if i!=j)
        bound = arb(0)
        for i,x in enumerate(momenta):
            basis = prod((pole-y)/(x-y) for k,y in enumerate(momenta) if k!=i)
            bound += abs(denominator(x)*basis/d)
        lower = min(lower,(floor/bound).lower())
    print("Three squared masses",*[lambdas[j].str(14) for j in chosen],
          "one residue above",floor.str(8),"finite Euclidean sample lower",lower.str(12))
    assert lower > arb("2.65e-44")


def threshold_profile(phase):
    weights = static["attaining_weights"](phase)
    normalized = acb_mat(19,1)
    for j,u,w in zip((13,17,16,1),(1,1,phase,phase.conjugate()),weights):
        normalized[j,0] = u*w
    derivative = (-(G*G*normalized))[1,0]
    assert derivative.imag > arb(".01144") and derivative.real < -arb(".00109")
    det = (-L).charpoly()
    # Beginning with I*s^18 in adj(L+sI), use the coefficient recurrence.
    column = normalized
    for k in range(18,0,-1):
        assert column[1,0].real < 0 and column[1,0].imag > 0
        column = det[k]*normalized-L*column
    assert abs(column[1,0]) < arb("1e-10")
    print("Static optimum, slope at zero",derivative.str(15),
          "all eighteen nonconstant numerator coefficients in the second quadrant")


def six_vanishing_coefficients(phase):
    support = [6,9,17,21,24,28,32,35,38,41,42,48,50,51]
    factors = [(1,phase,phase.conjugate())[j//19] for j in support]
    matrix = [[arb(int(j<19)) for j in support],[arb(int(j>=19)) for j in support]]
    power = G
    for k in range(6):
        row = [power[13,j%19]*u for j,u in zip(support,factors)]
        matrix.extend([[v.real for v in row],[v.imag for v in row]])
        power = power*G
    masses = [[arb(5)/26],[arb(21)/26]]+[[0]]*12
    values = arb_mat(matrix).solve(arb_mat(masses)).entries()
    assert all(w>arb(".00022") for w in values)
    reverse = sum(w for j,w in zip(support,values) if j>=38)
    next_coefficient = sum(power[13,j%19]*u*w for j,u,w in zip(support,factors,values))
    assert arb(".38442495474")<reverse<arb(".38442495476")
    assert abs(next_coefficient)>arb("1.49e-11")
    print("Order-six zero, reverse weight and coefficient",reverse.str(16),next_coefficient.str(13))


if __name__=="__main__":
    ctx.prec = 192
    vec = static["selector"]["prepare"]()
    phase = (vec[1,0]/abs(vec[1,0]))**5
    certify_residues(phase)
    threshold_profile(phase)
    six_vanishing_coefficients(phase)
