"""An exact finite-response test that distinguishes two global kinetic prescriptions."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
integers=run_path(str(root/"ch39-vote-weighted-integers.py"))
base=integers["base"]
B,round_out,multiply,absolute=integers["B"],integers["round_out"],integers["multiply"],integers["absolute"]


def main():
    rec=base["data"]("ch39-cutoff-threshold-inverses.json")[1]
    w=base["exact_fit"](rec)
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    def ordinate(i):
        h=prefix[i-1][0]
        return round_out(B(Q(2**800,(h+1)**2),Q(2**800,h*h)))
    z=[ordinate(i) for i in range(1,23)]
    A=integers["A"]
    tol=Q(1,2**272)
    W=[multiply(*(B(1)+B(a)*v for a in A)) for v in z]
    # Delta Gram entries are the finite measurement errors times the partial
    # fraction coefficients. Cauchy-Schwarz bounds all cross terms by the
    # largest diagonal absolute-coefficient sum. The reference Gram dominates
    # diag(w_1,...,w_22) on the square-root-rescaled cardinal basis.
    sums=[Q(0)]*22
    for a in A:
        den=prod(abs(1-Q(k,a)) for k in A if k!=a)
        for i,v in enumerate(z):
            other=z[:i]+z[i+1:]
            val=multiply(*[(t+B(Q(1,a)))/absolute(v-t) for t in other])
            sums[i]+=(W[i]*(val**2)).hi/den
    assert 22*max(sums)<2**201
    assert 22*max(sums)*tol<min(w[i].lo for i in range(1,23))
    mean_max=sum(v.hi for v in z)/22
    # The same universal mean bound as in the full mass replay.
    ref=sum((round_out(ordinate(i)*v) for i,v in w.items()),B(0))
    mu_min=ref.lo-Q("4e-55")
    eta=Q(1,128)
    lower=mu_min/(mean_max+eta*z[0].hi)
    assert lower>Q("2.06")
    # In the vote-weighted model sum(w*z/mu)=1 and each positive
    # finite-denominator correction strictly decreases the observable.
    print("At least 22 distinct positive heights are required even with 2^-272 time errors.")
    print("At eta=1/128 the normalized secular difference is >2.06 for the arithmetic mean, <1 for the vote-weighted mean.")


if __name__=="__main__":
    main()
