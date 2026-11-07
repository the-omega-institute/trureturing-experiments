# Bala's five totient-square period conjectures

`python3 check.py` uses Python 3 and the standard library. `output.txt` is its exact output. No floating point arithmetic is used.

The source statements are Peter Bala's conjecture comments on [A281180](https://oeis.org/A281180), [A281181](https://oeis.org/A281181), [A281183](https://oeis.org/A281183), [A281184](https://oeis.org/A281184), and [A281440](https://oeis.org/A281440). Their source revisions are respectively #81, #98, #56, #37, and #13. A281180 starts at n=1; the other four start at n=0.

The program reproduces all displayed terms: 19, 19, 17, 17, and 17 terms respectively. It independently computes rational ordinary coefficients from C'=C^q S, S'=C^(q+1), C(0)=1, S(0)=0, for q=3 and q=7. These agree, after multiplication by the appropriate factorial, with integer polynomial-operator iteration. It also reproduces A281180's displayed 18-term block modulo 7 and verifies that block's least cyclic period is 18.

For each positive modulus k<=60 and each member, it checks 400 consecutive equalities a(n+phi(k)^2)=a(n) modulo k from N(k)=1+max{p^(2r-1):p^r exactly divides k}, with N(1)=1. This gives 120,000 comparisons and generates at least 400 terms per case; the longest generated sequence reaches n=3823. It checks the characteristic-p operator identity delta^p=c delta for p in {2,3,5,7,11,13,17,19} and P=(1+y^2)^2 or (1+y^2)^4. The output also gives nonintegral logarithmic coordinates that prevent direct application of the screened integral-logarithmic-weight theorems.

The modular calculations use finite rings with stable ideals. For A281180, D=(1+y^2)^2 d/dy acts on (Z/kZ)[y]/(y^k). For the other members, E=D^2 acts on (Z/kZ)[c]/(c^k-1), with

E(c^j)=j(j+q)c^(j+2q)-j(j+q-1)c^(j+2q-2).

Evaluation at y=0 or c=1 gives the required sequence terms. The initial terms are also compared modulo each k against the exact source values.

A281440's modulo-13 example begins at n=1: its displayed [1,3,5,...] omits the offset-0 term. The actual offset-0 prefix modulo 13 is [1,1,3]. This indexing difference does not change the eventual-period conjecture.

These are finite checks supporting preregistration. They do not prove eventual periodicity, least eventual periods, or literature novelty. No Lean proof is part of this report.
