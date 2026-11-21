# Erdős–Nicolas 1989: reciprocal divisor-gap sum

Companion program for P. Erdős and J.-L. Nicolas, *On functions connected with prime divisors of an integer*, Number Theory and Applications (Banff 1988), NATO ASI Series C 265, Kluwer, 1989, pp. 381–391, conjecture on p. 390: for n > 5040, Ĥ(n) = Σ_{d<e, d|n, e|n} 1/(e−d) < τ(n). The Lean proof is in the main repository under `D5/S3/Arith/Divisors/ErdosNicolas/`.

    python3 h7560_pair_sums.py    # exact h(7560) = Ĥ(7560)/τ(7560) in (0.99275449, 0.99275450); integer floor sums Σ ⌊10¹⁰/(e−d)⌋ for 2¹⁰3⁴5³7⁴ and 2¹⁰3⁴5⁴7³

The two floor sums certify Ĥ(2¹⁰3⁴5³7⁴) − Ĥ(2¹⁰3⁴5⁴7³) ∈ (0.04146, 0.04159) (each sum has rounding error below C(1100, 2)), so moving exponents to smaller primes does not always increase Ĥ.
