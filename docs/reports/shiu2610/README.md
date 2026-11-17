# Shiu's prime-disc comparison in the field of discriminant −11

Source: Daniel Shiu, *Admissible sets of integers in imaginary quadratic fields of class number 1*, [arXiv:2610.34556v2](https://arxiv.org/abs/2610.34556v2), Equation (3) and Section 4.

Run with Python 3.12 or newer, using only the standard library:

```sh
python3 docs/reports/shiu2610/check.py
```

The deterministic output is in `check.out`; success ends with `ALL_CHECKS_OK` and exit 0. Arithmetic uses integers and `fractions.Fraction` throughout.

Coordinates `(a,b)` denote (a+b\alpha), where (\alpha=(1+\sqrt{-11})/2), (\alpha^2=\alpha-3), and (N(a+b\alpha)=a^2+ab+3b^2). The script enumerates norm discs using the integer identity (4N=(2a+b)^2+11b^2), and checks the enumerations independently in broad coordinate boxes. It computes the full finite quotients by each principal ideal using adjugate-lattice signatures and verifies that every nonzero residue class has a multiplicative inverse. This certifies prime elements of norm four as well as those of prime integer norm.

Source fidelity checks reproduce the Section 4 origin counts 4, 6, 10 at bounds 3, 4, 5, their disc sizes 7, 9, 13, and the failure of the radius-three admissible-set sieve (at most three survivors over all 36 omitted-residue choices). A separate source example in discriminant −163 gives two primes in a five-element disc at bound four.

The counterexample check gives (\pi(0,3)=4<5=\pi(1,3)). The five translated prime elements are (\alpha,1-\alpha,1+\alpha,2-\alpha,2). Their principal-ideal quotients are fields of sizes 3, 3, 5, 5, 4. The cluster covers all three classes modulo ((\alpha)), so it is nonadmissible. Norm integrality extends the same discs to every real (3\le A<4), and negation supplies center −1.

These are exact computational checks and a mathematical certificate description. No Lean proof of Equation (3) or its negation is included. The source's separate question about arbitrarily large bounds is outside this report.
