# Ferreira–Simon Mills derivative quotients

Source: Rui A. C. Ferreira and Thomas Simon, *Convolution of beta prime distribution*, arXiv:2108.09244v2, Section 5.5, Proposition 10, the Sampford corollary and Conjecture 3; published in *Transactions of the American Mathematical Society* 376 (2023), 855–890, DOI:10.1090/tran/8748.

Run from the repository root with Python 3 (standard library only):

```sh
python3 docs/reports/ferreira-simon-mills-derivative-quotients/check.py
```

`output.txt` is the program's output, ending with `ALL_CHECKS_OK`. All evidence uses integer arithmetic and `fractions.Fraction`. Closed interval endpoints have denominator 10^240 and every operation rounds outwards. Machin's identity and alternating arctangent remainders bound pi; integer square roots bound sqrt(pi/2); Taylor remainders bound the defining Gaussian integral and exponential. The Mills ODE generates derivative jets, which are squared and divided with interval arithmetic. Scalar division rounds directly, including division by large factorials.

Fidelity checks certify the source's Sampford bound at x=0, 1/100, 1, and the alternating signs in Proposition 10 at z=1 for n=0..4 and derivative orders k=0..12. These are finite checks of the source's results.

For the uncomposed quotient F_n(x)=−r^(n)(x)^2/r^(2n+1)(x), the program certifies negative alternating Taylor coefficients at x=1/100 for n=0..4, with orders 40, 70, 96, 121, 146 respectively. It also reproduces the strict bounds −2995/10^29 < F_0^(40)(1/100)/40! < −2994/10^29 and −4604/10^55 < F_1^(70)(1/100)/70! < −4603/10^55.

These are finite exact arithmetic certificates, not Lean proofs. They do not certify the proposed analytic extension to every n, every real center and arbitrarily high derivative orders. That extension uses entire-function growth, finite zero cancellation and Pringsheim's theorem in the preregistration's mathematical outline.
