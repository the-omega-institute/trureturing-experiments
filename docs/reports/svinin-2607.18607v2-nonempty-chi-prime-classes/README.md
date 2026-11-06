# Bernoulli χ prime classes: exact finite checks

Source: Andrei K. Svinin, *On a certain arithmetic function defined via Bernoulli numbers*, [arXiv:2607.18607v2](https://arxiv.org/abs/2607.18607v2), Section 3 and Section 4.

Run with Python 3 and SymPy:

```sh
python3 check.py > output.txt
```

The program evaluates the original rational denominator formula directly for every odd prime below 2000. It checks agreement with the source's gcd formula and with

\[
F(m)=\prod_{\substack{q\mid m,\ q\ \mathrm{prime}\\q-1\mid m}}q^{v_q(m)},\qquad \chi(p)=F(p-1).
\]

The Bernoulli convention is B₁ = −1/2. All evidence uses integers and exact rationals; the program performs no floating-point computations.

The output contains 302 prime checks through 1999, 832 checks of the divided-Bernoulli valuation condition, 87 entries of the source's table (five entire rows), and six stated least representatives within the range. Prime 2 is outside the source's domain. The output ends with `ALL_CHECKS_OK`.

These checks support source fidelity and the proposed arithmetic reduction. They do not prove infinitude, a density statement, or the reduction at every index. The unbounded argument belongs to the preregistered proof obligation.
