# Copeland A131758 polynomial content

Run `python3 docs/reports/copeland-a131758/check.py` from the repository root.
Python 3's standard library suffices. All arithmetic uses integers or exact
fractions. The stored output is `check.out`; success ends with `ALL_CHECKS_OK`.

The check reproduces the eight displayed source rows and compares two independent
representations through row 20: the finite Laguerre sum and the inverse Eulerian
binomial transform. At row 8 the content is 1, certified by
−6417·6209 − 962·(−41417) = 1. Both proposed divisor counts D(8) and D(9) are 2.
The nine coefficients were already published in the source b-file; they are not
newly discovered data.

The infinite obstruction is parity. For even n, every nonleading coefficient of
(−1)^n n! L_n(x) is even, so its value at 1, and therefore the row content, is odd.
For n = 2^(12s+7), D(n) = 6s+4 is even. Also n+1 = 3m with 3 not dividing m and
m ≡ 3 modulo 8, so divisor pairing makes D(n+1) even. The program checks finite
instances and residues of this argument; it does not certify an infinite theorem
or any Lean proof.

Sources: [OEIS A131758 revision 48](https://oeis.org/search?q=id:A131758&fmt=text),
[its b-file, indices 36–44](https://oeis.org/A131758/b131758.txt),
[A038548](https://oeis.org/A038548), [A060872](https://oeis.org/A060872), and
[DLMF 18.5.12](https://dlmf.nist.gov/18.5.E12) with α=0 for the ordinary Laguerre
normalization. Source data are subject to the [OEIS end-user agreement](https://oeis.org/LICENSE).
