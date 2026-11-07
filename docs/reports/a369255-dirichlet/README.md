# A369255: two Dirichlet convolution identities

Run `python3 check.py` with Python 3. The program uses integer arithmetic only and has no external dependencies. `output.txt` is its output; successful execution ends with `ALL_CHECKS_OK`.

The source is [OEIS A369255](https://oeis.org/A369255), revision 15. The two conjectures are attributed there to Sequence Machine and recorded by Antti Karttunen. The source reports agreement through its first 2^20 terms. The check reproduces that entire range using literal divisor counts for A038548 and A007875, divisor summation for A140773, the published fourth-power Möbius sum for A219009, and the multiplicative prime-power definition of A355837. Six embedded published prefixes provide a separate interpretation check.

It also checks the proposed prime-exponent support classification and the ordered-triple/square-divisor counting formula through 2^20, then tests 74,613 exponent multisets with at most six positive exponents, each at most 16. These are finite checks, not a proof of either unbounded conjecture. In particular, the parity sequence is not multiplicative: a(2)=a(3)=0, while a(6)=1.

Source data attribution: OEIS Foundation Inc. and the contributors to A038548, A140773, A007875, A219009, A355837 and A369255. The short embedded data prefixes are used under the [OEIS End-User License Agreement](https://oeis.org/LICENSE). The original check program is part of this repository under its applicable license.
