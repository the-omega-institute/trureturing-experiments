# Reciprocal Basel remainder floor

For positive integers n, OEIS A248121 defines D_n = 1/n − π²/6 + Σ_{h=1}^n 1/h² and a_n = floor(1/D_n). Chai Wah Wu's two FORMULA conjectures are the order-five recurrence and the rational generating function. Source: [OEIS A248121](https://oeis.org/A248121), revision 18, with 48 DATA terms reproduced in the program (OEIS data attribution; CC BY-SA 4.0).

Run `python3 check.py` with Python 3 and SymPy. The output is in `check.output.txt`. The script uses integers, rational intervals from Machin's formula and the alternating arctangent series, and exact SymPy algebra. It checks the original analytic expression for every n from 1 through 1000, reproduces all 48 source DATA terms, and checks both rational barrier difference identities, the three residue classes, the recurrence and the generating-function algebra. No floating-point value is used as evidence.

The expected all-index formula is a_n = 2n² + floor(2n/3). The finite interval checks and symbolic algebra support the preregistration's proof outline; they do not supply an infinite-tail proof or Lean kernel verification. The planned main result is a general analytic argument, with utility kind `none`.
