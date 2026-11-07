# Svinin forward gcd experiment

Exact checks for Andrei K. Svinin, arXiv:2611.02309v1, Section 6, the example with initial pairs (1,1), (1,1), (2,3) and coefficients (a,b,c,d)=(4,2,1,4).

Run `python3 check.py --limit 2000`. Dependencies: Python 3, SymPy 1.14.0 (BSD); optional gmpy2 2.3.2 (LGPL-3.0) supplies exact integer acceleration, with a standard-library integer fallback. The script does not copy third-party implementation source. No floating-point value is used as mathematical evidence.

The output reproduces all three source tables at indices 0 through 8. Independently generated integer coordinates pass exact division, positivity, complete gcd and nondecreasing exponent checks at every index 0 through 2000. A separate rational implementation and scalar Somos-5 recurrence agree through index 100. Symbolic invariant/error identities and all 32 rational residue certificates also pass.

For n=10k+r, 0<=r<10, the checked valuation formulas are X=5k^2+rk+a[r] and Y=X+2k+c[r], where a=(0,0,1,0,2,1,2,3,3,5) and c=(0,0,-1,4,-1,1,1,0,4,0). The checked gcd is 2^min(X,Y). The final measured values are X_2000=200000, Y_2000=200400. The source's backward continuation has y_-1=5/2; all gcd assertions here concern nonnegative indices.

`check-output.txt` contains the complete successful output, ending ALL_CHECKS_OK. These are finite experimental results supporting a proposed unbounded argument; they are not a proof-assistant settlement or a substitute for that argument.
