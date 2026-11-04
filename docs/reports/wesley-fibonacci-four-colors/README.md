# Wesley's Fibonacci four-color clause

Supplementary exact checks for the first clause of Conjecture 11 in William J.
Wesley, *Improved Ramsey-type theorems for Fibonacci numbers and other sequences*,
arXiv:2211.05167v1, renumbered Conjecture 1 in Integers 25, A29.
The target is Delta(F,3;4)=infinity. The arithmetic-progression clause is excluded.

Run with Python 3; no third-party packages are required:

```sh
python3 docs/reports/wesley-fibonacci-four-colors/check.py
```

The stored output is `check.output.txt`. The program reproduces the source's
Delta(F,2;4)=9 by canonical graph-color backtracking, including an avoiding
coloring on [8] and exhaustion on [9]. It also realizes the stated lower bound
Delta(G,3;2)>50000, with G={F_(3n)/2:n>=1}.

For a positive multiple a of 4, put u_0=0, u_1=1,
u_(n+2)=a*u_(n+1)+u_n, lambda=(a-sqrt(a*a+4))/2, and
alpha=(1-a+sqrt(a*a+4))/4. Define chi(m)=floor(2*frac(alpha*m)).
The program checks the exact recurrence and phase identity for a=4,8,...,40 and
indices 1 through 64. Quadratic algebra uses rational coefficient pairs and
exact sign comparisons. Integer square roots compute all coloring values;
no floating point is used.

For a=4, the four colors are C(n)=1+2*(n mod 2)+chi(floor(n/2)+1).
A dynamic program examines every allowed positive gap and finds maximum
monochromatic F-diffsequence length 2 in each color on [100000]. All 24 distinct
Fibonacci gaps below 100000 are used, and adjacent gaps may differ. Additional
finite two-color checks cover a=4,8,12,16 on [10000].

These finite checks support the proposed construction and expose indexing or
arithmetic errors. They prove no universal statement and are not a Lean proof.
The intended settlement uses a general recurrence identity and a circle-interval
argument. The source's printed Delta(F,3;2)>100000 in its discussion of the
four-color lift is a color-count typo; this report does not endorse that sentence.
