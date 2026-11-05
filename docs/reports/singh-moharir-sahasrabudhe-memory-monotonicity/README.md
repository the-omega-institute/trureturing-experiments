# Heterogeneous memory monotonicity: exact checks

Source: Somya Singh, Sharayu Moharir and Neeraja Sahasrabudhe,
[Opinion Dynamics with Memory Loss and Communication Delays,
arXiv:2608.10491v1](https://arxiv.org/abs/2608.10491v1), Conjecture 1.

Run with Python 3, without third-party packages:

```sh
python3 check.py
```

`output.txt` contains the resulting output. All arithmetic uses integers and
`fractions.Fraction`; every solved linear system has its residual checked exactly.

The program constructs the augmented history matrix directly, solves its full
fixed-point equation, and compares it with a separately constructed reduced
N-dimensional system. The shifts are those forced by the source's concatenated
history state. Its printed block-matrix schematic ends in an all-zero row; that
row cannot literally describe the shift of the last history block. Likewise the
time-domain recurrence and indicator-block indexing differ by one time step.
The check uses the indicator-block numbering in the augmented-system display and
the shifts prescribed by the history tuple. It does not use the source's printed
two-node closed form, whose cross-coefficients do not match its heterogeneous
row factors. These source-display defects do not enter the stationary comparison.

Fidelity reproduces the source's homogeneous formula
p* = c1/(1-M c2) on 1,640 full augmented systems. The exhaustive heterogeneous
check covers N in {1,2}, T in {1,2,3}, every row-stochastic matrix with entries
in {0,1/2,1}, every c_i in {1/8,1/4}, every d_i in {1/16,1/8}, and every
ordered-pair delay set of size M in {1,...,T}, for all 1 <= M <= T. This is
26,108 systems and 592 parameter/memory comparisons, with minimum coordinate
gain 1/105. Since every architecture has the same exact solution at each M,
the comparison covers every pair of delay architectures, including nonnested
sets and differing depths. Another 105 exact three-agent systems use identity,
cycle and reducible networks with T=5.

The numerical coefficients c and d are held fixed while M varies. In an urn
parameterization those coefficients can themselves depend on M; that is a
different question. These finite checks support a preregistration and do not
prove the conjecture for arbitrary dimensions, coefficients or memory sizes.
