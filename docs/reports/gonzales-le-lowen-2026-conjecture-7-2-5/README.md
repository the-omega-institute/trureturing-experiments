# Gonzales–Le–Lowen, Conjecture 7.2.5

Source: Teddy Gonzales, Tuong Le and Chayim Lowen, *A quasisymmetric analog of Grassmannian Schubert varieties*, [arXiv:2609.30257v1](https://arxiv.org/abs/2609.30257v1). The conjecture asks for an isomorphism of specified group objects in the weak homotopy category between the toric direct limit and the loop suspension of infinite complex projective space.

Run with Python 3.12 or later, using only the standard library:

```sh
python3 check.py
```

The stored output is `check.out`. Every calculation uses integers, including negative coefficients in the free associative algebra.

- **Source fidelity:** Section 2.4 states the ribbon identity $R_sR_t=R_{s0t}+R_{s1t}$. The script checks it in the free integer algebra on complete generators for all 4,097 pairs of positive degrees with total degree at most 10. It independently checks the ribbon convention against direct semistandard word enumeration in 189 finite-alphabet cases.
- **Necessary comparison identities:** Through half-degree 10, the stable ranks of the alternating snake filtration, computed from distinct subsequences as in the source's Theorem IV and Remark 7.2.2, agree with the composition ranks of the James tensor algebra. Empty words contribute the degree-zero unit.
- **Finite product model:** For every $N=0,1,2,3$ and $m,n=0,1,2,3,4$, it checks basepoint deletion on the cellular homology bases of $(\mathbb{CP}^{N})^m\times(\mathbb{CP}^{N})^n$, degree preservation, a group-level padding section, and agreement with the proposed wedge-of-smash-powers graded ranks. This covers 131,908 source tuples and 35,440 reduced target basis elements.
- **Multiplication and coalgebra:** It checks ordered tensor multiplication and the complete-generator coproduct through half-degree 10, including a noncommutativity control.

These are finite checks of necessary algebraic identities for the proposed positive settlement. The padding section in the script is a section of free abelian groups. It does not construct a continuous section after suspension. No map into the source's toric direct limit, homotopy, homotopy equivalence, or weak H-group comparison is certified. Passing this program does not settle the conjecture or replace its binding topological statement.

The planned mathematical argument uses strict-unit rectification, a finite James-product suspension section, suspension detection into an H-space, and homological Whitehead. Its main content is an unbounded topological argument; this program supplies auxiliary checks, not a positive exhaustive-enumeration settlement.
