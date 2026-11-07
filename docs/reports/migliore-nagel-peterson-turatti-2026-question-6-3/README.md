# Skew-line Rao modules on a smooth quadric: exact checks

The target is Question 6.3 in J. Migliore, U. Nagel, C. Peterson and
E. Teixeira Turatti, *Weak and Strong Lefschetz Properties for Hartshorne-Rao
Modules of Curves in P³*, [arXiv:2605.19434v1](https://arxiv.org/abs/2605.19434v1).
The question asks whether every union of skew lines on a smooth quadric has
a Hartshorne–Rao module with the strong Lefschetz property.

Run `python3 check.py`; Python 3 and SymPy are required. `check.out` contains
the exact output. The final line is `ALL_CHECKS_OK`.

The program uses the algebraic cokernel stated in the proof of Proposition
4.2 of the source:

\[
M(C)_t=\operatorname{coker}\left(R_t\longrightarrow
\bigoplus_{i=1}^r \mathbb Q[v_0,v_1]_t\right).
\]

Each actual line is parametrized by
\((z_0,z_1,z_2,z_3)=(u_{i0}v_0,u_{i0}v_1,u_{i1}v_0,u_{i1}v_1)\),
on the smooth quadric \(z_0z_3-z_1z_2=0\). There are two configurations for
each \(r=2,\ldots,6\): all chosen points of the ruling are affine, or one is
at infinity. All matrices have integer entries, and all ranks are exact over
\(\mathbb Q\). These ranks remain the same after a characteristic-zero field
extension.

The fidelity checks reproduce Hilbert-function symmetry (Lemma 4.1), initial
dimension and support (Proposition 4.2), and WLP (Proposition 4.4). The two
linear forms have coefficient vectors \((1,0,0,1)\) and \((2,3,5,7)\). All 200
tested induced power maps have maximal rank. The program also verifies that
multiplication preserves the restriction image, compares ranks with the
proposed differential model, and extracts Jordan blocks from successive
ranks. With \(d=r-2\), the blocks have length \(d-2j+1\) and multiplicity
\(d-2j+1\), for \(0\le j\le\lfloor d/2\rfloor\).

Negative controls use a rank-one coefficient matrix (a tangent plane) and
the zero linear form. They detect failures for every tested \(r\ge3\).
The single-line cokernel is also checked to vanish in the tested degrees.

These finite computations support a proposed general argument. They do not
prove the universal question, classify all configurations, or certify prior-art
completeness. No Lean proof is included.
