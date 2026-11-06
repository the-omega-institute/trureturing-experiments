# Order observations for Eidesen Question 11.1

Run `python3 docs/reports/eidesen-q111/check.py`. Python 3's standard library is sufficient. The checked output is in `output.txt`; success requires exit 0 and final line `ALL_CHECKS_OK`.

The source is Jonas Eidesen, *Projective error models: Stabilizer codes, Clifford codes, and weak stabilizer codes*, arXiv:2506.01843v3, Definitions 3.2, 5.2, 5.8, 6.1; Propositions 4.4, 6.7, 8.1 and 10.1; Question 11.1. Source URL: https://arxiv.org/abs/2506.01843v3. The program is an independently written implementation of the displayed matrices and constructions; no external source code or paper text is included.

For t=1,2,3,4, put r=2t, m=r^2+1 and N=4m. A monomial matrix is represented by a permutation p and integer phase exponents e, meaning M e_j=z_N^{e_j} e_{p(j)}. Roots of unity are formal: equality is equality of exponents modulo N, so all evidence uses exact integers. The four factors are A=C_2 x D_4 on C^4, B=D_m on C^2, Q=D_4 on C^2, and F=C_m semidirect C_4 on C^4 (action a -> a^r). Here |D_n|=2n.

For every tested parameter the program checks each factor's projective multiplication law for every pair, projective injectivity, and a certificate of irreducibility: a diagonal generator has distinct weights, and the permutation generators connect all coordinates. It checks the induced-code certificate by disjoint orthogonal translates partitioning the ambient coordinate basis and irreducibility of the restricted action. These facts reproduce the source's Proposition 8.1 at n=2, including |L(U_A)|=8 and |S(U_A)|=1 for U_A=span(e_0,e_1). They also check the other induced blocks.

In the shared product model on C^64, the codes are

- W_good=V_A tensor V_B tensor V_Q tensor span(e_0 in V_F);
- W_bad=span(e_0,e_1 in V_A) tensor span(e_0 in V_B) tensor V_Q tensor V_F.

Both have dimension 16. At m=5 the program constructs all 25,600 tensor matrices, checks their projective distinctness, counts L and S directly from the projections onto the coordinate supports, and computes the maximal common eigenspaces for the characters actually acting on each code. It obtains L=6400, S=5 for both; the eigenspace dimensions are 16 and 32. For m=17,37,65 the group orders and code invariants are computed from the fully enumerated factors using the source's product formulas. These larger product models are not fully enumerated.

The source's weak-stabilizer criterion identifies the code with its maximal scalar-character eigenspace. Thus the finite checks support opposite weak-stabilizer status with the same order observations. They do not establish the infinite family by exhaustive computation, authenticate the literature, or constitute Lean verification. The general family requires the mathematical argument and subsequent kernel checking specified by the preregistration.
