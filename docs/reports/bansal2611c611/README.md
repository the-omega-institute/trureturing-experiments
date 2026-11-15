# Forest Bruhat shellability certificate

`check.py` is an original, standalone Python standard-library program using exact integer arithmetic. Run `python3 check.py`; `output.txt` contains its output. The final line is `ALL_CHECKS_OK`.

Source: Ryan Bansal, *A Quasisymmetric Analogue of the Bruhat Order*, arXiv:2611.02210v1, Conjecture 6.11. Forests are reconstructed from their codes. Forest polynomials are computed directly by enumerating the source's bounded labelings. Monk products are expanded by lexicographic triangular elimination, independently of the trimming-operator coefficient computation.

The fidelity checks reproduce the three displayed Monk identities for the top (2,0,1), then construct and separately verify pure shellings for all 14, 42 and 132 tree tops of ranks 4, 5 and 6. Maximum facet counts are 48, 422 and 5030.

For E=(1,0,2,0,0,1) and Z=(2,2,1,1,1), the program computes all 14 displayed forest polynomials, individual cover labels, the six-by-six weighted matrix and both coefficient calculations yielding 50. It also verifies the seven intermediate trimming identities. Complete enumeration of the standard window P_7 independently finds that (E,Z) consists of two disjoint six-cycles. This is the link of the displayed saturated chain ending at E in the proper-part order complex below Z.

For m=0,1,2, the top Z_m=(m+2,2,1,1,1) is a positive-support zigzag tree with monomial polynomial, and multiplication by h_1 gives exactly P_(Z_(m+1)). The general link obstruction and propagation argument are mathematical obligations of the future Lean delivery; this program is not a Lean proof.

All sparse-polynomial output records are pairs [exponent-vector, integer coefficient], with trailing zero exponents omitted. Certificate vertex indices are one based. The calculation uses the published forest-basis/Monk facts and the support containment reduction to the finite window. No third-party program or Lean source was copied.
