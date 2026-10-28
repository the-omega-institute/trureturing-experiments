# Sahasrabudhe 2015: universal-input morphism parity check

`check.py` uses exact SymPy algebraic arithmetic for the two-dimensional
Sturmian construction in Lemma 1 of Sahasrabudhe's paper, checking the value
`d + 1 = 3` for lengths 1 through 20 on a long exact orbit sample.

It then exhausts all two-letter codebooks with codeword lengths one or two and
integer entries in `{-1, 0, 1}`.  For every codebook satisfying the necessary
common-slope condition, it enumerates all source words of the required finite
length and checks the complementary-factor involution, its fixed point, and
odd factor-sum cardinality at every small admissible multiple of the length
gcd.  This is a finite regression check for the registered proof mechanism;
the quantified theorem is not delegated to computation.
