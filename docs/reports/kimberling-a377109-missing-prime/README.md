# OEIS A377109: missing prime 241

The sequence is the coefficient independent of sqrt(2), sqrt(3), and sqrt(6) in
(2 + sqrt(2) + sqrt(3))^n, starting at n = 0. The source is
[OEIS A377109](https://oeis.org/A377109), revision 25. Its first COMMENTS conjecture
asserts that every prime divides infinitely many terms and that the successive
zero-index gaps are eventually periodic.

Run `python3 check.py`. Python 3 standard-library integer arithmetic suffices;
there are no external dependencies. The saved output is in `check.out`.

The check reproduces all 24 displayed source data terms using an independent
multinomial expansion, four-coordinate multiplication, and the fourth-order
recurrence. Modulo 241 the full four-coordinate state has a zero-free first
coordinate throughout its 240-step cycle and returns to (1,0,0,0). The nominated
240 residues agree with both a scalar recurrence and a separate split-conjugate
calculation of every coordinate. Primality and the two modular square roots are
also checked exactly.

The cycle is evidence for the planned all-index refutation, conditional on the
formal coefficient-identification and iteration-periodicity bridges. This report
contains no Lean proof and does not claim a kernel-verified settlement. The
conditional gap-periodicity clause is not separately refuted.

The source has discrepancies: its displayed exponent-3 example is actually the
exponent-2 expansion; the recurrence lists two initial indices incorrectly; and
the reported p=2 gap (4) disagrees with the defining sequence, whose zero indices
through 20 are 1,3,...,19. The check records these differences explicitly and uses
the defining expression and data, not those lines, for fidelity.
