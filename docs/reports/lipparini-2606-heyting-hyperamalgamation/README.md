# Heyting hyperamalgamation candidate

`check.py` uses only exact integer operations.  For every chain size `m = 3,…,30` and every interior `e`, it checks the Heyting implication preservation of the graph embedding `x ↦ (x, χ(x))`, the coordinate product strong amalgam, the forced mixed inequality, and the stalled transfer invariant through lengths `1,…,30`.  The generated JSON records the counts and the terminal line of the program is `ALL_CHECKS_OK`.

The construction is the family described in the preregistration for Lipparini, *A hyperamalgamation property*, arXiv:2606.09329v1, Problem 6.6(a).  It is a research-stage exact check; it does not claim a Lean proof.
