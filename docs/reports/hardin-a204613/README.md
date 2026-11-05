# OEIS A204613 exact finite check

This report independently enumerates the 64 binary rows and the weighted boundary states for R. H. Hardin's OEIS A204613. The program uses integer arithmetic only. It checks the 134-state/1114-transition compression, reproduces the first five published values, matches all 210 values in the revision-7 b-file, and evaluates the published order-23 recurrence at every index 27 through 210. The output records the counts and ends with `ALL_CHECKS_OK`.

The finite computation is evidence for fidelity and the pre-probe certificate only. It is not a proof that the recurrence holds at every index; the proposed settlement uses the source-faithful state construction, a certified characteristic polynomial, and Cayley–Hamilton for that general step.
