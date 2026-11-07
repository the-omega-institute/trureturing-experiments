# Reciprocal Rado numbers at odd prime powers

Run `python3 docs/reports/gr2607p52-reciprocal-rado-prime-power/check.py` from the repository root. Python 3's standard library suffices. `output.txt` is the complete deterministic output; a successful run exits 0 and ends with `ALL_CHECKS_OK`.

The source is Collier Gaiser and Mojtaba Ramezanpour, *A sharp lower bound for some reciprocal Rado numbers*, arXiv:2607.04373v1, Theorem 3.1, Remark 3.2 and Problem 5.2. Repetitions of denominators are permitted.

The check reproduces the avoiding coloring on [1,34] in Remark 3.2 by direct enumeration, and the Theorem 3.1 coloring at k=3,5. It searches the proposed extension through 3k²+2 for p=3,5,7 and m=1,2, with k=p^m. Every target and both colors are covered. The search enumerates sorted denominator multisets with repetitions. Exact sum bounds prune impossible branches. A denominator-valuation obstruction prunes further branches; where exactly one available denominator has maximal positive valuation, all its possible multiplicities are enumerated. No floating-point arithmetic is used.

The search algorithm is cross-checked against direct enumeration in 240 small cases. A separate exhaustive check verifies the LCM/reduced-denominator principle on 21,942 multisets. Additional checks cover the two exact tail identities for k=3..100 and residual valuations for p=3,5,7,11,13 and m=1..4.

Negative controls recover the old block coloring's monochromatic repeated-denominator solution at 3k². Extending the proposed coloring through 150 at k=7 produces 1/21=2/140+5/150. This is a boundary of that coloring, not a proof that every coloring of [1,150] has a solution.

These are exact finite checks supporting a preregistration draft. They neither prove the uniform bound nor certify the reported source table or any Lean theorem. The planned main result uses a uniform arithmetic argument, not positive finite enumeration.
