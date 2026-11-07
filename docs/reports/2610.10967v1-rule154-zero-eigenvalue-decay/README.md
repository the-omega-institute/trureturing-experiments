# Rule 154 zero-eigenvalue ratio

Exact checks for Q. Zhu and S. Richard, *Elementary cellular automata: spectral properties of Koopman operator for finite periodic realizations*, arXiv:2610.10967v1, Section 4.2, equation (4) and Figure 13. Source: https://arxiv.org/abs/2610.10967v1.

Run with Python 3, using only its standard library:

    python3 check.py > output.txt

The program checks the Wolfram truth table, reproduces the Rule 74 example at size 4, and enumerates the entire Rule 154 functional graph for every source-defined size from 4 through 22. Repeated deletion of vertices of indegree zero leaves exactly the cycle vertices; the deleted count is the algebraic zero multiplicity of the pullback matrix. No proposed parity-sector reduction is used in these full graph counts.

All counts agree with the candidate identity: for odd n the nonperiodic count is zero; for n = 2m it is 2(2^m − 2^(m−d)), where d is the largest power of 2 dividing m. The parity-sector compression is separately checked for 2 ≤ m ≤ 16. Distinct-input collision localization and invariance of the sector union are checked for 4 ≤ n ≤ 12. The formula gives ratios below 1/100 over the source's complete reported range, 13 ≤ n ≤ 27, and matches the eight-decimal Figure 13 annotations for n = 22 and 24. The compressed graph checks extend to n = 32. Formula-only inequalities are checked through n = 512.

These are finite exact checks, not a proof for unbounded sizes. Counts beyond n = 22 concern the compressed model and depend on the proposed reduction when interpreted as full Rule 154 counts. Every numerical assertion uses integers or rational numbers; the printed source decimals are parsed as exact rational annotations and compared with an exact rounding interval. Successful output ends with `ALL_CHECKS_OK`.
