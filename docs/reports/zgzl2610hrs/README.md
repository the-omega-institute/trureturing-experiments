# HRS deep-hole checks

Run `python3 docs/reports/zgzl2610hrs/check.py` from the repository root.
The program uses only Python's standard library and integer arithmetic in
the prime fields with 3 and 5 elements. Its stored output is `check.out`.

For Zhu, Gu, Zhang and Li, arXiv:2610.16643v1, Theorems 3.1 and 3.3, the
fidelity check enumerates every received word for `(q,s,k)=(3,2,2)`. It
computes distances directly to all nine codewords, obtains covering radius
4 and 18 deep holes, and checks the polynomial remainder criterion on all
729 words. It also checks that Hasse evaluation is bijective there.

The classification check covers all 30 admissible triples with
`q` in `{3,5}`, `s` in `{2,3}` and `max(2,s) <= k <= sq-1`.
It enumerates every bounded multiplicity profile. Necessary linear
constraints on a candidate functional follow from complete affine pencils
of admissible products and, over the field of five elements, an independently
enumerated block of all 125 quadratic functionals. Exact Gaussian
elimination retains all possible solutions; every remaining functional is
then tested against every original profile. The largest residual search has
15,625 functionals. Thus this is exhaustive constraint elimination, with no
random sampling or invocation of Segre's theorem. The program also checks
the source remainder matrix against the complementary-product pairing and
checks that its triangular coordinate transformation is invertible.

Every case has exactly `q-1` deep cosets and `(q-1)q^k` deep words. All 18
overlapping counts agree with sheet1, columns B/C/E/H, of the source's
`NRT deep hole.xlsx`, DOI
[10.5281/zenodo.23279452](https://doi.org/10.5281/zenodo.23279452).
The other 12 cases, with `q=5,s=3`, are absent from that workbook.
The workbook SHA-256 is
`7fbd89fe134d5e4e05199e6c967bc1c068e17654a17be09536f1ebe751363e03`.
The literal overlapping counts are embedded in the program. The workbook
is data licensed CC-BY-4.0 by Zhihao Zhu; no third-party program was read,
copied or ported to produce this check.

These results verify finite instances and the proposed reduction on those
instances. They provide no Lean proof or universal settlement. The
preregistration's proposed universal argument has general algebraic content;
this computation serves only as an independent check.
