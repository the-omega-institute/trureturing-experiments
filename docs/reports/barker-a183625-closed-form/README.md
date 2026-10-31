# Barker A183625 closed-form checks

Run `python3 docs/reports/barker-a183625-closed-form/check.py` from the repository root. Python 3's standard library is sufficient. The stored output is `output.txt`; successful execution ends with `ALL_CHECKS_OK`.

The program generates labelled three-column arrays from the original four-entry block equations. It reproduces the first twelve printed terms and all 200 entries of Hardin's b-file, checks the closed form and rational-series product through degree 200, and independently enumerates all ambient arrays of heights one through four. It also checks both directions of odd-row complementation for heights one through three and the profile multiplicities for widths one through eight. Arithmetic uses integers and `fractions.Fraction` throughout. These bounded checks do not prove the unbounded conjectures.

Source: [OEIS A183625 revision 10](https://oeis.org/A183625/internal), revision timestamp 2018-03-31 17:01:05; Hardin's native count and [b-file](https://oeis.org/A183625/b183625.txt), and Barker's FORMULA conjectures. `source-terms.txt` preserves the 200 source pairs retrieved on 2026-10-31. OEIS material is subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE). The check program is independently written and uses no source implementation.
