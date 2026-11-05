# Exact recovery-profile and relaxation checks

Run `python3 check.py` from this directory. Python 3 with its standard library is sufficient. The output is stored in `check.out`; successful execution ends with `ALL_CHECKS_OK`.

The source is D. Bar-Lev, H. H. López, G. L. Matthews and F. Salizzoni, *A Matroidal Framework for Random Access in Linear Codes*, arXiv:2610.28820v1, Definitions 3, 6 and 7 and Section V. The check uses the source's sum over subset sizes **1 through n−1**, including for matroids with loops.

The program verifies the source's sparse-paving formula for a_i(s), and its equivalence between recovery balance and equal circuit-hyperplane degrees, on 403 isomorphism representatives. It checks every circuit-hyperplane relaxation of every database representative of rank 2 through n−2 for n ≤ 8: 3,583 pairs. It checks basis exchange, the change of exactly one subset rank, both exact rational profile changes, and decoding the relaxed hyperplane as the set of coordinates attaining the smaller profile value. There are 218 pairs with a balanced source. The rank and profile checks cover arbitrary matroids in the database, including those which are not sparse paving.

The finite counts below are **labelled counts of all ranks**, obtained by weighting each isomorphism representative by n!/|Aut(M)|. Automorphism orders are computed with integer arithmetic by enumerating permutations of clone classes and checking a canonical lift against the full basis family. The program independently cross-checks automorphism orders by enumerating all permutations, and both labelled counts by enumerating every basis-family candidate, for n ≤ 5. The total labelled counts also agree with OEIS A058673; the OEIS values are a cross-check, not program inputs.

| n | balanced labelled | all labelled | exact ratio |
|---|---:|---:|---:|
| 0 | 1 | 1 | 1 |
| 1 | 2 | 2 | 1 |
| 2 | 3 | 5 | 3/5 |
| 3 | 4 | 16 | 1/4 |
| 4 | 8 | 68 | 2/17 |
| 5 | 6 | 406 | 3/203 |
| 6 | 127 | 3,807 | 127/3807 |
| 7 | 68 | 75,164 | 17/18791 |
| 8 | 57,698 | 10,607,540 | 28849/5303770 |

The source assumes positive n; n=0 is an extra sanity check. Finite counts and checks do not establish the asymptotic conjecture. Database completeness is an external input. The proof outline's further external input, Pendavingh–van der Pol Theorem 5.3 (arXiv:1602.04763v3), is not tested by this finite experiment.

## Database source and license

The files in `data/` are the required rank-at-most-half, n ≤ 8 subset of `matroid-database` 0.3, distributed under GPL-3.0; see `DATA_LICENSE.txt`. The database originates with Yoshitake Matsumoto's *Database of Matroids* and the package also incorporates output from `gmou3/matroid-generator`. Higher ranks are obtained by duality; ranks zero and n are generated directly. Each input is a full basis/nonbasis string in colexicographic order. No package code is vendored. The program prints SHA256 hashes for every data file it consumes.

- Package: https://pypi.org/project/matroid-database/0.3/
- Upstream: https://github.com/gmou3/matroid-database
- Wheel: https://files.pythonhosted.org/packages/py3/m/matroid_database/matroid_database-0.3-py3-none-any.whl
- Wheel SHA256: `85d0304575784ceb4797014fb5a761443d2a0bbf01aafc38ef293d2c64a1b5ce`
- Original database: https://www-imai.is.s.u-tokyo.ac.jp/~ymatsu/matroid/index.html
- Source paper: https://arxiv.org/html/2610.28820v1
- Counting input: https://arxiv.org/abs/1602.04763v3
