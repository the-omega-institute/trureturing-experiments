# Tierz BPS matrices — 19 October 2026

Full sparse integer matrices from the exact rational diagonal-action computations
for (p,N) = (5,3), (5,4), and (7,4). These are the unchanged matrix files from
`computation-supplement_20261019.zip`, published separately to keep the email
supplement small and readable. The supplement retains the reproduction code,
exact block ranks, summaries, tables, and independent verification reports.

Paths below are relative to this directory. To run
`python3 review/check_matrices.py` from the extracted supplement, download the
files into the matching `results/p*-n*/` paths. The code can also regenerate
them; see `code/README.md` in the supplement for dependencies and commands.

Each gzip file contains JSON Lines, one matrix block per record, with its
exterior-monomial bases and exact rational rank. Sparse entries use zero-based
`(row, column, value)` triples. For stacked maps, `domain1` precedes `domain2`
and `target1` precedes `target2`. The matrices represent Q/p and the stacked
maps for the diagonal creators a and a+b; the supplement defines conventions.

| File | Bytes | Matrix blocks | SHA-256 (compressed file) |
| --- | ---: | ---: | --- |
| `results/p5-n3/matrices.jsonl.gz` | 3,822 | 260 | `4d3623e4388c307b7d22d528264a0dee40749db769e90de1918dd6c21f08759c` |
| `results/p5-n4/matrices.jsonl.gz` | 2,285,552 | 5,937 | `9cf4d84c24a827a9a9e6612a84e7d07d16478c376e1a32c280f27bd1300a8d1e` |
| `results/p7-n4/matrices.jsonl.gz` | 1,411,970 | 5,511 | `755bd482f83721cca27c51970d644f915319876b1267bf917ac56265dea7e285` |

These files contain computational evidence, not a new Lean formalization.
The model and source comparison data are from Miguel Tierz,
*BPS spectra of Tr[Psi^p] matrix models for odd p*, arXiv:2604.27164v2.
The matrices were constructed independently of the source's rank arrays.
