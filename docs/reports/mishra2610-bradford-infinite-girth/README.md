# Bradford sparse lamps and the four-generator construction

`check.py` uses exact integer coordinates and determinant-one matrices modulo
integer moduli. Run it with Python 3.12 or newer. Its output is `output.txt`.

The source is Henry Bradford, *Quantifying lawlessness in finitely generated
groups*, arXiv:2112.08875v2, Section 4, published in Journal of Group Theory 27
(2024), 31–59, DOI 10.1515/jgth-2022-0113. The target is Pratyush Mishra,
*Infinite Girth, Finite Quotients, and Nielsen Connectivity in Groups*,
arXiv:2610.10925v2, Question 6.12, first sentence.

The fidelity checks reproduce the distinct-difference conditions on finite
prefixes of Bradford's sparse positions and his word-evaluation formula. The
matrices are `(1,2;0,1)` and `(1,0;2,1)`. The program verifies directly that
every nonempty reduced word of length at most six survives in their finite
matrix group.

For radii two and three, the construction check forms every required
single-site commutator word from aligned original lamps, constructs the lamp
Y from those words, recovers the original generators, and evaluates every
reduced word in the four new generators. It also compares the actual wreath
product value with the proposed free-word label calculation. The moduli are
`2 * 3**ell + 1`, with `ell = r * (2*m + 2)` and `m = 3*(2*r+1)`.

These are finite truncations with a common finite lamp group. They test the
construction and signs, rather than asserting that these finite models are
lawless or that they are complete Bradford groups. In particular, the check
does not enumerate every word of length `ell`, prove Sanov's theorem, verify
an infinite sequence of detecting pairs, or settle Mishra's question. The
infinite-family argument remains a mathematical proof obligation.

The check program is original code under Apache-2.0; no third-party program
source is incorporated.
