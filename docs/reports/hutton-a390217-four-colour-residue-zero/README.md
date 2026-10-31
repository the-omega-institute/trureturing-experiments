# Hutton A390217: four-color partitions of total color zero modulo four

Run `python3 docs/reports/hutton-a390217-four-colour-residue-zero/check.py` from the experiments repository. It uses only Python 3's standard library and exact integers. The stored output is `check.out`; success ends with `ALL_CHECKS_OK`.

The source is Thomas Hutton's [OEIS A390217](https://oeis.org/A390217), revision 29, with both generating functions labeled conjectured. The embedded 57 values are Sean A. Irvine's [b-file](https://oeis.org/A390217/b390217.txt), indexed from 0 through 56. Source data attribution: OEIS contributors, [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

An unbounded multiplicity dynamic program counts size/color multisets by energy and color residue. A separate implementation computes ordinary partition numbers by Euler's pentagonal recurrence and obtains two-color and four-color counts by convolution. Literal multiset enumeration for n from 0 through 8 checks the DP independently; the four residue-zero multisets at n = 2 agree individually with the source example. All 57 source values agree, including all 35 displayed entry values.

For every n from 0 through 128, the check verifies

$$4a(n)=p_4(n)+\mathbf 1_{2\mid n}p_2(n/2)+2\mathbf 1_{4\mid n}p(n/4).$$

It also verifies the unsigned, half-turn and quarter-turn coefficient evaluations separately. In particular, the real quarter-turn coefficient is the one-color partition count at n/4 when four divides n, and zero otherwise. These are finite checks; they do not prove the identity for every n or certify a Lean settlement.
