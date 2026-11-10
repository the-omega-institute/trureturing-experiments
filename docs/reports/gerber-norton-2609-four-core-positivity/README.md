# Gerber–Norton four-core positivity

Source: Thomas Gerber and Emily Norton, *Enumerating cores of charged multipartitions*, arXiv:2609.03738v1, Table 2, level 4, d = 1, e = 4, charge (0,0,1,2). The table marks positivity at every nonnegative size as conjectural. The count concerns total size, not rank.

Run with Python 3.12 or newer:

```sh
python3 check.py > check.out
```

The program uses integers and integer square roots throughout. It checks the size-16 example and the original 4-core bipartition example following Lemma 2.11. Independent direct enumeration reproduces Table 2's equal-charge series coefficients through size 12 and the stated zero at size 4 for charge (0,0,0,1). Target counts through size 12 are [1,2,2,3,4,3,5,7,2,5,7,4,7].

The proposed construction uses a charge-zero beta set B with residue-class thresholds t, then the rows B, B, B union {b_r}, B union {b_r,b_s}. All 2,772 choices with each t_i between -3 and 3, sum t_i = 0, and distinct r,s are checked by reconstructing partitions, computing charges, checking exact nesting, and comparing size with the quadratic identity. The eight possible sorted residue multisets for a vector with sum 1 modulo 4 are checked exhaustively. A search for three odd squares followed by the Hadamard and residue construction produces an original core at every size from 0 through 512, including four empty partitions at size 0.

Beta sets store all positions at least -96, with all smaller integers implicitly occupied. For direct enumeration, all partitions have at most 12 rows. For constructed sets, the program asserts thresholds well above the cutoff, reconstructs the zero tail, and checks equality after conversion back to beta sets. Thus the finite representation keeps every changed bead, including changes at negative positions.

Equation (5.5), the displayed Eureka dissection, has constant term 2 on its printed left side and 1 on its right side under the source's counting definitions. This discrepancy is reported and the equation is not used.

These are finite fidelity and construction checks. They neither prove universal positivity nor supply a Lean formalization of the three-square theorem. The planned universal argument needs the 8n+3 case of that theorem as proved content, not as an additional hypothesis.

The program is independently written from the mathematical definitions. No third-party program or Lean proof source is incorporated.
