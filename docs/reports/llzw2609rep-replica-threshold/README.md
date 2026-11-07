# Exact replica-threshold checks

Source: Guiwen Luo, Xiao Li, Yaqiang Zhang and Ke Wang, *Replica Thresholds for Stripeless Erasure Coding Based on Symmetric Block Designs*, arXiv:2609.29149v1, Appendix B and Section IV. The source matrix is transcribed as mathematical data; the program is independently written.

Run `python3 check.py`. Python 3 standard library only; every count, residue and rank calculation is exact integer arithmetic. The expected output is in `output.txt` and ends with `ALL_CHECKS_OK`.

The proposed identity is q*(A,p) = Gamma(A,p) + 1, where Gamma is the maximum over nonempty failure sets F of size at most p and nonempty H contained in F of the minimum, over i in H, of the number of targets t of i satisfying t in F or d_H(t) >= 2.

The program checks the following results:

- All 1,704 small integer degree vectors in the stated ranges decompose into distinct-target subsets.
- All 23 zero-diagonal regular incidence matrices on two through four nodes satisfy the primary-type peeling formula. Every obstructing pair and every replica count up to its core value has a completed legal witness: 2,540 witnesses, including zero replicas.
- All target choices and all legal parity partitions with two objects per primary on the two- and three-node regular representatives are enumerated independently. Exact GF(2) ranks agree with obstruction existence in every case (23 completed placements in total).
- Appendix B's complement-of-Fano representative has Gamma(A,2) = 3. A fully completed q=3 placement has lost-object observation rank 56 < 64; exhaustive primary-type peeling over all 28 nonempty failure sets of size at most two succeeds at q=4. This reproduces the source's threshold 4 and compares it with Gamma+1.
- A point-plane incidence representative of PG(3,3) is constructed using exact finite-field arithmetic and a nonincidence matching with a prescribed private 3-cycle. All 10,700 nonempty failure sets of size at most three give Gamma = 8. No separated triple exists. A completed q=8 ambiguity uses exactly 338 objects per primary and 13,520 objects overall, agreeing with the proposed threshold 9.

Every constructed witness is checked independently for distinct admissible targets, exactly one use of each replica, complete parity groups, equal object counts at all primaries, and identical surviving observations under a nonzero data change. The checks are finite evidence for the proposed general argument, not a proof of the universal identity and not Lean kernel evidence. The PG(3,3) computation does not establish the whole prime-power family.
