# Coset description of anti-phase stable matchings

Independent exact checks for Yoshiteru Ishida, arXiv:2607.10306v1, Problem 2.
Run `python3 check.py`; the final line must be `ALL_CHECKS_OK`. Python standard library only.

The program constructs multiplication tables for representatives of every group isomorphism type of order at most 8: cyclic groups, C2 x C2, S3, C4 x C2, C2 x C2 x C2, D8 and Q8. It validates the group laws exactly. It checks every ordering through order 5 and 40 deterministic orderings for each larger representative (497 cases). Every matching is tested by brute force from the preference lists. All arithmetic is integral; the pseudorandom generator only chooses test inputs.

Fidelity checks reproduce the source's C4 ordering (0,2,1,3) with eight stable matchings and its Klein ordering with ten stable matchings and six rotations. Every case also checks the source's rank-sum identity, maximal-coset lemma and counting lower bound. The C6 ordering (0,2,4,1,3,5), with nonnested adjacent cyclic subgroups, has 24 stable matchings and 11 rotations.

For the proposed answer, a separate coset calculation enumerates ideals of the adjacent-overlap poset. Its image is compared as a set of complete matchings with the brute-force stable set, and every male-dominance comparison is checked. Rotations are independently computed from the standard first-accepting-woman successor in each stable matching, without receiving the coset poset or any ideals. Elimination histories are then traversed to recover rotation precedence; every relation is compared with coset-poset reachability. `output.txt` contains the complete aggregate results and a hash of all per-case integer results.

These finite checks support the preregistered universal structural argument; they do not prove it. No third-party program source was read or copied.
