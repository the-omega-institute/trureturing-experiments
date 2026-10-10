# Cha–Lee: the canonical allocation is MS-GC optimal for balanced QFT distributions

`check.py` checks the conjecture of arXiv:2501.11816v3 (Section 4, "QFT circuits") that the canonical allocation π*(q_i) = p_⌈i/m⌉ is optimal for MS-GC.

1. It implements the source's definitions literally: QFT event times in the conventional order (final SWAPs omitted), candidate migrations, and home and joint coverage with the latest-slot rule. On 2400 random migration sets it checks that coverage of all non-local gates is equivalent to R_i ∩ L_j ≠ ∅ for all i < j.
2. For (k,m) ∈ {(3,1),(3,2),(4,1),(4,2),(5,1)} it computes the exact MS-GC optimum of every balanced allocation, up to module relabelling, by exhaustive search, and checks that none is below OPT(π*).

Run: `python3 check.py` (exit 0, final line `ALL_OK`).
