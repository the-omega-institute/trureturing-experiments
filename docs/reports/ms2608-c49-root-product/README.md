# Square-root Demazure root-product check

Source: Eric Marberg and Travis Scrimshaw, *Square root crystals and the square root of B(infinity)*, arXiv:2608.11009v1, Conjecture 4.9. This independent integer-only program implements the source's standard-set crystal, tensor rule and bottom-to-top column reading. No external program source is used. Run `python3 check.py`; the expected final line is `ALL_CHECKS_OK`.

The source-fidelity checks reproduce the standard S_3 character and all 17 tableaux in Figure 3.1. For w = 3412, direct singleton-tableau closure agrees with signed isobaric operators for lambda_m = (4m,3m,2m,m), m = 1,...,6. Their dimensions are 30, 184, 638, 1645, 3536, 6720. Direct square-root closure, including all intermediate set-valued tableaux, gives 179 and 1629 tableaux for m = 1,2; its excess-zero character agrees with singleton closure. The stable formula min(B+1,A+C+1) is checked on the stated coefficient boxes. All 64 root subsets fail the six coefficient requirements. Every consecutive four-coordinate placement is checked in ranks 4 through 7, including square-root closure, and 27 positive three-gap patterns are checked.

These are finite exact checks. They do not prove the cofinal-limit bridge, the uniform small-string classification or the infinite-rank family. Those obligations belong to the subsequent proof. The reported square-root string enumeration at abc is explicitly bounded; the proof must justify its uniform bound. The source's tensor convention, including its strict inequality for lowering, is essential.

Producer: independent codex-cli research seat, flight 60b1ca2b85aa9deff142d91c; prior disclosure repo-prior-exposed. Search candidate: codex-search-a-r366. This check precedes the Lean proof probe.
