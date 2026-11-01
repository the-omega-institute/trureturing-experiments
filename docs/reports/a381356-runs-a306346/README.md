# A381356 and the RUNS transform of A306346

`check.py` performs deterministic finite checks using Python 3 integers and the
standard library. Run it from this directory with `python3 check.py`. The saved
output is `check.out`; a successful run ends with `ALL_CHECKS_OK`.

The fixtures come from Paul D. Hanna's [A381356](https://oeis.org/A381356),
revision 11 (2025-03-04), [A381587](https://oeis.org/A381587), revision 36
(2026-03-26), and [A306346](https://oeis.org/A306346), revision 36 (2026-04-03).
A381587 was authored by Neal Gersh Tolunsky and A306346 by Michel Lagneau;
Hanna supplied the identity and the cited row-generation programs.
OEIS material is subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

The program reproduces all complete displayed initial rows, the displayed
sequence prefixes, the 9,259 entries in the first 20 flattened rows, and the
source's 10,797-entry stable prefix of row 25. It independently implements RUNS
using both equal-value groups and run endpoints. It compares the flattened
triangle construction with the row-only recurrence through row 27 and checks
the proposed seed-replacement and reversal relations for stages 1 through 24.
Completed runs are checked against every later tested stage. The final run of
a finite word is excluded from the stability assertion because it can extend.

These finite checks support fidelity and test the proposed identities. They
do not establish existence of the infinite RUNS transform or convergence of
the rows. The planned proof must establish that the A306346 prefix limit has
no eventually constant tail, then pass to each finite collection of complete
runs. No boundedness conjecture about the sequence's values is assumed.

Producers: codex-cli research seat, flight `fbdeed60ade0d7baf714c072`, using
`formal-thinking-and-answer`; candidate supplied by search seat
`codex-search-a-r214` for orchestrator session
`07b67e50-10a4-41b3-9622-660f37c654a5`. Prior exposure: `repo-prior-exposed`.
