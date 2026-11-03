# Monotone LAB execution: exact finite checks

The implemented finite checks pass. Strict acceptance remains incomplete:
`--require-all` exits 2 because the complete ML-A evaluator execution and
its full simultaneous resource account have not been instantiated. Host
law arithmetic and the acquired successor kernel do not discharge these
requirements. No theorem counterexample was found in the checked scopes.

Python 3.12+, standard library, assertions enabled, the pinned local Git
objects, and an explicit theory candidate checkout are required. The
candidate's theory and bibliography bytes must match `theory-candidate.json`.
`pins.py` reads those actual candidate files and refuses mismatches before
running checks. It never substitutes the sealed old draft. Paths are
explicit inputs; no personal checkout location is assumed.

Set `WORKTREE_PROTOCOL` to the repository's existing
`tools/scripts/worktree/worktree_protocol.py`, and `RRO_THEORY_ROOT` to the
actual theory candidate checkout. From this experiment checkout run:

```sh
python3 -B "$WORKTREE_PROTOCOL" --source "$RRO_THEORY_ROOT" with \
  --path "$RRO_THEORY_ROOT" --read docs/develop/theory \
  --read Library/Computability/nearywoodsmurphyglaschick2014wang.md -- \
  python3 -B "$WORKTREE_PROTOCOL" --source "$PWD" with --path "$PWD" \
  --read rro-structural-recovery --read rro-positive-certificates \
  --write rro-lab-execution -- \
  python3 -B rro-lab-execution/exact_checks.py \
  --theory-root "$RRO_THEORY_ROOT" --output rro-lab-execution/exact-results.json
```

For a deterministic comparison use `--read rro-lab-execution` in place of
its write scope and `--check rro-lab-execution/exact-results.json` in place
of `--output`. Add `--require-all` for full acceptance. Default exit 0 means
only that the implemented checks passed. No randomness or floating point
is used. The result binds all candidate/supplier bytes and program hashes.

| Source | Executed scope |
| --- | --- |
| ML-I/C/L/B, §§1–3 | `lab.py` uses bound object references, physical terminal/unary/binary nodes, guarded P/S, two comb arms, saved heads and actual ancestor routing. An independent integer-tape interpreter checks 1,722 program/input instances, program lengths 0–3, all legal jumps, input lengths 0–2 and bounds through 8 commits: 7,433 commits and 490 bounded nonhalting runs. All six MLX instruction/input examples remain checked. Audits are outside the modeled driver ledger. |
| NWMG pp.6–9 | `nwmg.py` implements the original thirteen-instruction blocks, 10/11 encoding and terminal M. 112 one-transition cases cover all six non-erasing rule forms, blank/existing neighbors and a multistate bilateral run through P/S. |
| ML-N's monotone metadata and decoder | `copier.py` generates an ordinary binary non-erasing transition table from any finite erasing binary table, including immediate halt. Its eight-bit blocks retain head/value triples, consumed and next-head flags, guards and active/retired handoff. `copier_checks.py` compares 28 erasing-TM cases to an independent interpreter, covering both moves, erasure, finite loop prefixes and copied window growth. One actual received-input erasure executes 738 non-erasing steps and 5,668 B commits through LAB and returns decoded `000`. This generator is generic; it contains no law answers or source-specific tables. |
| ML-T successor kernel with ML-R/Q | The one-state little-endian successor machine consumes an actually acquired bit at original sample/coordinate `(0,0)`. Its full copier/B/LAB run takes 2,516 non-erasing steps and 19,324 B commits, returning little-endian `01` from the actual copied tape. Original `LL`, `LL`, `LR`, `LL` calls retain four authentic receipts, the old index and the unchanged source port. This executes an arithmetic kernel, not the complete UJB evaluator. |
| ML-F | `framing.py` checks 511 V words and proper prefixes, 257 integers, six rationals, all 32,768 single-slot masks with two value patterns, and 181 framed host-law answers. A leading cylinder flag is part of the candidate's explicit grammar. These host-law answers are distinguished from compiled execution. |
| ML-A and UJB16.2–17.1 | `laws.py` checks exact old-array label DP, indexed conjunctions, rational priors 0, 1/7, 1/3, 1, partial masks at N=4,5,7, repeated/distinct-index labels, positive and null cases. The N=4 repeated complete-pair schedule visits all 225 partition pairs, prunes 224 and runs all 972/81 surviving G0/G1 type maps. A full square distinguishes positive G0/G1 laws and literal fourth-edge continuation. These are host exact-arithmetic checks. |
| MLX-capacity | Literal common-neighbor counting tests all 59,640 unordered left triples per source for K_(3,10), obtaining lift counts 0 and 92,575,499,208,192,000. N=73 distinct complete left labels gives zero evidence. Full Bell(30)^2/Bell(73)^2 schedules are not executed. |
| ML-R/Q | The provider privately traverses actual links of original P/S-prepared F1(z) objects. 8,248 calls cover four sources and 2,062 operands each, all four replies, long absences, slot suffixes and fixed marker suffixes. N=4 sessions complete sample 1, retain partial samples 2/4 and repeated receipts, leaving sample 3 unread. A null-law N=9 promised source still makes and receipts 136 actual calls. Static suppliers check the provider without supplying its returns. |
| ML-Cost | Work/routing, receipt/cache, framing, widening/old-new overlap, grants, cleanup and partial retained/peak coordinates are measured. Per-root censuses stay valid after later jobs in the same arena. The acquired successor session has a simultaneous **partial** peak of 1,320,626 bits, retaining 45,493 work/routing objects and four receipts. Complete installed code/control, framing/conversion/decoder counters/cursors/phases and the full evaluator job remain unvalidated; these partial figures cannot certify ML.16. |
| ML-V | Physical support/boundary round trips and annotated heads/control are checked. `(R,M)` on input `00` retains identical bare work geometry/routing while changing head and next action. No operational minimum or actual serial-history recovery is claimed. |

The content manifest names sealed lineage
`0434b0714812f29d6799f8e4ebd958f5eab2f987`; current candidate hashes,
not that commit's old theory bytes, bind the checks. AA/UJB/PG remain at
`c27e6de785aefbf37c0789a3f4cf63c6667ff3d7`, the exact structural supplier
at `865f38257ead1b4e857ff915c23b0b04276bd7f3`, static Read supplier at
`b490f4f17c97645ce7c6495ca8635bce3d3a1af7`, and unrestricted finite
indexed host evaluator at `bfd14367a5cbab3ee85cc0ce13db70dd797213b3`.
The NWMG primary PDF is https://arxiv.org/pdf/1304.0053v2; its copying
overview is distinguished from ML-N's explicit metadata construction.

After checkpointing unchanged candidate bytes, the caller can replay
`--check` using any checkout containing these same files and pinned Git
objects. No manifest change or result regeneration is needed merely
because HEAD changes. Changed theory content must receive a new manifest
and a fresh result; a saved result cannot establish an edited candidate.

The two unresolved obligations are **ML-N-evaluator-execution** (a complete
instantiated ML-T/W_U run returning ML.5, preserving every N/mask/prior
input domain) and **ML-Cost-full-evaluator-session** (that evaluator and
actual original-source continuation, with every installed and mutable
resource simultaneously charged). The strict flag preserves both.

These are exact simulations under paid ML-I premises, not hardware,
Lean verification, physical universality, a C117 frontend or independent
acceptance. Finite checks do not prove the universal ordinary theorems.
