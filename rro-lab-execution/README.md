# Sealed monotone LAB execution: exact finite checks

This is assignment2 experimental evidence against theory commit
`0434b0714812f29d6799f8e4ebd958f5eab2f987`, not acceptance of the draft.
The implemented checks pass. **The complete compiled UJB evaluator and its
full simultaneous resource account remain unexecuted.** `--require-all`
retains this acceptance failure with exit 2. No finite counterexample to a
theorem was found in the implemented scopes; an unexecuted interface is not
reported as a disproved theorem. Unit3 must reconcile these obligations.

Python 3.12+, standard library, assertions enabled, and local immutable Git
objects are required. Run from this checkout through its cooperative scopes:

```sh
python3 -B /Users/auricstudio/trureturing/tools/scripts/worktree/worktree_protocol.py \
  --source /Users/auricstudio/trureturing-experiments with --path "$PWD" \
  --read rro-structural-recovery --read rro-positive-certificates \
  --write rro-lab-execution -- \
  python3 -B rro-lab-execution/exact_checks.py \
  --output rro-lab-execution/exact-results.json
```

For a deterministic comparison use `--read rro-lab-execution` instead of its
write scope and `--check rro-lab-execution/exact-results.json` instead of
`--output`. Add `--require-all` to enforce all approved execution obligations.
The default exit 0 means the implemented checks passed; it does not override
`acceptance_complete: false` in the result. No randomness or floating point
is used. The result contains all source and program SHA256 values.

| Sealed source | Program correspondence and exact tested scope |
| --- | --- |
| ML-I/C/L/B, §§1–3 | `lab.py` uses bound object references, physical terminal/unary/binary nodes, guarded P/S, saved heads, two comb arms, real ancestor headers/cells and edit/append pending phases. An independent integer-tape interpreter checks 1,722 program/input instances: lengths 0–3, all legal jumps, all binary inputs of length 0–2, and bounds through 8 commits. There are 7,433 commits and 490 nonhalting bounded runs. All six MLX instruction/input examples are checked, including a 65-commit loop prefix. |
| NWMG pp.6–9, ML-N's supplied compiler stage | `nwmg.py` implements the original 13-instruction blocks, 10/11 encoding and terminal M. 112 one-transition fixtures cover all six rule forms and blank/existing neighbors; a four-transition bilateral run also executes through physical P/S. This does not instantiate the ML-T copying machine or the full evaluator's `W_U`. |
| ML-F, §§4.1–4.2 | `framing.py` implements U/V/List, signed rationals, original-index masks, cylinders, service replies and Receive/Print/End/Ack/Return. It checks 511 V words and all proper prefixes, 257 integers, six rationals, all 32,768 single-slot masks with two value patterns, and 181 framed host-arithmetic replies. Cylinder conflict flags use a leading bit, a layout left unspecified by the prose. |
| ML-A, ML.6–7, §8.3; UJB16.2–17.1 | `laws.py` independently implements old-array label DP and indexed conjunctions. For N=4 repeated complete pairs it runs all 225 partition pairs, prunes 224, then enumerates all 972/81 surviving G0/G1 type maps, keeping all four incidence occurrences and the denominator `1944**4`. Other cases reuse the old inversion supplier only at at most four informative indices. Unread indices integrate to one under the same shared injections. N=4,5,7, represented priors 0,1/7,1/3,1, all 15 coordinates, same/contrary/repeated/distinct-index restrictions, and null ratios are checked. A full N=4 square also distinguishes the positive G0/G1 laws and checks its conditional fourth-edge continuation against literal common-neighbor counts. |
| MLX-capacity | Independent literal common-neighbor counting tests all 59,640 unordered left triples per source for K_(3,10), obtaining lift counts 0 and 92,575,499,208,192,000. N=73 distinct complete left labels gives zero evidence. These are exact formula checks, not executions of the enormous full partition schedules. |
| ML-R/Q, §§5,8.4 | `session.py` prepares original F1(z) objects using the original P/S preparation, preserving its routing overhead. The private provider traverses physical links at the same bound root; the session obtains only actual original replies. 8,248 provider calls cover four fixtures and 2,062 distinct operands each, including all four replies, all slot suffixes through length 4, marker suffixes through length 5, all root words through length 7 and long absent operands. Two N=4 sessions keep every receipt, complete sample 1, reread LL three times, and retain partial samples 2 and 4 with sample 3 unread. A separate N=9 legal promised source forces K_(9,1), has zero mass under both laws, and still performs 136 authentic calls and retains all indices. |
| ML-Cost, §§6,8.5 | Primitive counters and object/field censuses check ML.10–12 and ML.15, receipt/cache sizes, the two-link head bound, ancestor visits, actual reference widening/overlap, buffers, linked free storage, pending/unlinked receipts, call framing, output, acknowledgement and cleanup for implemented runs. Provider primitive charges and inherited preparation are separate; provider observations are diagnostic counts, not an extra summand. Audits are excluded from the modeled driver ledger. Fixed driver/code installation, all parser/loop counter backing, the full arithmetic job and their simultaneous peak are not instantiated, so these are partial resource measurements. |
| ML-V, §7; AA17.3 | Actual physical support and terminal boundary round trips, comb parsing and annotated head/control comparisons are checked. The reachable `(R,M)` on input `00` retains identical bare work tree/routing but changes head and next action. No reverse minimality or actual serial-history claim is tested. |

The theory file SHA256 is
`0f93b9f2634965902d7ea35ee675e1886d2309ec10b017e926ff756f3136fcd2`.
`pins.py` reads fixed Git blobs, never current theory worktree bytes. It pins
AA/UJB/PG at `c27e6de785aefbf37c0789a3f4cf63c6667ff3d7`, the old exact
supplier at `865f38257ead1b4e857ff915c23b0b04276bd7f3`, the static Read
oracle at `b490f4f17c97645ce7c6495ca8635bce3d3a1af7`, and the general
finite indexed evaluator at experiment base
`bfd14367a5cbab3ee85cc0ce13db70dd797213b3`. The static oracle checks the
new physical provider; it is not used to supply its returns. The intrinsic
recovery supplier is preserved; its earlier geometry checks are not rerun.

These are exact simulations of the conditional installation, not physical
hardware, private-source access, Lean verification, or exhaustive tests of
infinite domains. Host rational evaluation and printing do not establish
ML-A → ML-T → copier → B → LAB evaluation. Tape steps, arithmetic operations,
bit widths, material and navigation are recorded on their own scales.

All files are confined to `rro-lab-execution/`. This path is outside the
existing sparse cone and needs caller-owned sparse-checkout adjustment.
No branch, index, commit, push, PR, source theory, supplier or unrelated cone
is changed by these programs.
