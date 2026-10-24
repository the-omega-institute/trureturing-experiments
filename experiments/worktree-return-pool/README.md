# Worktree return and claim comparison

A single pair measured return plus claim at 34.184 seconds and remove plus create at 187.673 seconds. These are observed command times on a shared host, with the initial cache provisioning and both later cache/build phases reported separately.

Both subjects started at `fff23e51861ff10d2f991e911056ad6971361896`. Both requested that exact base, used `leanprover/lean4:v4.35.0`, had equal manifest bytes (SHA256 `8550e44937f85a9be196a2917a075eeb8c22bf6c15baa83e72017c6fc4437a8e`), and completed canonical `make lean-cache-ensure` before comparison. Each initial `.lake` contained 334,292 files. The return arm synchronized to freshly confirmed dev before claiming the requested base; the recreate arm created directly at that base. Both finished at the requested commit.

The authorized delivered subjects were PR #21147 for return and PR #21160 for normalization/removal. The second subject was normalized through canonical removal and creation, preserving its delivered task ref. The program verifies MERGED state and original branch identity before normalization. After comparison it explicitly removes only the two named subjects and preserves their task/recovery refs. No other tree is a measurement subject.

| Phase | Wall seconds | Files created | Files modified | Files removed | Available KiB delta | stdout Built count | olean files written |
|---|---:|---:|---:|---:|---:|---:|---:|
| return | 21.889 | 555 | 29 | 522 | -286008 | 0 | 0 |
| claim (one winner, two concurrent requests) | 12.295 | 522 | 44 | 555 | -841764 | n/a | 0 |
| reused-cache-door | 4.995 | 0 | 26 | 0 | -297912 | 0 | 0 |
| reused-build | 11.239 | 0 | 0 | 0 | -40676 | 0 | 0 |
| remove | 120.493 | 0 | 0 | 665169 | 1747456 | 0 | 0 |
| create | 67.180 | 330747 | 0 | 0 | -3077004 | 0 | 0 |
| recreated-cache-door | 42.040 | 334422 | 15 | 0 | -1542808 | 0 | 21376 |
| recreated-build | 9.628 | 0 | 0 | 0 | 624 | 0 | 0 |

Counts compare entry inventories immediately around each command, including ignored caches and the tree's Git pointer, excluding Git directories and metadata outside the tree. Modified means size, mtime or inode changed. These counts measure changed entries, not individual filesystem calls; entries can appear in several phases. Inventory collection is outside command wall time. Available-space readings are raw `df -k` observations on the shared filesystem; concurrent host activity contributes, so the deltas do not establish causal physical savings.

| Cache door | Status | Method | mathlib olean state | Project olean state |
|---|---|---|---|---|
| initial-cache-reuse | seeded | clonefile | warm | warm |
| initial-cache-recreate | seeded | clonefile | warm | warm |
| reused-cache-door | present | none | warm | warm |
| recreated-cache-door | seeded | clonefile | warm | warm |

Both canonical `make lean LEAN_TARGETS=D5.S0.Carrier.TraceConjugation LEAN_SKIP_LOCK=1` commands exited 0. The same documented global-lock option applies to both build checks; their cache pin/freshness and per-tree operation doors remain active. The interrupted default-lock sample records exit 2 and 469.509 seconds of command wall time while waiting for the shared lock; it is excluded from build verification and timing comparison. The default-lock reuse sample is also retained separately in the raw data. The wrapper uses error-level Lake logging; zero stdout Built lines alone is not a build count. The olean entry inventories provide the accompanying observed write count. Raw cache receipts and build outputs are in `data.json`.

Two exact-path claims started through a barrier with different task names. Exactly one succeeded; the other returned `WORKTREE_RETAINED busy_scope:operation:idle-pool`. The winning claim's wall time is the claim phase above. This trial establishes the observed exclusion behavior for this pair of requests.

`implementation.json` binds the Python pool source and built CLI binary by SHA256, captured before measured return/claim. Their digests match source/binary at main delivery commit `7bb4719417553728e26601b46c9981d59debd168`. `data.json` separately records the checkout HEAD observed at driver startup. The measured pool source and CLI binary match the delivery commit above.

This is one pair and one concurrent trial. The ordered phases, shared filesystem traffic, warm donor and cache inventory overhead limit timing generalization. Neither a population speedup estimate nor a general physical-space saving is claimed.

Run from the experiments repository root with two explicitly authorized delivered subjects:

```sh
python3 experiments/worktree-return-pool/compare.py \
  --implementation /Users/auric/trureturing-219aba10-21130 \
  --reuse-tree /Users/auric/trureturing-219aba10-21131 \
  --recreate-tree /Users/auric/trureturing-219aba10-21129 \
  --return-pr 21147 --recreate-pr 21160 \
  --module D5.S0.Carrier.TraceConjugation \
  --output experiments/worktree-return-pool/data.json
```

The corresponding build completion command is:

```sh
python3 experiments/worktree-return-pool/compare.py \
  --implementation /Users/auric/trureturing-219aba10-21130 \
  --reuse-tree /Users/auric/trureturing-219aba10-21131 \
  --recreate-tree /Users/auric/trureturing-219aba10-21129 \
  --return-pr 21147 --recreate-pr 21160 \
  --module D5.S0.Carrier.TraceConjugation \
  --output experiments/worktree-return-pool/data.json \
  --resume-builds --skip-build-lock
```

For a fresh comparison with the reported build configuration, use `--skip-build-lock` on the first command. Resumption verifies both subject commits, clone identity, toolchain and manifest before building and preserves the preceding samples.

The program records every exact phase command, cwd, exit code, command wall time, cache receipt and filesystem reading in `data.json`. The built implementation CLI dispatches return, claim, remove and create; cache provisioning and builds use the subjects' canonical make entrances.
