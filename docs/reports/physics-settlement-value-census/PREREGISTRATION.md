# Fixed census reproduction and source audit

This registration fixes E0/E1 before execution of the analyzer or the label
audit. It does not register a value function, a global research policy, or
an automatic quota. E2 prospective trials and E3 criterion backtests are
outside this experiment.

## Bound inputs

- Main-repository observation snapshot:
  `53ec26edbed43dc4f4862d2db303ad80cb3dbd92`.
- Experiment input commit:
  `c42809bac1b923b5fc78aaa669125516c3a65ac7`, retained by base
  `9567f90e9a8333f0d75d7a04c7e37a57f33dea9c`.
- Fixed input: `physics-settlements-36.batch.json`, SHA256
  `8a26793329161a9b0654d5976e187c8e2d167b72cebd273fab57e51d6d90004d`.
- Candidate analyzer commit:
  `a122da84cd8e7966b67057aef266417aee3ed2df`, read and run only in
  `/Users/auricstudio/trureturing-0addbbde-l3`. Its tracked files must not
  change. Build failure is reported, not bypassed with another implementation.
- Original labels, requests, batch, Python preview and preview output are
  preserved. Audited/corrected data receive separate versioned files.

## E1-A success, falsification and stop

Success requires the actual analyzer batch entry point to exit zero, and
each order to report 36 states, 1,260 ordered distinct-state pairs, 102
finally inseparable ordered pairs and ten task-conflict ordered pairs.
The conflict identities must be exactly the two orientations of
{15354,16702}, {16702,16902}, {16702,17145}, {16702,17605},
{16971,17778}. `residual_sufficiency.sufficient` must be false.
Report the actual `determined_at_layer`, preserving JSON null semantics.

Compare the actual output with `replica.py` and `replica-preview.txt` at
every shared definition: readout identities and kernels, individual pair
identities, unique captures, leave-one-out readouts, ordered layer readings,
and task-on-fiber readings. Aggregate agreement alone is insufficient.
An unavailable preview field is explicitly distinguished from a mismatch.

Any mismatched count, pair identity or null meaning falsifies the registered
reproduction expectation. Locate whether the discrepancy belongs to the
input, the Python preview, or the analyzer; do not silently change the
expectation. Stop after one reconciliation of these fixed inputs. Do not
search for an input or analyzer that yields the registered counts.

## E1-B source audit

For each of the 36 settlements, keep the original self-reported label in a
separate column from the checked value and its status (verified, missing,
or inferred). Cite PR, commit and/or file location for proof shape,
admission basis, conclusion type, utility, search track, merge period and
the task: existence of an external D5 importer of an added D5 module.
Check module/header-derived fields against the fixed observation snapshot;
record the PR self-report as such rather than treating it as kernel proof.

Locate the proof-shape inference already mentioned in the census and give
its actual evidence. An unlabelled search track is missing information,
not a verified categorical track. Missing values must not be silently
imputed. Preserve a distinct status channel when materializing equality
readouts. Any inspection of live proof-path use would be a separate field;
it cannot replace the registered importer task.

Retain corrected labels and requests separately and recompute their actual
counts; 102/10 is not an expectation for corrected data. Report missing,
inferred, changed and checked rows individually. Stop after the fixed
36-row audit, leaving unavailable evidence explicit.

## Delivery and limits

Retain only formal findings, necessary programs/data and traceable source
references. Do not commit process logs, command transcripts or runner
receipts. This registration is a separate ancestor commit of result
commits. Push the experiment branch and open a PR to this repository's
`main`; do not merge it. Only the designated experiment worktree is a
write target, apart from the explicitly required runner result artifacts
and the authorized analyzer build products. No main-repository branch,
tracked file or PR is changed.
