# Physics settlement value census

A finite census of 36 merged physics settlement pull requests of
[trureturing](https://github.com/the-omega-institute/trureturing), created
2026-10-01 to 2026-10-16, presented as a typed-source arena for the main
repository's AURIC FIB analyzer (`tools/scripts/auric-fib-analysis.sh`). Each
settlement is one state. Six readouts and one task are equality kernels on
these states. The census asks which readouts separate which settlements and
whether the readouts determine the task.

Readings are taken at main-repository commit
`53ec26edbed43dc4f4862d2db303ad80cb3dbd92`.

## Population

All 39 merged physics settlement PRs created in the window, less #12598
(Bell polynomials, not physics) and #17345, #17379 (no `proof_shape` fields).

## Readouts and task

| Name | Meaning | Source |
| --- | --- | --- |
| `shape` | `proof_shape` of the settling theorem (`bind-only` / `content`) | PR self-report (`labels.json`) |
| `admission` | module `admission_basis` | PR self-report |
| `resolution` | `Proved` / `Refuted` / `none (library result)` | PR body |
| `utility` | strongest `utility:` header among D5 modules added by the merge (`certified-instance` > `checker` > `none`) | module sources at the commit above |
| `track` | search track named in the PR provenance (`quick` / `deep` / `unlabelled`) | PR body |
| `period` | merge date bucket, early 2026-10-01..09 or late 2026-10-10..16; exposure control | merge date |
| `downstream` (task) | some D5 module not added by the same merge imports an added module | import graph at the commit above |

## Files

| File | Content |
| --- | --- |
| `labels.json` | curated per-PR readouts with PR URLs and merge commits |
| `build_census.py` | computes `utility`, `period` and `downstream` from a main-repository checkout and writes the requests |
| `physics-settlements-36.request.json` | one typed-source request; readouts added in the order shape, admission, resolution, utility, track, period |
| `physics-settlements-36.batch.json` | the same request and a second order with `period` first |
| `replica.py` | Python preview of the analyzer's typed-source input checks and counts, plus a task-on-fiber check |
| `replica-preview.txt` | preview output for the batch |

## Reproduction

```sh
python3 build_census.py --labels labels.json --main-repo PATH/TO/trureturing --out .
python3 replica.py physics-settlements-36.batch.json
```

The main-repository checkout must contain every listed merge commit and its
first parent; in a shallow clone fetch each reported commit with
`git fetch --depth=2 origin <sha>`. The canonical counts come from the
analyzer in a main-repository checkout with the Lean toolchain:

```sh
make -C tools auric-fib-analysis-build
tools/scripts/auric-fib-analysis.sh --batch < physics-settlements-36.batch.json
```

## Preview readings

From `replica-preview.txt`; the analyzer itself has not been run on these
requests.

- Under all six readouts, 102 of 1,260 ordered settlement pairs remain
  inseparable.
- Unique capture: `track` 66, `shape` 24, `period` 12, `admission` 10,
  `resolution` 0, `utility` 0. `resolution` and `utility` are redundant given
  each other: every refutation carries a `certified-instance` header and no
  proof does. Adding `utility` after `resolution` is a collapsed layer.
- The readouts do not determine the task. Five unresolved pairs differ in
  `downstream`, in two fibers: {15354, 16702, 16902, 17145, 17605}, where only
  16702 is imported, and {16971, 17778}, where 17778 is imported by seven
  modules.
- The largest fiber has nine settlements: quick-track, bind-only,
  certified-instance refutations, all without importers. No current readout
  separates them.
- Two bind-only settlements are imported: 12833 by five modules and 17778 by
  seven.
- `downstream` holds for 6 of 14 content and 4 of 22 bind-only settlements,
  and for 6 of 10 early and 4 of 26 late merges. Six of the ten early merges
  are content settlements, so these 36 states do not separate shape from
  exposure time.

## Boundaries

- `shape` and `admission` are PR self-reports and were not re-derived from
  the Lean sources. The settling theorem of #15993 has no `proof_shape` row in
  its PR body and is recorded as bind-only. The resolution of #14159 is read
  from the PR body answer that the bound is tight for every n.
- An import is not a use on a live proof path, and `downstream` depends on
  time since merge.
- `track` is unlabelled for eight settlements.
- The analyzer reports counts over the supplied table rows; it does not
  certify that the rows represent the settlements.

## Settlements

| PR | merged | shape | admission | resolution | utility | track | period | downstream | external importers |
|---|---|---|---|---|---|---|---|---|---|
| 11696 | 2026-10-01 | content | open-problem-resolution only | Proved | none | unlabelled | early (2026-10-01..09) | no external D5 importer | 0 |
| 12188 | 2026-10-02 | content | escape-witness + open-problem-resolution | Proved | none | unlabelled | early (2026-10-01..09) | has external D5 importer | 1 |
| 12833 | 2026-10-04 | bind-only | open-problem-resolution only | Proved | none | unlabelled | early (2026-10-01..09) | has external D5 importer | 5 |
| 13018 | 2026-10-04 | content | escape-witness only | none (library result) | none | unlabelled | early (2026-10-01..09) | has external D5 importer | 3 |
| 13245 | 2026-10-05 | content | escape-witness + open-problem-resolution | Proved | none | unlabelled | early (2026-10-01..09) | has external D5 importer | 2 |
| 13409 | 2026-10-05 | bind-only | open-problem-resolution only | Refuted | certified-instance | deep | early (2026-10-01..09) | no external D5 importer | 0 |
| 13680 | 2026-10-06 | bind-only | open-problem-resolution only | Refuted | certified-instance | unlabelled | early (2026-10-01..09) | has external D5 importer | 1 |
| 14064 | 2026-10-07 | content | escape-witness only | Proved | checker | deep | early (2026-10-01..09) | no external D5 importer | 0 |
| 14159 | 2026-10-07 | content | open-problem-resolution only | Proved | none | deep | early (2026-10-01..09) | has external D5 importer | 1 |
| 14605 | 2026-10-09 | bind-only | open-problem-resolution only | Proved | none | quick | early (2026-10-01..09) | no external D5 importer | 0 |
| 15193 | 2026-10-10 | content | open-problem-resolution only | Proved | none | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 15354 | 2026-10-11 | bind-only | open-problem-resolution only | Proved | none | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 15530 | 2026-10-11 | content | escape-witness + open-problem-resolution | Refuted | certified-instance | quick | late (2026-10-10..16) | has external D5 importer | 3 |
| 15565 | 2026-10-11 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 15993 | 2026-10-12 | bind-only | open-problem-resolution only | Refuted | certified-instance | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 16103 | 2026-10-13 | content | escape-witness + open-problem-resolution | Refuted | certified-instance | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 16361 | 2026-10-13 | content | escape-witness only | none (library result) | none | unlabelled | late (2026-10-10..16) | has external D5 importer | 2 |
| 16702 | 2026-10-14 | bind-only | open-problem-resolution only | Proved | none | deep | late (2026-10-10..16) | has external D5 importer | 1 |
| 16709 | 2026-10-14 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 16727 | 2026-10-14 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 16879 | 2026-10-16 | bind-only | open-problem-resolution only | Proved | none | unlabelled | late (2026-10-10..16) | no external D5 importer | 0 |
| 16887 | 2026-10-14 | content | open-problem-resolution only | Proved | none | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 16902 | 2026-10-14 | bind-only | open-problem-resolution only | Proved | none | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 16949 | 2026-10-15 | content | open-problem-resolution only | Proved | none | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 16971 | 2026-10-14 | bind-only | open-problem-resolution only | Proved | none | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17103 | 2026-10-16 | content | open-problem-resolution only | Proved | none | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 17145 | 2026-10-15 | bind-only | open-problem-resolution only | Proved | none | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 17178 | 2026-10-15 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17487 | 2026-10-15 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17489 | 2026-10-15 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17605 | 2026-10-16 | bind-only | open-problem-resolution only | Proved | none | deep | late (2026-10-10..16) | no external D5 importer | 0 |
| 17628 | 2026-10-16 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17653 | 2026-10-16 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17760 | 2026-10-16 | bind-only | open-problem-resolution only | Refuted | certified-instance | quick | late (2026-10-10..16) | no external D5 importer | 0 |
| 17778 | 2026-10-16 | bind-only | open-problem-resolution only | Proved | none | quick | late (2026-10-10..16) | has external D5 importer | 7 |
| 17791 | 2026-10-16 | content | escape-witness + open-problem-resolution | Proved | none | quick | late (2026-10-10..16) | no external D5 importer | 0 |
