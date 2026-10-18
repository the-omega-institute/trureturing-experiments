# Native complete-continuation certificate data

These datasets accompany the main repository's [connected-interval theory](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/FIB_NATIVE_COMPLETE_CONTINUATION_CONNECTED_INTERVAL.md) and [two-parameter-strip theory](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/FIB_NATIVE_COMPLETE_CONTINUATION_TWO_PARAMETER_STRIP.md). Each payload file contains the original body bytes of one machine-readable code fence, including its final newline. Mathematical arguments and definitions belong to the theory text. The data provide rational interval enclosures and finite certificate inputs; they are not Lean certificates or kernel-verified proofs. Decimal display fields retain their original diagnostic status.

| Dataset | Source sections | Files | Payload bytes | Scope |
| --- | --- | ---: | ---: | --- |
| [connected-interval](connected-interval/manifest.json) | 12–16 | 10 | 36,837,685 | Section 13: 209 uniform cells, 1,672 arcs, 5,333 covering leaves; the remaining files include the source enclosure, event box and embedding margins, three seven-band joins, endpoints, seam response and feasibility data. |
| [two-parameter-strip](two-parameter-strip/manifest.json) | Appendices A–D | 67 | 5,022,436 | 50 uniform cells contain 400 arcs and 1,750 leaves. Two uniform edge boxes and two separate midpoint boxes add four roots, 32 arcs and 120 leaves: all 54 root boxes contain 432 arcs and 1,870 leaves. Independent enclosures and direct-evaluation summaries remain separate data. |

Each linked `manifest.json` gives every payload's fixed filename, source section and fence position, byte length, record count and complete SHA-256. A JSON object counts as one record; a top-level JSON array counts its elements; JSONL counts nonempty lines. The manifest's `root_boxes`, `arcs` and `leaves` count explicit objects and arrays, excluding scalar summary counts to avoid counting those data again.

The connected-interval source uses `text` fence labels, and the strip source uses `textJSON`; payload formats are determined by parsing their bodies. Connected-interval Section 13 has three JSONL files: 203 six-band cells, two final six-band cells and four seven-band cells. Its other seven files preserve all event, join, endpoint, response and source data.

The strip's 50 cell filenames retain the `final-strip/cell-NNNN.json` selectors used in independent-enclosure records. Its four edge and midpoint filenames retain their original selectors. The parameter-cover table uses cell basenames within `final-strip/`. Other `file`, `edge` and numeric `cell` fields identify these mathematical records and do not execute external programs.

Run the integrity and JSON-readability check from this repository root:

```sh
python3 experiments/theory/fib-native-complete-continuation/verify_data.py
```

To additionally compare every payload with its exact source fence in a local main-repository Git checkout:

```sh
python3 experiments/theory/fib-native-complete-continuation/verify_data.py --source-repo ../trureturing
```

The optional check reads the pinned source commit from Git, so it also works when the current theory text uses external data links. These checks establish byte integrity, parsability and stated finite record counts; they do not evaluate interval inequalities or recompute the experiments.

The payload data originate in the main repository under its Apache-2.0 license; the same license applies here through the [repository license](../../../LICENSE). The payloads include no third-party source code.
