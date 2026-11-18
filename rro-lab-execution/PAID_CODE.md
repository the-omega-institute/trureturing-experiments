# Installed-image decoder component

`paid_code.py` reads the complete uncompressed constructor-A block image with
a fixed register graph. It performs bit scans for arithmetic and restoring
division, follows actual page-directory references, resolves actual base and
binding records, and sends a 67-bit instruction reply through Print/End/Ack/
Return phases. `Decoder.result()` is an independent observer view; it must
not become a free instruction source for the LAB runtime.

`paid_code_checks.py --full` prepares A's unchanged compressed artifact
offline, installs all uncompressed image pages, and compares selected actual
fetches against A's reference reader. It retains the code image's hash and
separately counts pages, directory links, register fields, code nodes and
grants. It checks repeated use without allocation, pending reply states,
audit erasure and an independent primitive-operation census. The mutation
`--mutation omit-register-reads` must fail with
`[FAIL] register_bit_reads_census`.

This is a register-level component model. It does not supply a complete
installed ML-I primitive interpreter or its installation/selection scratch
and peak account. The high-level register graph is present in
`paid-code/decoder.graph.json`; the finite engine executing that graph still
needs its literal primitive realization and tariff correspondence. The
initializer's bulk field tariff is not a chronological installer execution.

No original-port, ML-F parser/loader, actual-tape output decoder or complete
LAB controller is replaced by this component. The original `lab.py`,
`framing.py` and `session.py` still have operational host values identified
in the dispatch inventory. Their obligations remain. In particular,
`records[-1]` and `Print.emitted` remain operational in that old path.
There is no complete shared source/observer account, widening-overlap
controller, full evaluator run, or B/C release.

Run from any working directory, using Python 3.12+:

```sh
python3 -B <experiment-root>/rro-lab-execution/paid_code_checks.py --full --output <result-path>
python3 -B <experiment-root>/rro-lab-execution/paid_code_checks.py --full --check <result-path>
```

Temporary uncompressed image files and reference-reader buffers are outside
the modeled fetch controller. They are fixed public installation inputs or
independent test data; deleting them after the pages have been installed
does not change the modeled fetch. Runtime decompression has not been
implemented or charged. A consumer choosing runtime decompression must
provide that additional state/operation correspondence.

The component's retained storage is not a full-session peak. The exact
acceptance suite retains both original failure obligations and a false
acceptance result. The long research goal and all original input domains
remain unchanged.
