# Fixed indexed-law evaluator constructor

`INTERFACE.json` seals the public program, table and complete encoded Wang word,
formats, input/output correspondence and paid-runtime dependencies. The owning
theory §§10–13 gives ordinary all-input proofs. This is no new Lean client,
new universality theorem, integrated LAB execution or complete resource account.

The fixed graph has 22,202 states and 361 Boolean-stack registers. It contains
only halt/goto/push/pop/peek. ML-F parsing, arbitrary finite integer/rational
arithmetic, vector access, partitions, old/new-array DP, capacities, all original
incidences, mixtures and output framing are expanded into those instructions.
The public incidence constants encode only the UJB15.2 source geometry.

The four-symbol table has 2,145,710 states. The fully materialized binary table
has 23,602,810 states and 47,205,620 rules (283,233,736 bytes). These tables are
losslessly gzip-packed and split into at most 32 MiB parts; the manifests hash
the original materialized bytes, complete packed stream and every part.

`W_U` has 496,162,138,109 instructions, represented by complete literal copier
templates and every resolved destination binding. Its complete compressed image
is 158,796,932 bytes, in five parts. The flat instruction array was not stored.
`constructor_image.Image.fetch` determines every instruction at every PC.
A consumer must pay for fetch/address search/decompression or actually install
and pay for a flat expansion. A Python image reader is no backed LAB controller.

Four full literal stack-reference runs returned independently checked laws:

| Input | Primitive steps | Received / output bits | Result |
| --- | ---: | ---: | --- |
| N4 full repeated pair; same and empty cylinders | 264,347,047 | 387 / 1308 | Z0=Z1=Z=1/240734712102912; ratios 1,0 |
| N4 fourth pair missing its last bit | 4,660,776,392 | 386 / 1560 | Z0=Z1=Z=281/61387351586242560; ratios 255/281,0 |
| N9 distinct complete left labels, one right label | 1,438,288,939 | 839 / 145 | All normalizers zero; NullEvidence |
| N4 full pair; unreduced prior 0/6 | 127,893,434 | 141 / 618 | Same positive normalizers; empty ratio list |

The runs execute the emitted graph with `constructor_reference.c`; expected
values are used only afterwards. `constructor-results.json` binds the input,
actual output and step census. They are full law evaluations on the stack
reference machine. Full binary/Wang/LAB law runs have not executed.

`audit-results.json` checks every one of the 23,602,810 resolved bindings,
every literal template transition, contiguous state bases, input encodings,
972 public incidence entries and the finite graph closure. It also includes
204 arithmetic/vector/framing probes, 180 compiler primitive instances and
16 copier/Wang image comparisons. `regression-results.json` records the DP
old/new-array mutation falsifier and portable invocation check. These finite
checks corroborate, and do not replace, arbitrary-word structural proofs.

Run with Python 3.12+, `-B`, assertions enabled, a C11 compiler for references,
and the repository cooperative worktree read/write scopes shown in the parent
README. All commands below are from the private experiment root, inside its
write scope. Supplier reads need `rro-structural-recovery` and
`rro-positive-certificates`. Use the actual candidate for the parent source
binding checks. Temporary unpacked files can be placed outside the repository.

```sh
python3 -B rro-lab-execution/constructor_checks.py --partial
python3 -B rro-lab-execution/constructor_compiler.py --output rro-lab-execution/constructor-tables
python3 -B rro-lab-execution/constructor_image.py --binary rro-lab-execution/constructor-tables/evaluator.binary.bin --output rro-lab-execution/constructor-code
python3 -B rro-lab-execution/constructor_pack.py --directory rro-lab-execution/constructor-tables --pack evaluator.four.bin evaluator.binary.bin
python3 -B rro-lab-execution/constructor_pack.py --directory rro-lab-execution/constructor-code --pack W_U.block.gz
python3 -B rro-lab-execution/constructor_audit.py
python3 -B rro-lab-execution/constructor_regressions.py
```

For sealed tables/images, use the unpacker without rebuilding:

```sh
python3 -B rro-lab-execution/constructor_pack.py --directory rro-lab-execution/constructor-tables --unpack evaluator.binary.bin --destination /tmp/rro-evaluator.binary.bin
python3 -B rro-lab-execution/constructor_pack.py --directory rro-lab-execution/constructor-code --unpack W_U.block.gz --destination /tmp/rro-W_U.block.gz
```

The image builder also leaves a duplicate `copier-bases.bin`; its exact bytes
are already contained in the block image and can be removed after verification.
The delivered image stores them once. Rebuilding tables/image requires space
for their actual unpacked sizes and packing temporaries. No full-run feasibility
or complete resource acceptance follows from these stored code sizes.

B consumes the fixed public objects and accepted received word. It must back
all code access, state, conversion, I/O and the unchanged original-port service.
C supplies the actual integrated full runs after that correspondence exists.
Both original strict acceptance failures remain in `exact-results.json`.
