# Licenses and source attribution

The original experimental programs/data and the retention code/documentation
are contributions under this private repository's Apache License 2.0; see
`../../LICENSE`. Source-derived statements and mathematical paraphrases are
attributed separately below and retain the source's CC BY 4.0 license.

Alexandr Kostochka et al., *Turán number of three disjoint cliques of arbitrary
sizes*, arXiv:2610.38430v1 (2026), Section 11, Conjecture C1:
<https://arxiv.org/abs/2610.38430v1> and
<https://arxiv.org/html/2610.38430v1#ThmConjecture1>.
License: Creative Commons Attribution 4.0 International,
<https://creativecommons.org/licenses/by/4.0/>.
The full author list is identified by that source citation. No endorsement by
the source authors is implied.

The source version/license and attribution were supplied in the authorized
input and prior documentary source metadata; the retention replay did not
fetch or independently reverify source/license pages. No full article or source
snapshot is redistributed. Changes consist of finite-graph notation,
paraphrase, encoding explanation, and independent replay descriptions.

The original `packing_core.py`, `packing-witnesses.json`, and
`clique-results.json` are copied without byte changes. `capacity_packing.py`
extracts the existing capacity-coloring function, removing orchestration and
wall-time wrappers and replacing its limit with NodeLimit. `encoding.md` is
rewritten from the existing explanation, explicitly clarifying that isolate
filling uses vertices already in the input graph. `fixed-inputs.jsonl.gz`
materializes only the original fixed stream. The standalone verifier, current
validation report, README, and checksums accompany these reusable artifacts.
