# Sources and reuse

Svante Linusson and Petter Restadh, *Rhombus Criterion and the Chordal Graph
Polytope*, arXiv:2305.05275v1:

- Original PDF: https://arxiv.org/pdf/2305.05275v1
- Original source: https://arxiv.org/src/2305.05275v1
- Conjecture 3.4 is on PDF page 5. The governing event uses every positive n,
  all labelled simple chordal graphs, every clique coordinate of size at least
  two, the ordinary real convex hull and ordinary one-dimensional faces.
- Source Delta definition and Lemmas 3.8, 3.9 and 3.11 govern component
  compatibility. Proposition 3.12, Lemma 3.15 and Theorem 3.16 supply the proper
  known cases skipped by this route.
- The acquired source archive has SHA256
  `500e2f40424382858e3eb2fe87e268216bc0da92bc7187e5f5722c31ab78e089`.

This implementation reuses these mathematical definitions and reductions. It is
originally authored Python; no upstream program, article text, diagram, data
set or solver implementation is copied. The author's finite computation
repository, https://github.com/PetterRestadh/rhombus-criterion at revision
`e9b46f9dfc28810f559afc875f78714afc11a329`, was qualified by the permitted prior
source conclusion. This program enforces endpoint exclusion and checks all
coordinates directly. It makes no claim about the correctness of historical
computations, global originality, priority or complete literature coverage.

The completed qualification's immutable native API catalog has SHA256
`18097f3ef592987461f04f1f8d70a02caca7e65d56e1ccc747f359b50629519a`.
It supplies only source-level bridge plans. No Lean declarations or dependencies
from that catalog are copied or compiled here. The canonical main-repository
observer remains read-only.

# Licensing

The authored Python program is licensed under Apache-2.0; the full license is
in `LICENSE`. The authored explanatory Markdown is licensed CC-BY-4.0
(https://creativecommons.org/licenses/by/4.0/). The generated computational
data is dedicated under CC0-1.0
(https://creativecommons.org/publicdomain/zero/1.0/). Cited works retain their
own copyrights and licenses; they are references, not redistributed material.
