# Finite positive catalogue of binary-relation properties

The target is the finite-domain extension left separate in Magnus Boman,
*A Catalogue of Properties of Binary Relations: Entailments, Incompatibilities,
and Independence Results*, arXiv:2609.13317v1, Remark 6.9 and Section 13.
This report supports a preregistration; it contains no Lean proof.

Run with Python 3, without dependencies:

```sh
python3 check.py > check.out
```

`check.py` uses only exact integer and Boolean operations. Two independently
implemented evaluators check the fifteen source definitions, allowing endpoint
interpolants in density and repeated successors in convergence. The program
reproduces the five finite maximal-consistent model profiles from the source,
all 188 profiles of the 66,066 relations on one through four labelled points,
and Proposition 6.8 on those relations. It checks the 54 proposed rules on
all these small relations.

For the proposed classification it evaluates 38 explicit matrices of size at
most nine, checks every maximal-omission list, and supplies actual model and
countermodel references for all 330 consistent closed masks. It checks all
32,768 antecedents and 491,520 antecedent/consequent pairs: 7,936 antecedents
are consistent with the rule system and 24,832 derive FALSE. The output ends
with `ALL_CHECKS_OK`.

The report also contains `symbolic-cover.json`, a compact propositional proof
forest for the fifteen omission cases. Its 123 excluded-middle splits, 137
profile leaves and one contradiction leaf are replayed by an independent
checker using only forward Horn steps and single-premise contraposition. The
maximum branch depth is nine. The interpretation of each node is valid for
an arbitrary rule-closed truth set; it does not rely on enumerating truth
assignments. A Lean implementation can replay it by cases and rule
applications, without `decide` on antecedent masks.

The finite calculation supports the proposed coverage and nine-point bound.
It does not prove universal rule soundness over every finite carrier or give
Lean kernel verification. The planned settlement must prove semantic rule
soundness and the general closure/model separation argument. Matrix checks
serve that argument; a positive exhaustive enumeration alone is insufficient.

## Data provenance

The original 51-rule manifest is from
`BinaryRelations/GlobalEntailmentCertificate.lean` in
[d3dd/binary-relations-catalogue](https://github.com/d3dd/binary-relations-catalogue),
revision `75278b521c476b55b39321c5809172dc8c0fc1ab`.
The three appended rules are Proposition 6.8:
P2+P3+P7 implies FALSE; P2+P7+P10 and P2+P7+P14 imply EMPTY.
Bit positions 0 through 14 denote P1 through P15; positions 15, 16, 17 and 18
denote EMPTY, NO_COMP, FALSE and UNIVERSAL.

Thirty matrices were independently selected by enumerating relations on at
most four points. The following larger matrices are copied from the pinned
source's `z3_witnesses/` and re-evaluated independently. These matrices,
including profile 8758, are existing ingredients, not claimed discoveries.

- Profile 182: `P3_P5_P8_not_P14_n5.json`.
- Profile 8758: `P3_P6_P10_not_P15_n8.json`.
- Profile 18998: `P3_P6_P10_P12_not_P14_n9.json`.
- Profile 19206: `P2_P3_P9_P10_not_P14_n6.json`.
- Profile 25102: `P2_P3_P4_P10_P15_not_P12_n6.json`.
- Profile 25142: `P3_P6_P10_not_P12_n8.json`.
- Profile 25350: `P2_P3_P9_P10_P15_not_P12_n5.json`.
- Profile 27190: `P6_P10_not_P13_n7.json`.

The immutable source archive is also available as
[Zenodo 10.5281/zenodo.22693298](https://doi.org/10.5281/zenodo.22693298).
A row integer's bit y records the edge from that row's point to point y.
The output includes all matrices, their exact profiles, maximal-omission lists
and a model/countermodel index for each consistent closed mask.

## Source code license

The copied rule manifest and larger adjacency matrices are covered by the
source's `LICENSE-CODE`, reproduced here:

MIT License

Copyright (c) 2026 Magnus Boman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
