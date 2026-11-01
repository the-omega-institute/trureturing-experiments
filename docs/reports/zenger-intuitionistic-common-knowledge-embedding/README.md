# Intuitionistic common-knowledge embedding: exact finite evidence

`python3 check.py` checks the proposed translation F(φ) = (C Dφ) → φ, where Dφ is the conjunction of ψ ∨ ¬ψ over ordinary syntactic subformulas ψ of φ. The guard is prior art: Santiago-Fernández, Fernández-Duque and Joosten, *The Complexity of the Constructive Master Modality*, [arXiv:2603.05131v1](https://arxiv.org/abs/2603.05131v1), Definition 6.1 and Theorem 6.5 (one unrestricted relation). The proposed extension concerns simultaneous agent relations, their union closure, and preservation of individual reflexive/S4 conditions.

The source is Lukas Zenger, *Intuitionistic Common Knowledge*, [arXiv:2605.00523v1](https://arxiv.org/abs/2605.00523v1), Sections 2 and 4. Worlds form a nonempty poset; atomic valuations persist; each agent satisfies w ≤ v and v R_a u implies w R_a u. Common knowledge uses the reflexive transitive closure of the union of all agent relations.

The program enumerates every labelled poset, every pair of admissible agent relations, and every pair of persistent atom valuations on one and two worlds. It checks 91 formulas (all syntax trees of size at most three over two atoms and bottom, plus the explicitly listed nested/modal cases), comprising 104 distinct subformulas. All calculations use Python integers, sets, and booleans; no floating point, sampling, or external packages are used. `check.out` is the complete output.

- Source fidelity: truth persistence (Section 2's Monotonicity lemma), and the displayed valid fixed-point law Cφ ↔ (φ ∧ E Cφ), checked throughout the finite domain.
- Candidate: agreement with classical truth on the guarded generated region; classical equivalence of Fφ and φ; transfer of failed translated implications; identity-order countermodel transfer; preservation of each relation's reflexivity/transitivity under restriction; and the quadratic size bound on the tested inputs.
- Model counts: 16 one-world models and 5,554 two-world models; the latter include 328 T models and 328 S4 models. The coincidence of T and S4 counts is specific to two-world relations. There are 1,012,284 formula/root checks and 5,295,348 guarded subformula comparisons, all successful.

These are finite semantic checks. They neither prove universal validity equivalence nor establish a complexity lower bound or literature novelty. The planned settlement uses structural induction for arbitrary formulas, arbitrary nonempty world sets, and every finite nonempty agent set. No Lean proof is supplied by this report.
