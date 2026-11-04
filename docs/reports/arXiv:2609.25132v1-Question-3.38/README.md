# Question 3.38: separate definable and external Ramsey degrees

Source: K. Krupiński, J. Lee and S. Moconja, *Ramsey theory and topological dynamics of 0-dimensional flows*, [arXiv:2609.25132v1](https://arxiv.org/abs/2609.25132v1), Definitions 3.1, 3.15–3.16 and Question 3.38.

Run `python3 check.py`. The program uses only integers, tuples and finite sets. Its output is in `check.out`; success means exit 0 and final line `ALL_CHECKS_OK`.

The proposed refutation uses the full integer shift over an alphabet of size r ≥ 2, with a dense-orbit basepoint and the fixed singleton clopen family Σ = {x : x(0) = 1}. Orbit parameters can approach the zero fixed point simultaneously. For m ≥ 2, m distinct periodic external parameters yield exactly m colors on every translate of {0,…,m−1}.

The exact checks cover:

- Fidelity: finite instances of the remark following Definition 3.16, that every ambit has separately finite weakly definable degree. All 79,170 representations on regular cyclic ambits of orders 1–4, clopen families of sizes 1–3 and color-vector lengths 1–4 satisfy the bound 2^|Σ|. Repetitions of clopens are included.
- External family: r = 2,…,6 and m = 2,…,6; translated membership is computed from the shift operation independently of the residue formula. For each m all translation residues are checked, with h from −3m to 3m, and the minimum and maximum color counts both equal m.
- Dense basepoint: concatenations contain every nonempty alphabet word of lengths 1–4 for r = 2,3,4.
- Orbit mechanism: 27 simultaneous zero-window checks on coordinates −2,…,2, covering arbitrary clopen Boolean truth tables supported in that window.
- Negative control: repeating a single periodic parameter yields two colors rather than m colors when m > 2.

These are bounded exact checks. They do not establish compactness, existence of a dense orbit, the arbitrary-clopen argument, or unbounded external degree. Those require the general mathematical argument and its later Lean proof. This report contains no Lean settlement proof.
