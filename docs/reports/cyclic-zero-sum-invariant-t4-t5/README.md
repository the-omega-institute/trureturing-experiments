# Cyclic zero-sum invariant: T4 and T5

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The output is in `output.txt`; success ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce the first instances of Theorem 5.1 in Wenkai Yang, [arXiv:2610.13103v1](https://arxiv.org/abs/2610.13103v1): its T2/T3 constructions, exact ceiling and B-value, and absence of a zero-sum subsequence exhausting a support element.

For the candidate parametrizations T4: d = 30030t + 2192, t >= 0; T5: d = 30030t - 2192, t >= 1; n = d² - 64, the program checks the two integer identities symbolically. It verifies the supplied first instances, coprimality, the ceiling formula, weak regularity, and both forbidden boundary residues for every integer parameter through 1000, using exact integers throughout. No enumeration of the enormous sequences is required: with support {e, xe}, a subsequence exhausting a support element lies on one of two faces, each decided by a unique residue modulo n.

These checks support the proposed lower bound N(C_n) >= 2n - 2d + 1. Equality also uses the published universal cyclic upper bound in Gao et al., *Representation of zero-sum invariants by sets of zero-sum sequences over a finite abelian group II*, Journal of Number Theory 241 (2022), 738–760, Theorem 1.1(1), DOI [10.1016/j.jnt.2022.03.016](https://doi.org/10.1016/j.jnt.2022.03.016), recalled by Jiang–Yang (2024) and Yang (2026). This program does not prove that upper bound or compute N by searching all weak-regular sequences.

The parametrizations agree with Jiang–Yang (2024), Lemma 5.2(4)–(5), Journal of Number Theory 260, 198–211, DOI [10.1016/j.jnt.2024.02.005](https://doi.org/10.1016/j.jnt.2024.02.005). The full journal text is available from [Tianjin University](https://cam.tju.edu.cn/en/research/downAchiev.php?id=798). Their literal formulas are T4 = {30030(30030t² + 4384t + 160)} for t >= 0 and T5 = {30030(30030t² - 4384t + 160)} for t >= 1; the symbolic checks also establish their equivalence to n = d² - 64.
