# Sutherland Remark 6.6

Run `python3 check.py` with Python 3.12 or newer. The check uses integer quaternion multiplication and bit vectors for Boolean factors. The saved output is in `check.out`.

The fidelity checks reproduce the codegree identity in §1.1 and count all permutations fixing the identity and preserving the paper's complete inverse-pair coloring on each tested group. The count of eight reproduces the quaternion–Boolean case of Theorem 1.1 on finite truncations. The obstruction check merges all forbidden difference colors and checks that the identity–central-involution transposition preserves every remaining color. Since each admissible connection set is a union of these colors, the check covers every such set on the listed finite groups. It also checks a generating connection set and an inverse-pair color that the transposition changes. The infinite argument is given in the preregistration draft.

Source [arXiv:2610.01049v4, §1.1 and Remark 6.6](https://arxiv.org/abs/2610.01049v4).
