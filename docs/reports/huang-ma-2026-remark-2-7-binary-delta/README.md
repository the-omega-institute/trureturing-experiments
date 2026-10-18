# Huang–Ma (arXiv:2602.03194) Remark 2.7: δ is binary for every fixed n and skew-symmetrizer

Source: M. Huang and Q. Ma, "A mutation invariant for skew-symmetrizable matrices", arXiv:2602.03194v1 (J. Algebra 712 (2027) 546–558). Definition 1.3: δ(B) = det S(B) mod 4. Remark 2.7: "We also conjecture that for any fixed n and any choice of skew-symmetrizer d₁, …, dₙ, the invariant δ is binary."

Theorem checked: with Π dᵢ = 2^e·u (u odd), δ(B) lies in {0, 2} if n is odd or e is odd, and in {0, (−1)^{n/2}·u mod 4} if n is even and e is even.

## Program

`check.py` (standard library only). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.

1. **The companion matrix.** Put Aᵢᵢ = 2, Aᵢⱼ = bᵢⱼ for i < j and Aᵢⱼ = −bᵢⱼ for i > j. Then det S(B) = det A for 400 random matrices (floating-point determinant of S(B), relative difference below 1e−9).
2. **Lemma 1.** det A ≡ (−1)^m·det B (mod 4) for n = 2m holds in 900 random cases.
3. **Residue sets, exhaustive for n ≤ 4.** Every pair parameter t ∈ [−2, 2], with dᵢbᵢⱼ = −dⱼbⱼᵢ = t·lcm(dᵢ, dⱼ). This covers 26 symmetrizers, up to 15625 matrices each.
4. **Residue sets, random for n = 5, 6.** 3000 random matrices for each of 12 symmetrizers.

In checks 3 and 4, every residue set has at most two elements and lies in the predicted set.

Reading: `run.log` (exit 0, `RESULT PASS`, 1.0 s).
