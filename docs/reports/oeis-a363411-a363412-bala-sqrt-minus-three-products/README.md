# OEIS A363411, A363412: four of Bala's conjectures on ∏ (1 + k√−3)

Write P_n = ∏_{k=0..n} (1 + k√−3) = A_n + B_n√−3, so A363411 = A_n and A363412 = B_n. The four conjectures checked are:
- **A363412.** 3 ∤ B_{3n+1}, and v₃(B_{3n}) = v₃(B_{3n−1}).
- **A363411.** v₅(A_{n+2}) = [5 | n], and the type 1 primes (dividing no A_n) up to 101 are 3, 11, 17, 23, 29, 41, 47, 53, 59, 83, 89, 101.

The proofs rest on four facts:
- **A divisibility lemma.** m(m+1)/2 divides every odd elementary symmetric sum e_j(1, …, m). Hence B_m = (m(m+1)/2)(1 + 3K), so v₃(B_m) = v₃(m(m+1)/2).
- **Inert primes.** For an inert prime p, ∏_{k mod p} (1 + k√−3) ≡ 2 (mod p), so P_{n+p} ≡ 2P_n.
- **Modulo 25.** ∏_{k mod 25} (1 + k√−3) ≡ 7 (mod 25).
- **The prime 3.** A_n ≡ 1 (mod 3).

`check.py` (Python 3 with sympy) verifies:
- the lemma for m ≤ L;
- the four statements for n ≤ N with exact integers;
- the residue-system products and the finite residue checks;
- the published data.

Run `python3 check.py 3000 300`. Exit 0 and `RESULT PASS` mean every check holds.
