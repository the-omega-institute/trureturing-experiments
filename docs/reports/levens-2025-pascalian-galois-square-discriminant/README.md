# Levens (arXiv:2511.03082), Conjecture 6.4: the odd-index Galois group is not hyperoctahedral at m = 24

Source: "For n = 2m + 1, the Galois group of P_n(z) is the mth hyperoctahedral group". Here P_n(z) = Σ_k C(n, ⌊(n−k)/2⌋) z^{n−k}, and for odd n, P_n(z) = (1+z) q(z²) with q(x) = Σ_{j≤m} C(2m+1, j) x^j.

**The discriminant.**
- **Closed form.** Disc(q) = (−1)^{m(m−1)/2} (2m+1)(m+1)^{m−1} c^{m−2}, where c = C(2m+1, m). This is Filaseta–Moy (Colloq. Math. 2018), Lemma 3, at n = 2m+1 and k = m. It follows from Scherbak's identity (1+x)q′ − (2m+1)q = −(m+1)c x^m.
- **The square case.** At m = 24 the closed form is the square (7·5²³·c¹¹)².

**The contradiction.**
1. Every Galois automorphism fixes √Disc, so it induces even permutations of the 24 roots of q.
2. The splitting field of P_49 adjoins at most 2²⁴ further square roots.
3. Hence |Gal(P_49)| ≤ 2²⁴·24!/2 < |B_24| = 2²⁴·24!.

**The family.** The same holds whenever 2m+1 and m+1 are both squares (Pell: m = 24, 840, 28560, 970224, …).

## Program

`check.py` needs sympy (`python3 -m pip install sympy`). Run `python3 check.py`; it prints `RESULT PASS` when every check holds. It checks:
- the factorization P_49 = (1+z) q(z²);
- the differential identity for m ≤ 30;
- the discriminant closed form for 2 ≤ m ≤ 30;
- that Disc(q_24) is the square above;
- that among 2 ≤ m ≤ 40 only m = 24 gives a square.

It also lists the Pell family.

Reading: `run.log` (exit 0, `RESULT PASS`).
