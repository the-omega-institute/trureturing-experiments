# Krebs–Leyva Conjecture 1.2

Run `python3 check.py` with Python 3.12 or newer. It uses the standard library and exact integer arithmetic. The output is in [check.out](check.out).

The source check reproduces the one-row gcd lemma in §2.1 by computing the chromatic numbers of the associated finite Cayley graphs for 624 two-column matrices. This includes the source’s convention that a loop precludes a proper coloring. [arXiv:2511.03028v1](https://arxiv.org/abs/2511.03028v1)

The family check uses the identity minor to eliminate pivot coordinates over the integers. For each sampled pair it verifies the quotient generators, a clique of size \(r+1\) and a compatible coloring modulo \(r+1\). The sampled ranks are \(3\) through \(12\) and the sampled values of \(k\) are \(1\) through \(8\).
