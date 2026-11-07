# Contiguous entropy maximization: exact finite checks

Source: Roberto Bruno and Ugo Vaccaro, *Constrained Maximum Entropy Contiguous Aggregations*, arXiv:2608.25533v1, Sections I, IV and V, DOI 10.1109/ISIT62367.2026.11654039. The source asks whether maximizing Shannon entropy of a contiguous aggregation under an entropy ceiling is NP-hard. Output length is unrestricted.

Run with Python 3 and its standard library:

```sh
python3 check.py
```

`check-output.txt` is the resulting output. Every arithmetic operation uses integers or `fractions.Fraction`. The source's merge-loss identity is checked in 660 positive rational cases by exact coefficients of logarithms of primes. No numerical approximation is needed for these identities.

For positive SUBSET SUM data, let n be the number of weights, S their sum, N=2n, M=S+1, w=(M+a₁,…,M+aₙ,M,…,M), W=Σwᵢ, B=nM+T, K=8N, c=1/(M·2ᴷ), xᵢ=cwᵢ/2 and b=(1−cW)/(N−1). The ordered vector has N twin pairs, with one b between successive pairs: (x₁,x₁,b,x₂,x₂,b,…,b,x_N,x_N). The program checks forced cardinality and the minimum weight sum at least B for all aᵢ∈{1,2,3,4}, 1≤n≤3, and every 1≤T≤S: 570 instances.

For the five listed examples, the program enumerates every contiguous partition. Logarithms are enclosed by rational partial sums of the atanh series with a geometric tail bound, then rounded outward to a fixed binary denominator. Every feasibility comparison and every strict comparison that excludes separator crossings is decided by these exact enclosures. Unresolved comparisons cause failure. Rational ceilings use an entropy enclosure midpoint h and R=h−cB+c/8; the program checks the required error bound. It also checks the rational entropy-window decision variant with lower endpoint h−cB−c/8.

| Original data | Target | Gadget B | Optimal merged weight U | Partitions |
|---|---:|---:|---:|---:|
| (1) | 1 | 3 | 3 | 16 |
| (1,3) | 2 | 12 | 13 | 1024 |
| (1,3) | 3 | 13 | 13 | 1024 |
| (2,3,7) | 5 | 44 | 44 | 65536 |
| (2,4,8) | 5 | 50 | 51 | 65536 |

All maximizing partitions in these examples respect the separators. These are finite checks of the proposed general reduction. They establish neither a uniform NP-hardness theorem nor any machine-level polynomial-time bound, and contain no Lean proof.
