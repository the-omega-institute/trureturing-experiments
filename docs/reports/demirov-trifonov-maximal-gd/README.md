# Maximal GD extensions of arithmetic and geometric progressions

Source: Georgi Demirov and Ognian Trifonov, *Integer Sequences which Are Closed with Respect to Multiplication and whose Sumset does not Intersect the Sequence*, arXiv:2607.09937v1, Problem 7 restricted to Examples 2 and 3. The perfect-power family of Example 1 is outside scope.

A GD set is an infinite set of positive integers closed under multiplication with no three members satisfying a+b=c, allowing a=b. The proposed construction chooses the largest divisor d of p-1 for which the d-th roots of 1 modulo p are sum-free, and takes all positive integers in those residue classes. The polynomial gcd criterion is checked independently against the additive residue condition.

Run with Python 3, using only the standard library and exact integer arithmetic:

```sh
python3 docs/reports/demirov-trifonov-maximal-gd/check.py
```

The stored output is `check-output.txt`. A successful run exits 0 and ends with `ALL_CHECKS_OK`.

- Fidelity: source Theorem 5(i), the proper extension of V_k by the residue classes +1 and -1 for every 4<=k<=100; multiplication and sum-freeness are checked on the entire finite quotient. Examples 2 and 3 are also checked for the displayed parameter and exponent ranges.
- Construction: all 46 primes p<=199 and all 354 divisor orders. Polynomial Euclidean gcds, direct additive tests, subgroup maximality and the first six periods of the increasing enumeration agree. Every excluded nonzero residue generates a subgroup that fails sum-freeness.
- Classification: all 2,603 subsets of the nonzero residues containing 1 for primes p<=13 are tested for multiplicative closure and compared with the complete divisor-root subgroup list.
- Lifting: 14,884 exact certificates for the integer lifting identity, for primes p<=43, using both canonical and shifted positive representatives.
- Family inclusions: V_k for 2<=k<=150 and W_k for 3<=k<=150, with exponents or progression indices 0..100.

These finite checks support a preregistration. They do not prove the universal claim or maximality against arbitrary nonperiodic infinite extensions. That step requires the proposed integer lifting argument and a later Lean proof. No floating point arithmetic or numerical approximation is used as evidence.
