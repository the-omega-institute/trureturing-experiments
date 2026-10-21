# Exact biased boundary counting identities

`verify.py` checks finite identities used by the source-specific global counting proof in `FIB_ATOM_BIASED_BOUNDARY_SPECTRAL_COUNTING.md` of the main repository. It uses Python 3, its standard library, and NumPy. Run from any directory:

```sh
python3 verify.py
```

The program writes `results.json` beside itself. Its checks cover:

- First mass crossings at thresholds `a^k`, `k=1,...,16`, for `a=phi^-1`. Stopping includes exact equality. Coefficient pairs in `Q(a)` use `a^2=1-a`, with rational arithmetic and exact comparisons; no floating point classification is used.
- Every entry of the residual Gram identity for all 62 nonempty edge words of length at most five, at thresholds `a`, `a^2`, `a^3`, `a^5`: 15,376 matrix entries, including 3,924 nonzero entries. The left side is computed from conditional cylinder intersections; the right side is computed from strict-descendant biased suffix covariances.
- Rejection of a fair-law substitution: the biased child variance is `-1+2a`, whereas the fair variance is `1/4`.
- Original occurrence trees `T_n`, `n=0,...,12`, including the two singleton sources, with all 254 candidate nonempty addresses of length at most seven. Direct descendant enumeration agrees with the reached-address guard and descendant source count.
- The normalized finite covariance and inverse correspondence for `n=0,...,8`. With literal unit-edge difference matrix `D`, root-path sum `J`, `L=D^T D`, `N=|V|`, and `S=N I-11^T`, exact arbitrary-precision integer matrices check `D J=I`, `(S J)(S J)^T L=N S`, and `(S J)^T(S J)=N[N J^T J-c c^T]`, where `c=J^T 1`. These are respectively path inversion, `C_e C_e^T(N L)=Pi`, and `C_e^T C_e=Sigma_n` for the normalized Euclidean representative `C_e=Pi J/sqrt(N)`. They also verify the constant zero mode.

These results test finite identities and normalization. They do not certify a limit, an infinite counting estimate, an exact Weyl coefficient, a pointwise heat kernel, a physical clock or a physical dimension. The infinite result rests on the ordinary proof in the main theory source.

The program is original repository code under the root Apache-2.0 license. No third-party source or paper text is bundled.
