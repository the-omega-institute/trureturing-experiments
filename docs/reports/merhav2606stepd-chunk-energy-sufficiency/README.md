# Chunk-energy sufficiency: exact finite checks

Source: Neri Merhav, *Grouped Reverse Importance Sampling for the Partition Function*, arXiv:2606.26748v1, Section 6; Proposition 2, equation (35), and Appendix C. The proposed settlement concerns the general measurable-space covariance criterion. These finite checks provide evidence for the proposed argument, not a proof of that unrestricted statement.

Run with Python 3, using only the standard library:

```sh
python3 check.py
```

The program uses integers and `fractions.Fraction` throughout. The integer parameter `r` represents `exp(beta)`, with integer energies, so Boltzmann probabilities are rational. No numerical logarithm or floating-point value enters any check.

The fidelity check reproduces Proposition 2's preservation of the normalizing constant and reduction of the non-overlapping second-moment ratio, together with the identities relating the ratio, variance and non-overlapping MSE constant in equations (8)-(10).

For each positive `k <= 4` on the binary alphabet with energies `(0,1)` and `exp(beta)=2`, and for each positive `k <= 3` on a three-state alphabet with energies `(0,0,1)` and `exp(beta)=3`, every divisor `d` is tested. There are 13 parameter cases. The local kernel is represented by a free coefficient vector; no finite selection of weights is substituted in the matrix check. Exact rational matrices verify that:

- one fibre-averaging projection preserves both the mean and the counting-measure normalizing constant for every local kernel;
- conditioning each cyclic window on all chunk energies applies the same local projection, including when distinct states have the same energy;
- covariance-pair counting equals the independently enumerated variance of a cyclic average on `2*k/d` chunks;
- the covariance sum before projection minus the covariance sum after projection equals a positive weighted Gram matrix of the residual window average.

The last identity tests covariance-sum contraction for all real kernels at those finite parameter cases. It asserts no separate monotonicity of each overlap covariance.

The nondegenerate binary example uses `k=4`, `d=2`, and `m(x)=1+x_1*x_3`. It reproduces `M=20`, `Z=3/2`, `mu=320/81`, the two covariance pairs, the strict MSE-constant improvement `1179/51200`, and the two cyclic-average variances at eight scalar samples. Full exact readings are in `output.txt`; successful execution ends with `ALL_CHECKS_OK`.
