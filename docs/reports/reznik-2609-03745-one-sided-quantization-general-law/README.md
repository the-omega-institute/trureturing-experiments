# One-sided quantization: exact checks

Run `python3 check.py`. It uses only Python's standard library and rational arithmetic and must end with `ALL_CHECKS_OK`.

The fidelity check reproduces Example 1 of Yuriy A. Reznik, *A Quantization Problem Posed by Adaptive Streaming*, arXiv:2609.03745v1: for uniform bandwidth and linear quality on [0,1], the pinned zero rung gives n+1 cells and optimal gap 1/(2(n+1)). The check evaluates the nonnegative variance identity on 8,568 rational partitions for n=0,...,5.

For finite atomic laws, dynamic programming and exhaustive codebook search agree on 1,815 law/budget cases: up to five equally spaced support points, all nonnegative mass vectors with denominator six, and every budget up to support size. Atoms at a rung use that rung. A ladder containing all positive support points has zero gap.

The mixed-law check uses Q(B)=B/(1+B) and bandwidth law

    mu(dB) = 1_[0,1](B)/(1+B)^2 dB + (1/2) delta_1(dB).

Its continuous and atomic masses are both 1/2. In quality coordinates it is (1/2) U[0,1/2] + (1/2) delta_(1/2). Write b=1/2, let h_0,...,h_n be consecutive cell widths, and let t=h_n. The loss is sum(h_i^2)/(4b)+t/2. For n>=1,

    loss - b/(4n)
      = sum_{i<n}(h_i-(b-t)/n)^2/(4b)
        + t(n-1)/(2n) + t^2(n+1)/(4bn).

Every term on the right is nonnegative. Equality is attained at quality rungs i/(2n), hence bandwidth rungs i/(2n-i). Thus G_n*=1/(8n). The candidate constant is (1/2)(integral_0^1 (1+B)^(-2) dB)^2=1/8. The program checks the identity on 12,869 rational partitions for n=1,...,8 and checks every equality ladder exactly.

These finite checks and the displayed elementary mixed-law calculation do not prove the universal general-law asymptotic claim and do not certify literature novelty. No Lean proof is included.
