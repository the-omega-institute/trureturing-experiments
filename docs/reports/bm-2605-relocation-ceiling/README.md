# Exact checks for the preferential-relocation ceiling

Source: Jean Bertoin and Martin Minchev, *Preferential relocations enhance survival for Markov chains with killing*, arXiv:2605.29030v1, final question in section 5. The question asks which matrices satisfy the uniform spectral-radius ceiling. The proposed answer is every matrix in the standing class.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/bm-2605-relocation-ceiling/check.py
```

The program uses only the Python standard library and exact integers and fractions. For three irreducible aperiodic two-state substochastic matrices with different row sums (including the source's displayed matrix and a matrix with a zero entry), it constructs the finite-memory transition matrix directly. It checks every relocation probability law supported in {0,1,2} whose common reduced denominator is at most four: 22 laws per matrix, 66 cases. For each matrix M and positive rational vector v, min(Mv/v) and max(Mv/v) enclose its Perron root. Iterating rational matrix-vector multiplication tightens those certified intervals; all displayed endpoints are rounded outward using integer division.

An independent rational probability p gives q=p sigma and a=p/q, with p sigma_a=p. Thus the source's supremum is at least exp(F(p)), where F(p)=sum p_i log(q_i/p_i). For a common denominator D, exp(D F(p)) is rational. Every comparison of the certified upper spectral bound with this entropy value is made by raising to D, without evaluating a logarithm, exponential or floating-point number. This uses only the finite Perron correspondence, not the proposed history-ceiling theorem or its optimizer identity.

Fidelity checks reproduce the source matrix's benchmark root (65+sqrt(193))/100, its equality to the transfer matrix for delta_0, and the source's Theorem 4.1(ii) strict enhancement in all 57 non-Dirac finite-memory cases. Additional checks test the boundary perturbation at all five finite-entropy simplex vertices and 1,644 exact prefix telescoping identities and tail-error inequalities. A zero-diagonal vertex has entropy minus infinity and is explicitly excluded from the finite logarithmic comparison.

The adjacent output.txt is the program's complete deterministic output and ends in ALL_CHECKS_OK. These checks support source fidelity and the proposed argument on the stated finite grid. They do not establish the universal infinite-memory claim, interior attainment, or prior-art status. The proposed proof uses a general entropy maximizer and finite-cylinder supersolutions; enumeration is supplementary evidence.
