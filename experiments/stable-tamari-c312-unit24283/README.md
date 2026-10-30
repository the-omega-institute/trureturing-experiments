# Stable Tamari Conjecture 3.12: finite experiments

Conjecture 3.12 remains unresolved. These are finite computations and an exact
obstruction to an auxiliary two-summand law, not a full solution or a KPI delivery.

## Source and literal target

H. Chau, S. Daugherty, S. Grate, M. Haiman, J. C. Martínez Mori, and A. Pun,
*The stable Tamari lattice*, [arXiv:2610.39257v1](https://arxiv.org/abs/2610.39257v1),
Conjecture 3.12, source license [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The version is fixed; no article HTML, transcripts, or source archive is included.
See [SOURCE-LICENSE.md](SOURCE-LICENSE.md) for attribution and adaptation details.

The source question is: “For every $\mathbf b\in\mathbb N^n$,
$f_{\mathbf b}(t)=f_{\operatorname{rev}(\mathbf b)}(t)$.”
Here $\mathbb N=\{0,1,2,\ldots\}$ and $n$ ranges over every positive integer.
For $\mathbf a\in\mathbb N^n$, the position relation $i\leq_{\mathbf a}j$
means $i=j$, or $i<j$ and $a_k>a_j$ for every $i\leq k<j$.
The stable Tamari order $\mathbf a\leq_T\mathbf b$ means both
$a_i\leq b_i$ for all $i$ and
$a_i-a_j\leq b_i-b_j$ for every $i\leq_{\mathbf a}j$.
The interval height $\mathsf{ht}(\mathbf a,\mathbf b)$ is the maximum number
of strict steps in any chain from $\mathbf a$ to $\mathbf b$, equivalently
the maximum number of chain elements minus one. The polynomial is

$$f_{\mathbf b}(t)=\sum_{\mathbf a\leq_T\mathbf b}
t^{\mathsf{ht}(\mathbf a,\mathbf b)},\qquad
\operatorname{rev}(b_1,\ldots,b_n)=(b_n,\ldots,b_1).$$

The target requires equality at every degree $d\geq0$, including the constant
coefficient and the longest-chain degrees, for all nonnegative vectors.

## Retained finite computations

`target_search.py` and `target_scan.cpp` enumerate the entire coordinate box
$\prod_i\{0,\ldots,b_i\}$ and filter by the literal stable order.
Height uses the source Theorem 3.7 forest formula: set
$p(i)=\min\{j>i:a_j<a_i\}$ when this set is nonempty; sum
$(b_i-a_i)-(b_{p(i)}-a_{p(i)})$ when $p(i)$ exists, and $b_i-a_i$ otherwise.
This formula's source fidelity was externally reviewed. That review is not an
independent replay of the large searches and is not Lean/kernel verification.
The data report direct longest-chain dynamic-program controls for all upper
vectors in boxes $(n,M)=(1,8),(2,5),(3,3),(4,2)$, totaling 9, 36, 64, and 81
vectors respectively; the control implementation is not among these retained inputs.

`target_probe_data.json` stores full coefficient lists in increasing degree order,
including degree zero and all interior zeros. Coefficients above the final entry
are zero. It also stores the following completed exhaustive ranges:

| n | maximum coordinate M | upper vectors checked | unequal reversal polynomials |
|---|---|---|---|
| 1 | 12 | 13 | 0 |
| 2 | 10 | 121 | 0 |
| 3 | 25 | 17,576 | 0 |
| 4 | 12 | 28,561 | 0 |
| 5 | 8 | 59,049 | 0 |
| 6 | 5 | 46,656 | 0 |
| 7 | 3 | 16,384 | 0 |
| 8 | 2 | 6,561 | 0 |

`target_random.cpp` uses `mt19937` seed 20261030, proposing dimensions uniformly
on 3..8 and each coordinate uniformly on 0..20. It rejects proposals whose
coordinate-box product exceeds 100000, and reports 20000 accepted samples with
no discrepancy. Accepted dimensions/vectors are **not uniform**; the number of
unique accepted vectors was not measured. Exact reproducibility of distribution
mapping can depend on the C++ standard-library implementation. C++ counters are
`uint64_t`; the recorded boxes bound the number of contributing vectors below
its capacity. Arbitrary enlarged command-line ranges have no overflow guarantee.

These observations supply neither a universal proof nor a literal finite
coefficient-fiber counterexample to Conjecture 3.12.

## Auxiliary exchange certificate

`context_exchange.py` constructs Theorem 4.17 rational summands from laminar
intervals $[h_r,r]$. In its encoding,
$a_i=\#\{r:h_r\leq i\leq r\}-1$, with $a_0=0$.
For each $j$ with $h_j=j$, and $a_{j-1}\leq\gamma\leq a_j$, the uncolored
factor is $1-x_j\cdots x_k$, where $k$ is maximal such that every
$a_j,\ldots,a_k\geq\gamma$. Colored factors are $1-tx_{h_r}\cdots x_r$;
the numerator is $\prod_i x_i^{a_i}$. All exponent vectors include the final
$t$ coordinate. The stored data contain complete weights and cleared polynomials.

The proposed auxiliary law was invariance of $W_{h^0}+W_{h^1}$ under the
adjacent variable swap $\sigma_p$. Its bounded test was $2\leq n\leq4$,
in increasing $n$, lexicographic $h^0$, ascending $p$, stopping at the first
nonzero cleared difference. Three eligible cases were inspected; $n=4$ was
not inspected. At $n=3$, $h^0=(1,2,3)$, $h^1=(1,1,3)$, $p=1$, put
$(x_1,x_2,x_3)=(x,y,z)$. The complete two weights are

$$W_{h^0}=\frac{1}{(1-xyz)(1-yz)(1-z)(1-tx)(1-ty)(1-tz)},$$
$$W_{h^1}=\frac{x}{(1-xyz)(1-x)(1-z)(1-tx)(1-txy)(1-tz)}.$$

For $D=W_{h^0}(x,y,z;t)+W_{h^1}(x,y,z;t)
-W_{h^0}(y,x,z;t)-W_{h^1}(y,x,z;t)$, the exact certificate is

$$QD=(x-y)(1-xyz)(1-z)(1-txy),$$
$$Q=(1-xyz)(1-yz)(1-z)(1-tx)(1-ty)(1-tz)
(1-x)(1-txy)(1-xz)(1-y).$$

Thus $D=(x-y)/((1-tx)(1-ty)(1-tz)(1-x)(1-y)(1-xz)(1-yz))\ne0$.
At $(x,y,z,t)=(1/5,1/7,1/11,2)$, the unswapped sum is
$2235233/590976$, the swapped sum is $990143/279936$, and their difference
is $326095/1329696$. This refutes only that auxiliary law.

`independent_review.py` is a separate fixed-witness implementation using direct
subset expansion, exact integer polynomial comparison, and rational arithmetic.
Its retained JSON supplies the full replay certificate. Its author had the supplied
outcome and permitted program/data available; this was not a blind test or kernel
verification. The replay does not enumerate upper vectors or search new hypotheses.
`reviewer_formula_diagnostic.py` and its JSON retain a separate exact $n=2$
diagnostic: the full two-term sum is swap invariant, while individual summands
need not permute; naive interval reflection can create an endpoint collision.

## Portability and limits

Python programs use only the standard library; C++ programs require C++17.
The C++ scanner differs from its supplied input only by the explicit missing
`<functional>` include; no computed result was changed by that repair.
Other portability changes remove machine paths and metadata, use a relative
certificate-data path, and omit timestamps and source-cache dependencies.
The scanner and random programs remain search programs; retaining them does
not authorize or imply a new run. The original finite searches were not rerun
for this archive. Only compilation, a one-vector scanner smoke replay, and the
fixed certificate/diagnostic replays were used for archival validation.

The remaining bridge is a valid mechanism for **all** whole-polynomial reversal
coefficient fibers; the refuted two-summand exchange cannot supply it. No new
successor hypothesis is selected here. No Lean closure, independent large-search
replay, global novelty determination, or full-resolution claim accompanies these files.
