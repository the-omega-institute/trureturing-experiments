# A chiral bosonic phase with stable and unstable bulk regions

This exact check supports a preregistration draft for the first open question in Section V of Mariam Ughrelidze, Lorenza Viola and Emilio Cobanera, *Decoupling band topology from criticality in bosonic systems*, arXiv:2607.01334v2, J. Appl. Phys. 140, 124401, DOI 10.1063/5.0350468. It is an independent mathematical construction, not a Lean proof or a claim of literature completeness. The underlying nonzero-frequency instability mechanism is established prior art, including Peano–Schulz-Baldes, arXiv:1708.08798; the proposed settlement concerns the continuous path in the specified gapped chiral phase.

Run with Python >=3.12 and SymPy 1.14.0:

```sh
python3 docs/reports/uvc2607-bosonic-stable-unstable-phase/check.py
```

`check.out` is the complete output. Every comparison uses exact symbolic algebra, integers, rational numbers or radicals. The first check reproduces the source interpolation model's determinant and its doubly degenerate bands. The remaining checks establish the displayed polynomial identities, symmetry relations and eigenbasis certificates. They do not numerically sample momentum. The inequalities and winding argument below apply to all momenta and all parameters in their stated domain.

Choose real a,c>0, a≠c, and 1/2<s<1. Set A=diag(a,-c), D=[[0,t],[t,0]], B=[[A,D],[-D,-A]], J=diag(I₂,-I₂), C=[[0,I₂],[I₂,0]], f=1-s+s exp(ik), b=fB, and g=[[0,conj(f)B],[fB,0]]. Since JBᵀJ=B and CBC=-B, these are exactly the source's bosonic chiral constraints. The physical Hessian diag(J,J)g is Hermitian. Both g and its Hessian have only Fourier modes -1,0,1 and describe a finite-range translation-invariant quadratic Hamiltonian with four physical orbitals per unit cell. The Nambu carrier and R=diag(-I₄,I₄) remain fixed.

Writing α=(a+c)/2 and β=(a-c)/2, the coordinate order (1,4,2,3) splits B into [[a,t],[-t,c]] and [[-c,t],[-t,-a]]. Their eigenvalues are ±α±sqrt(β²-t²), and det B=(ac+t²)²>0. For 0<t<|β|, the eigenvalues are distinct and real: α>|β| follows from α²-β²=ac>0. For t>|β|, the square root is imaginary and nonzero. At t=|β|, the blocks are defective with eigenvalues ±α, so this is an exceptional point at nonzero frequency under the source's terminology. It should not be called a diagonalizable Krein collision.

The exact tensor factorization g=[[0,conj(f)],[f,0]]⊗B transfers these stability conclusions to g for every momentum, even though g's bands are repeated. Its first factor has an explicit eigenbasis at eigenvalues ±|f|; tensor products of the two eigenbases give the required full eigenbasis in the stable region. Because

```text
|f|² = (2s-1)² + 2s(1-s)(1+cos k) >= (2s-1)² > 0,
det g = |f|⁸(ac+t²)⁴ > 0,
```

the entire interval [|β|/2,3|β|/2] is gapped at zero, including the defective collision. The open intervals (|β|/2,|β|) and (|β|,3|β|/2) are stable and unstable, respectively. The program checks symbolic eigenbases and the concrete values a=3,c=1,t=1/2 and t=3/2; at t=1 it verifies rank(B-2I)=rank(B+2I)=3 and det B=16.

The source index is the winding of det b, not det g. Here det b=f⁴(ac+t²)². Let r=(1-s)/s∈(0,1). Then f=s exp(ik)(1+r exp(-ik)), and the last factor has real part at least 1-r>0. Its principal argument ψ(k) is therefore continuous and 2π-periodic. The continuous argument lift θ(k)=4k+4ψ(k) of det b increases by 8π on [0,2π]. Thus the index is +4 for every t, with no appeal to a discrete cyclic-winding surrogate. The only division is by s>1/2. The discriminant square roots are used in their appropriate sign domains; at the collision their value is zero.

For an optional boundary calculation, the left half-line backward shift T gives F=(1-s)I+sT with kernel spanned by ρⁿ, ρ=-(1-s)/s and |ρ|<1. Its adjoint is injective, and invertibility of B makes the corresponding chiral half-line dynamical matrix have exactly four zero modes. The program checks the kernel recurrence only; square summability and completeness of this kernel follow from the recurrence and |ρ|<1. This boundary calculation is outside the preregistered core theorem and is not needed to answer the existence question.
