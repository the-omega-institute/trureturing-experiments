# Pure-state distances for three cyclic sites

Fix three distinct sites $x,y,z$ of $\mathbb Z_d$, with positive coordinate amplitudes of moduli $a,b,c$ and $a^2+b^2+c^2=1$. Let $A=ab,\ B=bc,\ C=ca$. Write $R$ for the minimum Euclidean distance from the uniform Fourier probabilities $u=(1/d,\ldots,1/d)$ over pure states with these fixed moduli. The minimum rank in the last column is over all density matrices with the specified diagonal and Fourier probabilities exactly $u$.

| Configuration, after relabelling | $dR^2$ | Minimum rank at $u$ |
|---|---|---|
| All six nonzero ordered differences distinct | $2(A^2+B^2+C^2)$ | 3 |
| Only $x,z$ antipodal, no progression | $2(A^2+B^2)$ | 2 |
| $x,y,z$ in progression, step order at least five | $2(A-B)^2+2C^2$ | 2 |
| Progression, step order four | $2(A-B)^2$ | 1 if $a=c$, otherwise 2 |
| Coset of the order-three subgroup | $2(2\max(A,B,C)-A-B-C)_+^2$ | 1 if the three lengths close a triangle, otherwise 2 |

Indeed, two coincident directed differences on distinct underlying edges force a three-term progression with their shared vertex at the midpoint. On a progression with step $h$, the remaining difference is $2h$. The classes $\pm h$ and $\pm2h$ merge at order three, whereas $2h=-2h$ at order four. A single edge can collide with its own reverse only at the cyclic group's unique nonzero involution.

For pure amplitudes $v_x$, write $b_j=\sum_{s-t=j}v_s\overline{v_t}$. Parseval gives $d\|q-u\|_2^2=\sum_{j\ne0}|b_j|^2$. A singleton edge and its reverse contribute twice the squared edge length. An antipodal edge contributes four times its squared real part, which vanishes in imaginary relative phase. With two consecutive edges in the same bin the squared modulus of their sum is at least $(A-B)^2$. Real amplitudes $(a,b,-c)$ attain this in a progression of step order at least five. At order four, $(a,e^{-i\pi/4}b,ic)$ simultaneously opposes the two consecutive contributions and kills the real part of the antipodal coherence. At order three the smallest resultant of the three edge vectors has length $(2\max(A,B,C)-A-B-C)_+$. Rotating all edge phases together preserves the resultant's modulus and can impose zero sum of the three phases modulo $2\pi$, so vertex phases exist.

At a uniform Fourier readout with six distinct differences the density matrix is forced to be diagonal. Whenever a difference collides, there is a nonzero Hermitian perturbation with all diagonal entries and all cyclic sums zero. An antipodal edge permits a purely imaginary perturbation. Two distinct edges at the same nonself-inverse difference permit opposite entries together with their conjugates on the reverse edges. Start at the positive diagonal on the three sites and move along this trace-zero direction to the first singular positive semidefinite matrix. Thus every collision permits rank at most two. The pure-state minimum specifies when rank one suffices.

## Flatness after leakage onto one additional site

Let the original two sites retain conditional probabilities $\alpha,1-\alpha$ with $0<\alpha<1$, and the third carry mass $0<t<1$. In an order-three coset put $s=\sqrt{\alpha(1-\alpha)}$. Factoring the pair amplitudes in the triangle inequalities gives
$
|\sqrt{\alpha}-\sqrt{1-\alpha}|\,\sqrt{t}\le s\sqrt{1-t}
\le(\sqrt{\alpha}+\sqrt{1-\alpha})\sqrt t .
$
Thus pure flatness holds exactly on the interval
$
\left(\frac{s}{1+s}\right)^2\le t\le\left(\frac{s}{1-s}\right)^2 .
$
For equal conditional weights the lower endpoint is $1/9$; for $\alpha=4/5$ both endpoints $4/49$ and $4/9$ are attained. In an order-four progression the two endpoint masses must be equal. If the third site is the midpoint, flatness therefore requires equal original conditional weights and then holds at every $0<t<1$. If the third site is an endpoint and the other endpoint's original conditional mass is $\alpha$, the unique leakage is $t=\alpha/(1+\alpha)$. All other configurations have positive radius for strictly positive masses.

The replay enumerates all 66,045 unordered triples through dimension 36 and checks the difference multiplicities independently from the midpoint test. It verifies the squared radii of explicit complex and real states directly from their cyclic correlations, including embeddings in larger dimensions and both endpoints for unequal leakage. Run with Python and SymPy. [Exact checks](ch08-three-site-pure-radius.py)
