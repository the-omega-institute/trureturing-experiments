# Prescribed poles and bounded inflection components

Run `python3 check.py` with Python 3 and SymPy. The program uses integer, rational and symbolic complex arithmetic exclusively. Its output is in `check.out`.

Source: Boris Shapiro and Guillaume Tahar, *Inflection curves of rational vector fields*, arXiv:2605.28266v1, Sections 3, 4 and 11, and Problem 12.4. Source URL: https://arxiv.org/html/2605.28266v1. The check reproduces the quotient derivative formula, the coordinate-axis curves for R=1/z and R=z+1/z, the polynomial examples R=z^d, and the leading derivative term at a pole.

The proposed answer is n bounded components when infinity is a pole of R, and n−1 when infinity is regular, where n is the number of distinct finite poles. The general proof requires topology of the spherical inverse image of the real projective line. The program supplies exact construction certificates; it does not compute arbitrary real algebraic connected components or verify the general topological upper bound.

The finite construction domain has n from 0 through 3, finite pole orders from 1 through 3, infinity order from 0 through 3, the zero divisor omitted, and finite pole positions 0 through n−1. For a polynomial background use R₀=i(z+1)^m. For a proper background retain R₀=i/z^s and insert the other principal parts. All disks have radius r=1/32. A single positive rational epsilon makes the perturbation smaller than the signed background margin on every boundary. The program verifies the resulting denominator, absence of cancellation at every prescribed finite pole, and the exact infinity pole order.

For the polynomial background, the margin follows from the binomial bound on (a+u)^k−a^k for |u|≤r. For the proper background, it follows from

    |(a+u)^(-k)−a^(-k)| ≤ ((a+r)^k−a^k)/(a^k(a−r)^k).

Here a>r>0, so all denominators are nonzero. Each inserted derivative term is bounded by its pole order divided by the appropriate lower distance to that pole raised to one more than its order. Strict margins certify that no inflection point occurs on a boundary. Once the pole-inclusion and connected-component isolation lemmas are proved, these certificates produce the required number of distinct bounded components. The exact component maximum additionally uses the general upper bound; the finite checks alone do not prove it.

Two further explicit checks use radius 1/4: i*z+(1/100)(1/z+1/(z−1)+1/(z−2)) has boundary margin at least 181/225 on its three disks; i/z+1/(100(z−1)) has a negative derivative-imaginary margin at most −28/625 on the disk about 1. Those are uniform inequalities, not point samples.
