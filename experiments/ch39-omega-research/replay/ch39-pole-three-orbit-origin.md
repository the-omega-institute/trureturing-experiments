# Pole-assisted recovery after one loss in each of three orbits

The nineteen-station shell has six ordinary orbits of size three and a pole. Fix a known common clock background and a common phase reference for the five charges $1,-4,2,-3,6$. Use one unit vector for the lowest line of each cyclic character and charge. With the pole retained and one station lost in each of three different ordinary orbits, the sixteen readings recover all fifteen complex amplitudes whenever the two clock angles have absolute value at most $10^{-6}$ radians and $0\leq\delta\leq10^{-8}$. Over this entire closed box the normalized smallest singular value exceeds $1.32\times10^{-5}$. Recovery on the full common-clock torus for these losses remains open. [Book research question](https://github.com/the-omega-institute/trureturing-books/issues/1856)

## Coverage and error bound

A common cyclic rotation acts by unitary permutations of the stations and by phases on the input lines, so it preserves singular values. Choose three of six orbits, place the loss in the first orbit at a fixed station, and allow each remaining loss three relative positions. The $20\cdot9=180$ representatives cover all 540 masks.

For each charge and cyclic block, write the angle perturbation of the zero-clock Hamiltonian as
$$
\sum_{j=1}^2\left[
\bigl(e^{iq\theta_j}-1\bigr)L_j
+\bigl(e^{-iq\theta_j}-1\bigr)L_j^*
\right].
$$
The three softened unit spokes form a star with nonzero Laplacian eigenvalues $1,1,4$, so the softening changes the Hamiltonian by at most $4\delta$ in operator norm. The replay uses interval Hermitian elimination to count one eigenvalue below $9.6$ plus a uniform perturbation bound at the origin. Weyl's inequality then puts every eigenvalue above the lowest line strictly above $9.6$ throughout the box.

The floating-point eigensolver supplies a trial vector $v$ and an energy $e<8.685$. Its residual at the origin is checked in Arb. The change in the residual on the box is evaluated on $v$. Convexity bounds the first-order clock term by the largest norm at a corner of the angle square. The two pairs of opposite corners have the same norm. The exponential remainder obeys $\lvert e^{is}-1-is\rvert\le s^2/2$ for real $s$.

Let $P$ project onto the lowest line at a background in the box. Dividing the residual upper bound by $9.6-8.685=0.915$ bounds $\|v-Pv\|$. The distance to the unit vector $Pv/\|Pv\|$ is at most
$$
\|v-Pv\| + \bigl|\|v\|-1\bigr|
+ \frac{\|v-Pv\|^2}{\|v\|}.
$$
Use a Frobenius bound for the five column errors in each character. Before dropping any station rows, both the input and output character spaces are orthogonal, so the largest of these three bounds controls the operator norm of the full difference.

For each mask the replay checks an approximate left inverse $L$ of the retained trial matrix $M$. The bound $(1-\|I-LM\|_F)/\|L\|_F$ exceeds $3.0030\times10^{-5}$ for every trial matrix. The transfer bound is below $1.6732\times10^{-5}$. Arb checks the strict $1.32\times10^{-5}$ margin on their difference. The choice of the fifteen column phases does not affect the result.

## Reproduction

The input model and its dependencies are fixed by the experiment commit below. Run the replay with assertions enabled. Python 3.12, python-flint 0.9.0, and NumPy 2.5.5 reproduced the bound.

```sh
git worktree add --detach ../omega-origin-input cb58c8e64ecea67184e0c1d68a690be098a7a384
python3 experiments/ch39-omega-research/replay/ch39-pole-three-orbit-origin.py \
  ../omega-origin-input/experiments/ch39-omega-research/replay
```
