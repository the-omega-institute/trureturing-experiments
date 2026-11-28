# Pole readout after three losses in two cyclic orbits

Retain the pole and omit two ordinary stations from one cyclic orbit and one ordinary station from another. For five signed charges and three cyclic characters, use unit vectors on the fifteen lowest lines as the input basis. The common two-clock background, station basis and cross-charge phase reference are known. For both clock angles on the entire torus and every spoke deficit in \([0,10^{-8}]\), the sixteen-row readout has smallest singular value strictly greater than \(5\times10^{-7}\). The 90 cyclic representatives cover all 270 layouts.

## Reflection

Let \(P\) exchange the first two coordinates of every site. It fixes the pole and conjugates the cyclic permutation to its inverse. The real dot-product term in the edge phase is preserved, and the oriented determinant changes sign. Both families of directed clock edges keep their orientations. Thus \(P\) followed by complex conjugation takes the shell at the original angles to the shell at the two negated angles. It carries the three softened spokes to the reflected orbit. Each possible perturbation has operator norm at most \(6\delta\), so the same transfer bounds apply to both orbits.

The two operations each interchange the nontrivial cyclic characters. Their composition preserves each character and transports a normalized lowest vector up to a unit phase. On the readout this gives complex conjugation, a permutation of surviving rows, and unit phases on columns. These operations preserve singular values. The set of omitted layouts is closed under reflection. It suffices to cover 128 closed rows of the 256-by-256 angle grid. The replay verifies the permutation, edge distances, ordered edge phase numerators and clock orientations.

## Transfer and coverage

The degree-14 Laurent trial vectors are proposed numerically. A 64-by-64 complex-ball transform encloses the degree-28 residual coefficients without aliasing. The uniform reference gap, the trial norm bounds \(0.4<\lVert v\rVert<1.5\), and the bound \(6\delta\) control the projection onto the softened ground line. The five residual columns in each cyclic character are combined with a Frobenius bound. The three character spaces are orthogonal, and deleting rows is a contraction.

Cubic interpolation uses four equidistant nodes on each axis. The normalized nodal product has maximum modulus \(1/81\), and the single-axis interpolation operator has norm at most \(7/4\). Absolute weighted Fourier coefficients bound the two fourth derivatives. The integer interpolation first rounds on the \(2^{-42}\) grid and then converts the controls to the \(2^{-30}\) grid. The combined transfer at depth zero is below \(9.83758294689\times10^{-7}\). Each dyadic refinement adds at most six integer units per real component, with the Frobenius bound \(\sqrt{480}<22\).

For every retained angular box an integer left inverse candidate \(L\) satisfies \(\lVert I-LB\rVert_F+\tau\lVert L\rVert_F<1\) on each of the sixteen Bernstein controls \(B\). Convexity gives the same inequality at all angles in the box. Applying it to a unit input yields a lower bound \(\tau\) for the interpolant's smallest singular value. The code subtracts the depth-dependent transfer and divides by \(1.5\) to account for normalization of the projected columns. The remaining bound is greater than \(5\times10^{-7}\). The ledger records 11,877,939 accepted leaves and zero unresolved leaves. For each angular row, the exact dyadic area sum is 256 cells times 90 mask classes.

## Replay

Use a checkout of the shell-model source at commit cb58c8e64ecea67184e0c1d68a690be098a7a384. The old shell sources can share the same repository as the current replay through a detached worktree. Install NumPy and python-flint with eigenvalue enclosure support. The replay was run with python-flint 0.9.0.

~~~sh
git worktree add --detach ../shell-model cb58c8e64ecea67184e0c1d68a690be098a7a384
replay=experiments/ch39-omega-research/replay
model=../shell-model/experiments/ch39-omega-research/replay
nets=two-orbit-control-nets
python3 "$replay/ch39-pole-two-orbit-torus.py" "$model" "$nets" --generate-only
python3 "$replay/ch39-pole-two-orbit-torus.py" "$model" "$nets"
~~~

To reproduce the faster exact-integer path, compile the optional kernel to a shared library and pass it with --native. A portable NumPy integer calculation is available through the invocation above. For a sharded run, use disjoint --start and --stop ranges covering rows 0 through 127. Every row must be present with zero unresolved leaves.

~~~sh
c++ -std=c++17 -O3 -shared -fPIC "$replay/ch39-pole-two-orbit-torus-kernel.cpp" -o two-orbit-kernel.so
python3 "$replay/ch39-pole-two-orbit-torus.py" "$model" "$nets" --native two-orbit-kernel.so
~~~

[Shell clock construction](https://github.com/the-omega-institute/trureturing-experiments/blob/cb58c8e64ecea67184e0c1d68a690be098a7a384/experiments/ch39-omega-research/replay/ch39-native-light-projector-pump.py) · [Reference residual check](https://github.com/the-omega-institute/trureturing-experiments/blob/1499a01569f845d3230d2120ad8c4c3a347b20c7/experiments/ch39-omega-research/replay/ch39-five-charge-full-torus.py) · [Coverage ledger](ch39-pole-two-orbit-torus.json)
