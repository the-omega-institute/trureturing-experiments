# Karam's fixed-centre ball-growth conjecture

Source: Steve Karam, “Ball growth at all scales in universal covers of short metric graphs”, arXiv:2610.20088v1, Conjecture D, Definition 1.1, Proposition 7.1 and Example 7.2.

Run `python3 check.py`. Python 3's standard library suffices. All evidence uses integers and `fractions.Fraction`; `output.txt` is the complete deterministic output. The program exits zero and ends with `ALL_CHECKS_OK`.

The ball evaluator unfolds the graph along nonbacktracking edge histories and sums the lengths of full and clipped lifted edges. It does not identify distinct lifts with the same projected vertex. Fidelity checks reproduce equilateral equality and the source's non-equilateral theta equality at the longest edge length.

The proposed counterexamples have vertices `(i, sign, j)` with `i` modulo 5, `sign` equal to ±1 and `j` modulo a positive integer `d`. Top edges form one cycle of length `5d` by advancing the sheet at the 4-to-0 seam. Bottom edges form separate 5-cycles, and rungs join the layers. The graph is connected, simple and cubic, with `10d` vertices, `15d` edges and Betti number `5d+1`.

Give top edge `i` length `scale * (1 + eps * a[i])`, bottom edge `i` length `scale * (1 - eps * a[i])`, and each rung length `scale`, where `a=(3,-8,10,-8,3)` and `0<eps<1/60`. Total length is `15d * scale`, exactly the conjecture's threshold. Every edge lies strictly between `5*scale/6` and `7*scale/6`.

The incident perturbation sums are `sign * t[i]`, with `t=(6,-5,2,2,-5)`. The second-radius coefficients are `sign * u[i]`, with `u=(-4,3,-1,-1,3)`. At radii `3*scale/2` and `5*scale/2`, the ball lengths are respectively `scale * (6 - sign*eps*t[i])` and `scale * (15 - sign*eps*u[i])`. Every product `t[i]*u[i]` is negative, so every lifted vertex falls strictly below the model at one radius.

For normalized edge lengths in `(5/6,7/6)`, one level is complete at radius `3/2` and two are complete at radius `5/2`; the next level has not finished. Summing those edge segments gives `6R-S(v)` and `12R-2S(v)-sum(S(w), w adjacent v)`. These identities and the local 5-periodic coefficient calculation prove the family statement for every positive `d` and every permitted real `eps`, independently of the finite checks. The graph projection to the `d=1` graph is locally bijective, so the local calculation is independent of the number of sheets.

The program checks the exact affine coefficients over the entire open perturbation interval for 1,260 vertices, and evaluates 11,340 vertex/parameter cases for `d=1..12,17,31`, three perturbations and three scales. These checks are independent arithmetic evidence; no Lean proof or kernel settlement is claimed.
