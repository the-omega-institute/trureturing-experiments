# Çevik Problems 10.1 and 10.2: finite dimensional cover separation

Run `python3 check.py` with Python 3.12 or newer; no third party packages are needed. The fidelity check reproduces the two-majorant meet construction of the pyramid cone in [arXiv:2611.12376v1](https://arxiv.org/abs/2611.12376v1), Example 6.4, on 625 integer inputs. The other checks use exact rational arithmetic: they compare the piecewise affine spikes at their breakpoints, verify selected quadratic majorants on the whole interval, and bound positive polynomial tests of degrees 2 through 8 by exact interpolation. The degree bounds tend to zero with the spike index even though the spike values at zero remain above one.

The output is finite corroboration. The argument for all finite dimensional extensions, the order infimum in the continuous function lattice, and the claims about every vector metric use compactness and continuity.
