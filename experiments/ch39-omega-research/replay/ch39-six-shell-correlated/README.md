# Correlated six-shell interval calibration

The model fixes the fifth-shell charge-six source and four rounds of directed chords. The second-round chord $(0,2)$ and fourth-round chord $(6,13)$ each get an independent gate. The six gates meet the native third-pivot targets $K_j=1$ and $W_j=10^{-122}$ simultaneously for shells 7 through 12. The run fixes 12 of the 24 real angles and solves the 12 residual equations for the other angles. [Original inputs, source, and derivation](https://github.com/the-omega-institute/trureturing-books/issues/1743#issuecomment-7235856614)

`save-native-intervals.py` regenerates the exact dyadic interval inputs. In `certify-correlated-shells.py`, `verify` preconditions the residual before taking absolute values. This preserves the common-source cancellation across the six shells. The second derivative estimate uses all 78 unordered angle pairs, including diagonal terms. The interval replay reproduced a maximum half-width below $2.452\times10^{-7}$, a derivative defect below $0.00395$, and a self-map ratio below $0.91301$. The resulting contraction gives a unique zero inside the verified box for the fixed angles.

Run in a writable copy of this directory. Extract the [pinned native replay ZIP](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/ch39-omega-math-replay.zip) into a child directory named `pinned`. Place the [pinned final-sign script](https://github.com/the-omega-institute/trureturing-experiments/blob/47db2f76c022b3afbcc128f31084608ac113f927/experiments/ch39-omega-research/replay/ch39-three-shell-final-sign.py) alongside `candidate.json`. With Python, `python-flint` 0.9.0, NumPy, and SymPy installed, run with assertions enabled.

```sh
OPENBLAS_NUM_THREADS=1 python save-native-intervals.py
OPENBLAS_NUM_THREADS=1 python certify-correlated-shells.py
```
