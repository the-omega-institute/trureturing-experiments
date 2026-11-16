# Polujan–Pott vectorial development check

`check.py` uses exact integer arithmetic in the field \(\mathbb F_4=\mathbb F_2[t]/(t^2+t+1)\). It exhaustively checks the proposed identity for every pair of points in \(K^6\), verifies the Maiorana–McFarland Walsh orthogonality calculation for every nonzero output multiplier and all linear characters, and checks degree witnesses. The output is in `check.out`; the final line is `ALL_CHECKS_OK`.
