# ddt2610bs0 exact check

`check.py` exhaustively enumerates the eight strictly upper triangular 3 by 3 matrices over the exact field `F₂`. It verifies the source's bounded-nilexponent identity `x³ = 0` and `x² y x² = 0`, constructs the source family `Y₃ = {b x² : b ∈ R¹, x ∈ R}`, checks every member is second-order, and finds the first vanishing product length for `Y₃`. It also checks a finite right `BS₀` compression bound and the nilpotence class of the example. All calculations are integer arithmetic modulo 2; no floating point is used.
