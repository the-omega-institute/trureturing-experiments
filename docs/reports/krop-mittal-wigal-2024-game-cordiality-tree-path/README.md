# Krop–Mittal–Wigal: game cordiality of trees versus paths

`check.py` evaluates the cordiality game of arXiv:2403.18060v1 exactly: Admirable starts and labels 0, Impish labels 1; an edge gets the sum of its endpoint labels mod 2; the value c_g is the minimax value of |e_1 − e_0|. It refutes the conjecture of Section 3, "For any tree T of order n, c_g(T) ≤ c_g(P_n)".

- c_g(P_10) = 1 and c_g(T) = 3 for the subdivided claw T with edges (i, i+1) for i < 7 and (0, 8), (0, 9) (arm lengths 7, 1, 1).
- Impish's pairing strategy (0,4), (1,5), (2,6), (3,7), (8,9) alone forces |e_1 − e_0| ≥ 3 against every Admirable play.
- The family T_r (2r leaves attached at vertex 0 of P_8) has c_g(T_r) = 3 > c_g(P_{8+2r}) = 1 for r = 1, 2.

Run: `python3 check.py` (exit 0; output lines `cg(P10)= 1 cg(T)= 3`, `pairing lower bound on T = 3`, `r=1 n=10 cg(T_r)=3 cg(P_n)=1`, `r=2 n=12 cg(T_r)=3 cg(P_n)=1`, `bad= 0`).
