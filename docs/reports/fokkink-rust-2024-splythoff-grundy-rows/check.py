"""Check of the observation of Fokkink-Rust, 'Queen reflections: a modification of Wythoff Nim' (arXiv:1904.08339v2;
Int. J. Game Theory 53 (2024) 1099-1117), Section 'Some further remarks on Splythoff Nim': 'every row or column contains
each Sprague-Grundy value once and no more.'

Moves from (x, y): single (reduce one pile), double (reduce both by the same k >= 1), split (a double move that empties
the smaller pile, then split the remainder r = max - min into two piles (c, r - c)).  G is the mex of follower values.

(1) sg.c computes G on [0, N]^2 (here N = 900); the computed rows 0, 1, 15, 16 agree with the source's Table 6.
(2) For rows x <= 100 and values g <= 150, the value g occurs exactly once in row x within [0, N], at a column
    y <= 5x + 2g + 3 (the bound of the proposed counting proof).  G is symmetric and every row is injective on [0, N].
"""
import subprocess, sys, os
here = os.path.dirname(os.path.abspath(__file__))
bad = 0
subprocess.run(["cc", "-O2", "-o", os.path.join(here, "sg"), os.path.join(here, "sg.c")], check=True)
r = subprocess.run([os.path.join(here, "sg"), "900", "100", "150"], capture_output=True, text=True)
print(r.stdout.strip())
if r.returncode != 0: bad += 1
# Table 6 rows transcribed from the source (rows 16, 15, 1, 0; columns 0..17)
table = {16: [16,17,15,19,20,14,21,12,22,0,5,8,6,24,10,9,13,25],
         15: [15,16,14,18,19,17,20,21,12,2,4,22,23,7,8,11,9,24],
         1: [1,2,0,4,5,3,7,8,6,10,11,12,9,14,15,16,17,13],
         0: list(range(18))}
N = 40
G = [[0]*(N+1) for _ in range(N+1)]
for s in range(2*N+1):
    for x in range(max(0, s-N), min(s, N)+1):
        y = s-x; o = set()
        o |= {G[x-k][y] for k in range(1, x+1)}; o |= {G[x][y-k] for k in range(1, y+1)}
        m = min(x, y); o |= {G[x-k][y-k] for k in range(1, m+1)}
        if m > 0:
            rr = max(x, y) - m; o |= {G[c][rr-c] for c in range(rr+1)}
        g = 0
        while g in o: g += 1
        G[x][y] = g
for x, row in table.items():
    if [G[x][y] for y in range(18)] != row: bad += 1; print("table mismatch row", x)
os.remove(os.path.join(here, "sg"))
print("bad=", bad)
sys.exit(1 if bad else 0)
