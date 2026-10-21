"""Check for Question 5 of arXiv:2504.02795v4 (Question 26 in v1): primitives of rank d - 1 in the greedy convolution.

The greedy trees T_d(n) of Definition 7 are built by processing n = 1, 2, ...: n joins the branch with the smallest
base k whose next vertex is n and whose length stays at most d, otherwise it starts a new branch.  The rank of a
primitive a is the final length of its branch (Lemma 10 (iii)); it is final once d*a has been processed.
The script reproduces the source's matrix A(4,20) (Example 2), lists the ranks that occur for d <= 20, and confirms
that rank d - 1 does not occur; it also tests Conjecture 1's first clause for d = 4 (rank-1 and rank-2 primitives
divisible by 8 or 9).
"""


def greedy(d, N):
    length = {}
    for n in range(1, N + 1):
        best = None
        for a in sorted(length):
            if a >= n:
                break
            m = length[a]
            if m < d and n == (m + 1) * a:
                best = a
                break
        if best is None:
            length[n] = 1
        else:
            length[best] += 1
    return length


ok = True
L4 = greedy(4, 200)
prim = [a for a in sorted(L4) if a <= 20]
src_row = [1, 5, 6, 7, 8, 9, 11, 13, 17, 19]
src_rank = [4, 4, 4, 4, 2, 1, 4, 4, 4, 4]
print(f"d=4 primitives <= 20: {prim}, ranks {[L4[a] for a in prim]}; source Example 2 A(4,20): {src_row}, ranks {src_rank}")
ok &= prim == src_row and [L4[a] for a in prim] == src_rank
for d in range(1, 21):
    N = 4000
    L = greedy(d, N)
    ranks = sorted({l for a, l in L.items() if d * a <= N})
    has = (d - 1) in ranks if d >= 2 else False
    print(f"d={d}: ranks of primitives a <= {N // d}: {ranks}; rank d-1 occurs: {has}")
    ok &= not has
L4 = greedy(4, 40000)
bad = [a for a, l in L4.items() if 4 * a <= 40000 and l in (1, 2) and a % 8 and a % 9]
print(f"d=4: rank-1 and rank-2 primitives a <= 10000 not divisible by 8 or 9: {len(bad)}")
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
