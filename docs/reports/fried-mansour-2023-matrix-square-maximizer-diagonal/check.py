"""Fried–Mansour arXiv:2308.00348v1, Section 2.2, condition (d): if s(A^2) is maximal among n x n matrices
whose entries are 1, ..., n^2, then 1 and n^2 lie on the main diagonal.

1. n = 2, 3: all (n^2)! labelings; every maximizer has 1 and n^2 on the diagonal, and after the simultaneous
   row-and-column permutation sorting the row sums, every row and every column is strictly increasing
   (the stronger statement proved in the issue; it is H. Pfoertner's observation in OEIS A368539).
   The maximum values are A368539(2), A368539(3) = 54, 761.
2. n = 4, ..., 7: restarted hill climbing with entry transpositions reaches the A368539 values
   5284, 24303, 85352, 248045; every maximizer found satisfies both properties.
3. Identity s(A^2) = sum_k r_k c_k on random matrices; tie-freeness of row sums at maximizers.
4. Control: non-maximal labelings can have 1 or n^2 off the diagonal (most random labelings do).
"""
import itertools, random, sys

random.seed(20261017)


def s2(A):
    n = len(A)
    r = [sum(A[k]) for k in range(n)]
    c = [sum(A[i][k] for i in range(n)) for k in range(n)]
    return sum(r[k] * c[k] for k in range(n))


def s2_direct(A):
    n = len(A)
    return sum(A[i][k] * A[k][j] for i in range(n) for k in range(n) for j in range(n))


def props(A):
    n = len(A)
    diag = {A[i][i] for i in range(n)}
    d_ok = 1 in diag and n * n in diag
    r = [sum(A[k]) for k in range(n)]
    if len(set(r)) < n: return d_ok, False
    p = sorted(range(n), key=lambda k: r[k])
    B = [[A[p[i]][p[j]] for j in range(n)] for i in range(n)]
    mono = all(B[i][j] < B[i][j + 1] for i in range(n) for j in range(n - 1)) and \
           all(B[i][j] < B[i + 1][j] for i in range(n - 1) for j in range(n))
    return d_ok, mono


def exhaustive(n):
    best, arg = -1, []
    for perm in itertools.permutations(range(1, n * n + 1)):
        A = [list(perm[i * n:(i + 1) * n]) for i in range(n)]
        v = s2(A)
        if v > best: best, arg = v, [A]
        elif v == best: arg.append(A)
    return best, arg


def climb(n, restarts, iters):
    best, found = -1, {}
    cells = [(i, j) for i in range(n) for j in range(n)]
    for _ in range(restarts):
        vals = list(range(1, n * n + 1)); random.shuffle(vals)
        A = [vals[i * n:(i + 1) * n] for i in range(n)]
        v = s2(A); stall = 0
        while stall < iters:
            (a, b), (c, d) = random.sample(cells, 2)
            A[a][b], A[c][d] = A[c][d], A[a][b]
            w = s2(A)
            if w >= v:
                stall = 0 if w > v else stall + 1; v = w
            else:
                A[a][b], A[c][d] = A[c][d], A[a][b]; stall += 1
        if v > best: best, found = v, {}
        if v == best: found[tuple(map(tuple, A))] = True
    return best, [list(map(list, k)) for k in found]


def main():
    fails = []
    for _ in range(200):
        n = random.randint(1, 6); A = [[random.randint(-9, 9) for _ in range(n)] for _ in range(n)]
        if s2(A) != s2_direct(A): fails.append("identity")
    print("identity s(A^2) = sum r_k c_k: checked on 200 random integer matrices")
    target = {2: 54, 3: 761, 4: 5284, 5: 24303, 6: 85352, 7: 248045}
    for n in (2, 3):
        best, arg = exhaustive(n)
        ok = [props(A) for A in arg]
        print(f"n={n}: max s(A^2) = {best} (A368539: {target[n]}), maximizers = {len(arg)}, "
              f"(d) holds for {sum(o[0] for o in ok)}, monotone after sorting for {sum(o[1] for o in ok)}")
        if best != target[n] or not all(o[0] and o[1] for o in ok): fails.append(("exhaustive", n))
    for n, (R, I) in {4: (60, 3000), 5: (60, 6000), 6: (40, 12000), 7: (30, 20000)}.items():
        best, arg = climb(n, R, I)
        ok = [props(A) for A in arg]
        print(f"n={n}: best s(A^2) found = {best} (A368539 bound: {target[n]}), distinct maximizers found = {len(arg)}, "
              f"(d) {sum(o[0] for o in ok)}/{len(arg)}, monotone {sum(o[1] for o in ok)}/{len(arg)}")
        if best == target[n] and not all(o[0] and o[1] for o in ok): fails.append(("climb", n))
        if best > target[n]: fails.append(("exceeds-oeis", n))
    # control: random labelings usually have 1 or n^2 off the diagonal
    n = 4; off = 0
    for _ in range(1000):
        vals = list(range(1, 17)); random.shuffle(vals); A = [vals[i * 4:(i + 1) * 4] for i in range(4)]
        off += not props(A)[0]
    print(f"control: {off}/1000 random 4x4 labelings have 1 or 16 off the diagonal")
    if off == 0: fails.append("control")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
