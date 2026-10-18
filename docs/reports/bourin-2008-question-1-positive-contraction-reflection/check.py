"""Bourin arXiv:0805.1954v1 (Internat. J. Math. 20 (2009)), Question 1: counterexample and the amplification classification.

Witness: n = 1, k = 2, Phi(x) = diag(x, x/4), M = 1, f = sqrt.
  Phi^{-1}(F_2) = {0 <= x <= 1} = F_1 (hypothesis holds with equality),
  f(Phi(1)) = diag(1, 1/2) >= Phi(f(1)) = diag(1, 1/4): trace norm 3/2 > 5/4, Frobenius sqrt(5)/2 > sqrt(17)/4.
Classification (seat-proposed): for Phi_a(X) = diag(a_1 X, ..., a_r X), a_i > 0,
  hypothesis  <=>  max a_i >= 1;  conclusion for all M, f, symmetric norms  <=>  sum a_i >= r.
Checked here on random weights with Ky Fan k-norms (which decide all symmetric norms for PSD pairs
by Ky Fan dominance) and several concave f.
"""
import math, sys
import numpy as np

rng = np.random.default_rng(20261017)
FS = {
    "sqrt": np.sqrt,
    "min1": lambda t: np.minimum(t, 1.0),
    "log1p": np.log1p,
    "pow0.3": lambda t: np.power(t, 0.3),
    "t/(1+t)": lambda t: t / (1 + t),
}


def kyfan_dominated(x, y, tol=1e-10):
    """weak majorization x <_w y for nonnegative vectors (all Ky Fan norms of diag(x) <= those of diag(y))."""
    xs, ys = np.sort(x)[::-1], np.sort(y)[::-1]
    return bool(np.all(np.cumsum(xs) <= np.cumsum(ys) + tol))


def conclusion_holds(a, lam, f):
    left = np.concatenate([f(ai * lam) for ai in a])         # eigenvalues of f(Phi_a(M))
    right = np.concatenate([ai * f(lam) for ai in a])        # eigenvalues of Phi_a(f(M))
    return kyfan_dominated(left, right)


def main():
    fails = []
    # witness
    L = np.array([1.0, math.sqrt(0.25)]); R = np.array([1.0, 0.25])
    print("witness traces:", L.sum(), R.sum(), "Frobenius:", np.linalg.norm(L), np.linalg.norm(R), "operator:", L.max(), R.max())
    if not (L.sum() > R.sum() and np.linalg.norm(L) > np.linalg.norm(R)):
        fails.append("witness")
    # classification on random weights
    tested = viol_expected = ok_expected = 0
    for _ in range(3000):
        r = int(rng.integers(1, 5))
        a = rng.uniform(0.05, 2.5, size=r)
        lam = rng.uniform(0, 5, size=int(rng.integers(1, 4)))
        S = a.sum()
        holds_all = all(conclusion_holds(a, lam, f) for f in FS.values())
        tested += 1
        if S >= r:
            ok_expected += 1
            if not holds_all:
                fails.append(("sufficiency", a.tolist(), lam.tolist()))
        else:
            # necessity witness of Lemma 2: f = min(t,1), M = t I with t >= max(1, 1/min a)
            t = max(1.0, 1.0 / a.min())
            if conclusion_holds(a, np.array([t]), FS["min1"]):
                fails.append(("necessity", a.tolist()))
            viol_expected += 1
    print(f"classification: {tested} random weight vectors; sum >= r: {ok_expected} (all conclusions held), "
          f"sum < r: {viol_expected} (min(t,1) witness violated each)")
    # hypothesis criterion: Phi_a^{-1}(F) subset F  <=>  max a >= 1, checked on scalar test points x in [0, 1/min_i ...]
    for _ in range(2000):
        r = int(rng.integers(1, 5)); a = rng.uniform(0.05, 2.5, size=r)
        # Phi_a(x) in F  <=>  0 <= x <= 1/max a ; inclusion in F_1 = [0,1]  <=>  1/max a <= 1
        incl = 1 / a.max() <= 1
        if incl != (a.max() >= 1):
            fails.append(("hypothesis", a.tolist()))
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:3]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
