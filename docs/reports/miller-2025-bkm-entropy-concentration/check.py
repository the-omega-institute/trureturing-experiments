"""Numerical check for the entropy-concentration conjecture of H. J. D. Miller, "Entropy-based random
quantum states", arXiv:2511.01988v2 (Section VI, last paragraph: "A rigorous proof of this concentration
of measure is left as an open conjecture").

Sampler (source Eq. (41), Appendix B): lower-triangular X with independent standard circular complex
Gaussian off-diagonal entries (real and imaginary variance 1/2) and diagonal e^{i phi_k}|Z_k|/sqrt 2;
rho = X X* / Tr(X X*) (the Haar factor U is omitted: entropy is unitarily invariant).
Proposed theorem (search seat): Var S(rho_N) <= pi^2 B_N / (2 (N^2 - 2)), B_N = (ln N)^2 + 2 ln N + 2.

Prediction registered before the run:
  (1) the sample variance of S is below the proposed bound for N = 4, 8, 16, 32, 64, 128 (with 2000
      samples; the bound is not tight, so the ratio should be well below 1) and decreases with N;
  (2) the sample mean of S agrees with the source's mean formula ln N - gamma - ln 2 + 1/2 up to
      O(ln N / N);
  (3) removing the diagonal phases does not change the entropy distribution (two-sample KS statistic
      below 0.05 at N = 16 with 4000 samples each);
  (4) the spectral lemma sum p_i (ln p_i)^2 <= B_N holds on every sampled spectrum and on the uniform
      spectrum (where it is (ln N)^2).
Negative control (version 1, failed as registered and kept in run-v1.log): one tenth of the bound was
predicted to be violated at N = 4; the bound exceeds the sample variance by a factor of about 110 there, so
that control did not fail.  It measured the looseness of the bound, not the power of the test.
Negative control (version 2, registered before its first run): a law that does not concentrate - the
pure state with probability 1/2 and the maximally mixed state otherwise - must have tail probability
P(|S - E S| > 1/2) = 1 for N >= 8, so the tail test used for the claim detects non-concentration.
Usage: python3 check.py   (numpy; seed fixed)
"""
import sys
import numpy as np

rng = np.random.default_rng(20261017)
OK = True


def fail(m):
    global OK
    OK = False
    print("FAIL", m)


def sample(N, phases=True):
    X = np.zeros((N, N), complex)
    il = np.tril_indices(N, -1)
    X[il] = (rng.normal(size=len(il[0])) + 1j * rng.normal(size=len(il[0]))) / np.sqrt(2)
    Z = rng.normal(size=N)
    ph = np.exp(2j * np.pi * rng.uniform(size=N)) if phases else np.ones(N)
    X[np.diag_indices(N)] = ph * np.abs(Z) / np.sqrt(2)
    G = X @ X.conj().T
    p = np.linalg.eigvalsh(G)
    p = np.clip(p, 0, None)
    p = p / p.sum()
    return p


def entropy(p):
    q = p[p > 0]
    return -(q * np.log(q)).sum()


def BN(N):
    return np.log(N) ** 2 + 2 * np.log(N) + 2


print("== (1),(2),(4) variance, mean and spectral lemma")
prev = None
gamma = 0.5772156649015329
for N in (4, 8, 16, 32, 64, 128):
    M = 2000 if N <= 64 else 800
    S = []
    lemma_ok = True
    for _ in range(M):
        p = sample(N)
        S.append(entropy(p))
        q = p[p > 0]
        if (q * np.log(q) ** 2).sum() > BN(N) + 1e-9:
            lemma_ok = False
    S = np.array(S)
    var = S.var(ddof=1)
    bound = np.pi ** 2 * BN(N) / (2 * (N * N - 2))
    mean_formula = np.log(N) - gamma - np.log(2) + 0.5
    print(f"  N={N:4d}: mean {S.mean():.4f} (formula {mean_formula:.4f}), var {var:.3e}, bound {bound:.3e}, ratio {var/bound:.3f}, lemma {lemma_ok}")
    if not var < bound:
        fail(f"variance bound N={N}")
    if prev is not None and not var < prev:
        fail(f"variance not decreasing at N={N}")
    if not lemma_ok:
        fail(f"spectral lemma N={N}")
    if not abs(S.mean() - mean_formula) < 3 * np.log(N) / N + 0.05:
        fail(f"mean formula N={N}")
    prev = var
if not abs(np.log(8) ** 2 - (np.ones(8) / 8 * np.log(np.ones(8) / 8) ** 2).sum()) < 1e-12:
    fail("uniform spectrum lemma")

print("== (3) diagonal phases do not change the entropy law (N = 16)")
a = np.sort([entropy(sample(16, True)) for _ in range(4000)])
b = np.sort([entropy(sample(16, False)) for _ in range(4000)])
xs = np.concatenate([a, b])
ks = max(abs(np.searchsorted(a, x, side="right") / len(a) - np.searchsorted(b, x, side="right") / len(b)) for x in xs)
print(f"  KS statistic {ks:.4f}")
if not ks < 0.05:
    fail("phase removal")

print("== tail probabilities P(|S - E S| > 1/2) for the BKM law")
for N in (8, 32, 128):
    S = np.array([entropy(sample(N)) for _ in range(1000 if N < 128 else 400)])
    print(f"  N={N}: {np.mean(abs(S - S.mean()) > 0.5):.4f}")
print("== negative control (version 2): a non-concentrating law")
for N in (8, 32, 128):
    S = np.where(rng.uniform(size=2000) < 0.5, 0.0, np.log(N))
    tail = np.mean(abs(S - S.mean()) > 0.5)
    print(f"  N={N}: tail {tail:.4f}")
    if not tail > 0.99:
        fail(f"control v2 N={N}")
print("ALL_OK" if OK else "SOME_FAILED")
sys.exit(0 if OK else 1)
