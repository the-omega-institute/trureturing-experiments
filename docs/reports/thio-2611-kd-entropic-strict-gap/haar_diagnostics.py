"""Floating-point exploration only; no float is used as mathematical evidence."""
import numpy as np


def entropy(p):
    p = p[p > 0]
    return -np.dot(p, np.log(p))


def objective(U, psi):
    return entropy(abs(psi)**2) + entropy(abs(U.conj().T @ psi)**2)


def optimize(U, psi):
    """Find a feasible upper bound to the pure minimum, without claiming global optimality."""
    for _ in range(350):
        b = U.conj().T @ psi
        p, q = abs(psi)**2, abs(b)**2
        la = np.zeros_like(p); lb = np.zeros_like(q)
        np.log(p, out=la, where=p > 0); np.log(q, out=lb, where=q > 0)
        g = -psi * (la + 1) - U @ (b * (lb + 1))
        g -= psi * np.vdot(psi, g).real
        old = objective(U, psi)
        step = 0.25
        for _ in range(35):
            trial = psi - step * g
            trial /= np.linalg.norm(trial)
            if objective(U, trial) < old - 1e-14:
                psi = trial
                break
            step /= 2
        else:
            break
    return objective(U, psi)


def main():
    print('FLOAT_DIAGNOSTICS_ONLY (not certificates; numerical Haar QR and local optimization)')
    rng = np.random.default_rng(261104009)
    for d in (2, 3, 4):
        for case in range(3):
            X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
            U, R = np.linalg.qr(X)
            diagonal = np.diag(R)
            U = U @ np.diag(diagonal / abs(diagonal))
            laws = [abs(U[i, :])**2 for i in range(d)] + [abs(U[:, j])**2 for j in range(d)]
            entropies = [entropy(p) for p in laws]
            idx = int(np.argmin(entropies)); h = entropies[idx]
            W = U if idx < d else U.conj().T
            i = idx % d
            p = abs(W[i, :])**2
            logs = np.log(p); m = np.dot(p, logs)
            V = np.dot(p, (logs - m)**2)
            a = np.eye(d, dtype=complex)[:, i]
            v = W @ (logs * (W.conj().T @ a)) - m * a
            t = 1e-6
            psi = (a + t*v) / np.sqrt(1 + t*t*V)
            fd = (entropy(abs(W.conj().T @ psi)**2) - entropy(p)) / t
            seeds = [np.eye(d, dtype=complex)[:, k] for k in range(d)]
            seeds += [U[:, k].copy() for k in range(d)]
            for _ in range(3):
                z = rng.normal(size=d) + 1j*rng.normal(size=d)
                seeds.append(z/np.linalg.norm(z))
            upper = min(optimize(U, z) for z in seeds)
            assert upper < h - 1e-9
            assert abs(fd + 2*V) < 1e-4
            print(f'FLOAT d={d} case={case} vertex_h={h:.12f} feasible_pure_F={upper:.12f} '
                  f'fd_B={fd:.9f} minus2V={-2*V:.9f}')
