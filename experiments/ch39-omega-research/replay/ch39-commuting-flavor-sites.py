"""Reconstruct one nonzero-J flavor pair with five, six and seven physical sites."""
import numpy as np


def five_site_frame(a, b, floor):
    """Subtract positive parent floors, whiten the sum, and complete its frame."""
    eye = np.eye(3)
    s = a + b - 2 * floor * eye
    mass, basis = np.linalg.eigh(s)
    assert mass[0] > 0 and mass[1] < mass[2]
    scale = mass[-1]
    root = (basis * np.sqrt(mass)) @ basis.conj().T
    inverse = (basis * (1 / np.sqrt(mass))) @ basis.conj().T
    fractions, common = np.linalg.eigh(inverse @ (a - floor * eye) @ inverse)
    assert 0 < fractions[0] <= fractions[-1] < 1
    rows = common.conj().T @ root / np.sqrt(scale)
    slack = np.sqrt(1 - mass[:2] / scale)[:, None] * basis[:, :2].conj().T
    frame = np.vstack((rows, slack))
    da = np.r_[floor + scale * fractions, floor, floor]
    db = np.r_[floor + scale * (1 - fractions), floor, floor]
    return frame, da, db


def check(a, b, frame, da, db):
    assert np.all((0 < da) & (da < 1)) and np.all((0 < db) & (db < 1))
    pull = frame.conj().T
    errors = [np.linalg.norm(pull @ frame - np.eye(3)),
              np.linalg.norm(pull @ (da[:, None] * frame) - a),
              np.linalg.norm(pull @ (db[:, None] * frame) - b)]
    assert max(errors) < 3e-14
    print(frame.shape[0], "sites, largest reconstruction error", max(errors))


def main():
    indices = np.arange(3)
    mixing = np.exp(2j * np.pi * np.outer(indices, indices) / 3) / np.sqrt(3)
    up = np.diag([1., 2., 3.]) / 8
    down = mixing @ np.diag([1., 2., 4.]) @ mixing.conj().T / 8
    j = (mixing[0, 0] * mixing[1, 1] * (mixing[0, 1] * mixing[1, 0]).conjugate()).imag
    assert abs(j - 1 / (6 * np.sqrt(3))) < 1e-15
    print("signed J", j)
    frame, da, db = five_site_frame(up, down, floor=1/16)
    check(up, down, frame, da, db)

    # Three distinguished sites; the last three form a rotation orbit.
    f6 = np.vstack((frame, np.zeros((1, 3))))
    check(up, down, f6, np.r_[da, da[-1]], np.r_[db, db[-1]])

    # Four distinguished sites; the invariant vector on the orbit is the fifth row.
    f7 = np.vstack((frame[:4], np.tile(frame[4] / np.sqrt(3), (3, 1))))
    assert np.array_equal(np.roll(f7[4:], 1, axis=0), f7[4:])
    check(up, down, f7, np.r_[da, da[-1], da[-1]], np.r_[db, db[-1], db[-1]])


if __name__ == "__main__":
    main()
