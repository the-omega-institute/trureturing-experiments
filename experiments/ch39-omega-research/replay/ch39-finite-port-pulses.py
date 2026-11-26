"""Certify finite directional pulse measurements on the particle-hole shell."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

model = run_path(str(Path(__file__).with_name("ch39-majorana-vacuum-energy.py")))
contacts = model["shell"]["columns"]
identity = np.eye(40)
imaginary = acb(0, 1)


def generator(height):
    return attach_ports(model["hamiltonian"](height))


def attach_ports(shell):
    coupled = acb_mat(40, 40)
    for i in range(38):
        for j in range(38):
            coupled[i, j] = shell[i, j]
    for i in range(2):
        for j in range(2):
            coupled[i, 38+j] = coupled[38+j, i] = contacts[i, j]
    return coupled


def port_outputs(coupled, length, *, impulse=False):
    matrix = np.array(coupled.tolist(), dtype=complex)
    cayley = np.linalg.solve(identity + 1j*matrix, identity - 1j*matrix)
    state = np.zeros((40, 2), dtype=complex)
    accumulated_error = arb(0)
    for tick in range(length):
        state[38:] = np.eye(2) if tick==0 or not impulse else 0
        before = acb_mat(state.tolist())
        state = cayley @ state
        after = acb_mat(state.tolist())
        # The inverse of I+iG is a contraction for Hermitian G.
        residual = after-before + imaginary*(coupled*(after+before))
        accumulated_error += model["frobenius"](residual)
        radius = arb(0, accumulated_error.upper())
        box = acb(radius, radius)
        yield after[39, 0]+box, after[38, 1]+box


def pulse_records(coupled):
    amplitudes = [acb(0), acb(0)]
    forward_intensity = arb(0)
    for length,(forward,reverse) in enumerate(port_outputs(coupled,8),1):
        amplitudes[0] += forward
        amplitudes[1] += reverse
        forward_intensity += abs(forward)**2
        p, q = [abs(a/length)**2 for a in amplitudes]
        yield p, q, (p+q)/2, forward_intensity/length


def verify():
    forward_winners = (11143, 11959, 11831, 11431, 11831, 11431, 11431, 11431)
    checks = {(i, 0): n for i, n in enumerate(forward_winners, 1)}
    checks.update({(4, 1): 10965, (4, 2): 11270, (4, 3): 11831})
    values = {key: {} for key in checks}
    for n, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        for length, record in enumerate(pulse_records(generator(zero.imag)), 1):
            for channel, probability in enumerate(record):
                if (length, channel) in values:
                    values[length, channel][n] = probability
    for key, winner in checks.items():
        scores = values[key]
        margin = scores[winner] - max(p.upper() for n, p in scores.items() if n != winner)
        assert margin > 0
        if key == (4, 0):
            assert margin > arb(".0000174")
        if key == (7, 0):
            assert margin > arb(".000412")
        print(key, winner, scores[winner].str(12), "margin >", float(margin.lower()))

    # One fixed-size head-to-head match per challenger in a single full pass.
    shots = 2**27
    failure = 1277 * (-shots * arb(".000412")**2).exp()
    assert failure < arb("0.000000163")
    print("1277 seven-tick comparisons, wrong final label probability <", failure.upper().str(6))


if __name__ == "__main__":
    verify()
