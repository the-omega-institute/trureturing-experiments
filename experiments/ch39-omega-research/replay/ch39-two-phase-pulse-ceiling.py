"""Enclose both two-tick transmission scores on the whole two-angle torus."""
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx

pulse = run_path(str(Path(__file__).with_name("ch39-finite-port-pulses.py")))
imaginary = acb(0, 1)
ports = acb_mat([[int(row == 38+col) for col in range(2)] for row in range(40)])


def generator_symbol():
    model = pulse["model"]
    def evaluate(z, w):
        blocks = model["shell"]["parent_blocks"]({724: z, 804: w})
        return pulse["attach_ports"](model["paired_hamiltonian"](*blocks))
    constant = evaluate(0, 0)
    symbol = {(0, 0): constant}
    for mode in ((1, 0), (0, 1)):
        cosine = evaluate(*mode)-constant
        sine = evaluate(*(imaginary*x for x in mode))-constant
        symbol[mode] = (cosine-imaginary*sine)/2
        symbol[tuple(-x for x in mode)] = (cosine+imaginary*sine)/2
    return symbol


def trial_polynomials(symbol, degree):
    modes = {k: np.array(v.tolist(), dtype=complex) for k, v in symbol.items()}
    size = 2*degree+1
    phase = np.exp(2j*np.pi*np.arange(size)/size)
    unit = np.eye(40)
    source = np.array(ports.tolist(), dtype=complex)
    samples = np.empty((2, size, size, 40, 2), dtype=complex)
    for i, j in product(range(size), repeat=2):
        g = sum(v*phase[i]**k*phase[j]**l for (k,l),v in modes.items())
        matrix = unit+1j*g
        first = 2*np.linalg.solve(matrix, source)-source
        retained = first.copy()
        retained[38:] = 0
        second = 2*np.linalg.solve(matrix, retained)-retained
        samples[:, i, j] = first, second
    coeffs = np.fft.fft2(samples, axes=(1,2))/size**2
    return [{(k,l): acb_mat(array[k,l].tolist())
             for k,l in product(range(-degree,degree+1), repeat=2)} for array in coeffs]


def residual_bound(symbol, before, after):
    zero = acb_mat(40, 2)
    residual = {}
    for k,l in before.keys() | after.keys():
        old, new = before.get((k,l), zero), after.get((k,l), zero)
        residual[k,l] = residual.get((k,l), zero)+(new-old)
        for (u,v), g in symbol.items():
            tag = (k+u,l+v)
            residual[tag] = residual.get(tag, zero)+imaginary*(g*(new+old))
    return sum(pulse["model"]["frobenius"](value) for value in residual.values())


def amplitude_ceiling(coeffs, amplitude_error, size):
    # Bilinear interpolation error is at most h^2/8 times the sum of the
    # pure second derivative bounds. Convex interpolation preserves a norm cap.
    curvature = sum((k*k+l*l)*abs(v) for (k,l),v in coeffs.items())
    interpolation_error = (2*arb.pi()/size)**2*curvature/8
    columns = []
    for k in sorted({k for k,l in coeffs}):
        row = [acb(0)]*size
        for (a,l),v in coeffs.items():
            if a == k:
                row[l % size] = v
        columns.append((k,acb.dft(row)))
    upper = arb(0)
    for j in range(size):
        row = [acb(0)]*size
        for k,col in columns:
            row[k % size] = col[j]
        upper = max(upper, max(abs(v).upper() for v in acb.dft(row)))
    return (upper + interpolation_error + amplitude_error)**2


def verify():
    ctx.prec = 144
    symbol = generator_symbol()
    first, second = trial_polynomials(symbol, 12)
    first_error = residual_bound(symbol, {(0,0): ports}, first)
    retained = {k: acb_mat(v) for k,v in first.items()}
    for v in retained.values():
        for i,j in product((38,39),range(2)):
            v[i,j] = 0
    second_error = residual_bound(symbol, retained, second)+first_error
    error = (7*first_error+24*second_error)/25
    print("uniform pulse amplitude error", error.upper().str(8), flush=True)
    limit = arb(".00742951")
    for row,col in ((39,0),(38,1)):
        coeffs = {k:(-7*v[row,col]+24*second[k][row,col])/25 for k,v in first.items()}
        bound = amplitude_ceiling(coeffs, error, 512)
        print("port score ceiling", bound.upper().str(10), flush=True)
        assert bound < limit
        if row==38:
            assert bound < arb(".007300684")
    return limit


if __name__ == "__main__":
    verify()
