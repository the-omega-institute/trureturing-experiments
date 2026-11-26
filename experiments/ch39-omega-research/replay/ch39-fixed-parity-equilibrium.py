"""Certify both fixed-physical-parity ground states of the 38-mode shell."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

model = run_path(str(Path(__file__).with_name("ch39-majorana-vacuum-energy.py")))
shell = model["shell"]


def vacuum_sector(data):
    bulk, ports = data
    coupling = ports * shell["source"].transpose()
    response = -coupling * bulk.solve(coupling.conjugate().transpose())
    pfaffian = bulk.det().real * response[0, 1].imag
    assert 0 not in pfaffian
    return 1 if pfaffian > 0 else -1


def verify():
    winners = {1: 11177, -1: 11618}
    margins = {1: arb(".00015"), -1: arb(".000028")}
    thresholds = {}
    for sector, index in winners.items():
        height = acb.zeta_zero(index).imag
        assert vacuum_sector(shell["parent_data"](height)) == sector
        thresholds[sector] = model["vacuum"](model["hamiltonian"](height)) + margins[sector]

    for index, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        energy = model["vacuum"](model["hamiltonian"](zero.imag))
        if all(energy > cutoff for cutoff in thresholds.values()):
            continue
        data = shell["parent_data"](zero.imag)
        vacuum_sign = vacuum_sector(data)
        for sector, cutoff in thresholds.items():
            if index == winners[sector] or energy > cutoff:
                continue
            assert sector != vacuum_sign
            # The zero of the decreasing upper Schur eigenvalue is the positive gap.
            required_gap = (cutoff - energy).upper()
            assert 0 <= required_gap < arb(".025")
            assert shell["sign_at"](data, required_gap)

    # Raising k bulk modes in a fixed sector costs at least k*(1-.025).
    beta = arb(12)
    thermal_error = 18 * (-arb(".975") * beta).exp() / beta
    assert thermal_error < min(margins.values())
    print("Fixed Pfaffian orientations + and - select 11177 and 11618.")
    print("Certified rival margins > 0.00015 and 0.000028; beta >= 12.")


if __name__ == "__main__":
    verify()
