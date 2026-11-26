"""Bound all sufficiently long coherent pulses and the zero-energy contact threshold."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

shell = run_path(str(Path(__file__).with_name("ch39-majorana-portal-gap.py")))


def verify():
    values = {}
    for index, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        bulk, contacts = shell["parent_data"](zero.imag)
        coupling = contacts * shell["source"].transpose()
        lifted = bulk.solve(coupling.conjugate().transpose())
        assert sum(abs(z)**2 for z in lifted.entries()) < arb(".127")
        values[index] = abs((coupling*lifted)[0, 1].imag)

    peak, second = values[11431], values[10965]
    assert peak > second
    assert all(second > v for n, v in values.items() if n not in (11431, 10965))
    threshold = (peak*second).sqrt()
    width = abs(shell["columns"].det())
    assert width > threshold
    def amplitude(y):
        return 2*width*y/(width**2+y**2)
    delta = amplitude(peak)-amplitude(second)
    assert amplitude(peak)**2-amplitude(second)**2 > arb(".00004102")
    gram = shell["columns"].transpose()*shell["columns"]
    # The two columns are normalized shell normals.
    lower = 1 - abs(gram[0, 1])
    stored_bound = arb("1.127")/lower
    assert stored_bound < arb("6.8")
    assert delta > arb("27.2") / (2**17)
    print("Contact determinant threshold", threshold.str(14), "native value", width.str(14))
    print("Native steady probability", (amplitude(peak)**2).str(14),
          "amplitude lead", delta.str(12), "winner 11431 at all lengths >= 131072")


if __name__ == "__main__":
    verify()
