"""Enclose the many-body vacuum energy of the full particle-hole shell at each candidate."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

shell = run_path(str(Path(__file__).with_name("ch39-majorana-portal-gap.py")))
unit = acb_mat([[int(i==j) for j in range(38)] for i in range(38)])


def hamiltonian(height):
    return paired_hamiltonian(*shell["parent_data"](height))


def paired_hamiltonian(m, f):
    f = (f*shell["source"].transpose())/arb(2).sqrt()
    h = acb_mat(38,38)
    for i in range(18):
        for j in range(18):
            h[2+i,2+j],h[20+i,20+j] = m[i,j],-m[i,j].conjugate()
        for a in range(2):
            h[a,2+i],h[a,20+i] = f[a,i],-f[a,i].conjugate()
            h[2+i,a],h[20+i,a] = f[a,i].conjugate(),-f[a,i]
    return h


def frobenius(m):
    return sum(abs(v)**2 for v in m.entries()).sqrt().upper()


def vacuum(h):
    # The floating eigensolver proposes a basis. All error bounds use ball arithmetic.
    _,approx = np.linalg.eigh(np.array(h.tolist(),dtype=complex))
    v = acb_mat(approx.tolist())
    gram = v.conjugate().transpose()*v
    v = v*(3*unit-gram)/2
    adj = v.conjugate().transpose()
    defect = frobenius(adj*v-unit)
    assert defect < arb("1e-20")
    hv = h*v
    ray = adj*hv
    diag = [ray[i,i].real.mid() for i in range(38)]
    residual = acb_mat([[hv[i,j]-v[i,j]*diag[j] for j in range(38)] for i in range(38)])
    # For the polar unitary Q, ||Q-V|| <= defect and ||V|| < 2.
    error = 2*frobenius(residual)+defect*(frobenius(h)+3*max(abs(z).upper() for z in diag))
    filled = sum(sorted(diag)[:19])
    return (filled+arb(0,19*error))/2


def verify():
    from flint import acb
    energies = {n:vacuum(hamiltonian(z.imag)) for n,z in enumerate(acb.zeta_zeros(10815,1278),10815)}
    first,second = energies[11618],energies[11619]
    assert all(second < e for n,e in energies.items() if n not in (11618,11619))
    assert first < second
    assert second-first > arb(".000028")
    assert energies[11431]-first > arb(".00624")
    beta = arb(10)
    parity_error = 36*(-beta).exp()*(beta*arb(".025")).sinh()/beta
    assert parity_error < arb(".000042") < arb(".0000928")
    beta = arb(32768)
    assert (arb(2).log()+18*(-beta).exp())/beta < second-first
    print("vacuum winner 11618, second 11619, separation",(second-first).str(10))
    print("vacuum energy excess at maximum quasiparticle gap 11431",(energies[11431]-first).str(10))


if __name__ == "__main__":
    verify()
