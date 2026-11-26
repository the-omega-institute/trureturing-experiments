"""Density of the two existing SU(3) color-clock matrices at zeta zero number eight."""
from pathlib import Path
import runpy

from flint import acb, acb_mat, arb

root = Path(__file__).resolve().parent
sampler = runpy.run_path(str(root/'ch39-shell-standard-model-haar.py'))
gamma = acb.zeta_zero(8).imag
phases = [gamma*sampler['theta'](d) for d in sampler['color']]
generators = [sampler['su3'](phases[i:i+8]) for i in (0,8)]

def identity(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])

def adjoint(g):
    return g.conjugate().transpose()

def hilbert_inner(a,b):
    return sum(a[i,j].conjugate()*b[i,j] for i in range(a.nrows()) for j in range(a.ncols()))

def traceless_basis():
    basis = []
    for i in range(3):
        for j in range(3):
            if i != j:
                e = acb_mat(3,3)
                e[i,j] = 1
                basis.append(e)
    for diag in ((1,-1,0),(1,1,-2)):
        e = acb_mat(3,3)
        length = arb(sum(d*d for d in diag)).sqrt()
        for j,d in enumerate(diag):
            e[j,j] = d/length
        basis.append(e)
    return basis

def commutation_map(g):
    # Columns are the coefficients of GE_ij-E_ijG in End(C^8).
    n = g.nrows()
    out = acb_mat(n*n,n*n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                out[k*n+j,i*n+j] += g[k,i]
                out[i*n+k,i*n+j] -= g[j,k]
    return out

def certify_irreducible():
    basis = traceless_basis()
    gram = acb_mat(64,64)
    for generator in generators:
        images = [generator*b*adjoint(generator) for b in basis]
        action = acb_mat([[hilbert_inner(a,v) for v in images] for a in basis])
        linear = commutation_map(action)
        gram += adjoint(linear)*linear
    for i in range(8):
        for j in range(8):
            gram[9*i,9*j] += 1
    det = gram.det()
    print('augmented commutant Gram determinant lower',det.real.lower().str(12),flush=True)
    assert det.real > 0 and det.imag.contains(0)

def certify_small_element():
    g = generators[0]
    for power in range(1,2049):
        word = g**power
        # Powers of a single unitary admit direct binary powering at full precision.
        delta = word-identity(3)
        bound = hilbert_inner(delta,delta).real.sqrt()
        if bound<arb(1)/2:
            off = abs(word[0,1])
            if off>0:
                print('small noncentral power, Frobenius bound and off-diagonal lower',
                      power,bound.upper().str(12),off.lower().str(12),flush=True)
                return
    raise RuntimeError('no small noncentral power certified')

if __name__=='__main__':
    for g in generators:
        assert abs(g.det()-1)<arb('1e-40')
    certify_irreducible()
    certify_small_element()
