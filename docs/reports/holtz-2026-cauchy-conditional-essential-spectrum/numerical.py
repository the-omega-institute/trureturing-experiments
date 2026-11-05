#!/usr/bin/env python3
"""Floating point diagnostics only; check.py contains the exact evidence.
Nyström matrices discretize the actual Q_eta after its unitary conjugation.
Finite-section eigenvalues of this nonnormal operator do not certify a spectrum.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import sys
import numpy as np
import mpmath as mp

mp.mp.dps = 65
print('Python', sys.version.split()[0], 'NumPy', np.__version__, 'mpmath', mp.__version__)
print('FLOATING POINT DIAGNOSTICS; not proof or exact evidence')


def log_c_scaled(eta, y):
    """log C_eta(eta*y) + pi*eta*abs(y)/2, evaluated before conversion."""
    a = mp.mpf(str(eta))
    x = a * mp.mpf(str(abs(y)))
    return float((a-2)*mp.log(2)-mp.log(mp.pi)-mp.loggamma(a)
                 +2*mp.re(mp.loggamma((a+1j*x)/2))+mp.pi*x/2)


def matrix(eta, nodes, weights):
    d = nodes[:, None]
    y = nodes[None, :]
    same = d*y > 0
    # Stable same-sign evaluation avoids subtracting nearly equal pi*y/2 terms.
    gap = np.pi/2 - np.sign(y)*np.arctan(d)
    gap = np.where(same & (np.abs(d) >= 1), np.arctan(1/np.abs(d)), gap)
    lc = np.array([log_c_scaled(eta, z) for z in nodes])
    logA = (np.log(eta)+lc[None, :]+np.log1p(y*y)/2
            -(eta+1)*np.log1p(d*d)/2-eta*np.abs(y)*gap
            +np.log(weights[:, None]*weights[None, :])/2)
    return np.exp(np.maximum(logA, -745))


def full_section(eta, length, n_tail, n_core=48):
    h = length/n_tail
    s = (np.arange(n_tail)+0.5)*h
    tails = np.exp(s)
    c = -1+(np.arange(n_core)+0.5)*2/n_core
    nodes = np.concatenate((-tails[::-1], c, tails))
    weights = np.concatenate(((tails*h)[::-1], np.full(n_core, 2/n_core), tails*h))
    A = matrix(eta, nodes, weights)
    vals = np.linalg.eigvals(A)
    opnorm = np.linalg.svd(A, compute_uv=False)[0]
    closest = vals[np.argmin(np.abs(vals-1))]
    return len(nodes), opnorm, closest


def tail_diagnostics(eta, shift, length, h=1/12):
    # Input packet on [shift,shift+length], output on the entire positive tail
    # to shift+length+6, rather than compression to the input window.
    outmax = shift+length+6
    s = (np.arange(int(round(outmax/h)))+0.5)*h
    nodes = np.exp(s)
    B = matrix(eta, nodes, nodes*h)
    v = s[None, :]-s[:, None]
    logk = eta*np.log(eta)-float(mp.loggamma(eta))+(eta+0.5)*v-eta*np.exp(v)
    W = np.exp(np.maximum(logk, -745))*h
    cutoff = (s >= shift) & (s <= shift+length)
    window = np.zeros(len(s))
    window[cutoff] = np.sin(np.pi*(s[cutoff]-shift)/length)
    print('TAIL eta=%g shift=%g length=%g h=%g N=%d' % (eta,shift,length,h,len(s)))
    idx = np.where(cutoff)[0]
    diff = B[np.ix_(idx,idx)]-W[np.ix_(idx,idx)]
    print('  Q_tail-W finite-window HS difference = %.8g' % np.linalg.norm(diff))
    for xi in [0, 0.5, 1]:
        p = window*np.exp(1j*xi*s)
        p /= np.linalg.norm(p)
        phi = complex(mp.gamma(eta+mp.mpf('0.5')+1j*xi)/(mp.gamma(eta)*mp.power(eta,mp.mpf('0.5')+1j*xi)))
        residual = np.linalg.norm(B@p-phi*p)
        model_residual = np.linalg.norm(W@p-phi*p)
        ratio = np.linalg.norm(B@p)
        print('  xi=%g phi=(%.8g,%+.8g) residual_Q=%.8g residual_W=%.8g ||Qp||=%.8g' %
              (xi,phi.real,phi.imag,residual,model_residual,ratio))

for eta in [0.5, 1, 3]:
    M = float(mp.gamma(eta+mp.mpf('0.5'))/(mp.sqrt(eta)*mp.gamma(eta)))
    print('ETA %g proposed essential norm/radius %.12g' % (eta,M))
    for length,n in [(6,96),(10,160)]:
        N,norm,val=full_section(eta,length,n)
        print('FULL Q eta=%g tail_length=%g core=48 tail_nodes=%d N=%d norm=%.10g eigenvalue_near_1=(%.10g,%+.3g)' %
              (eta,length,n,N,norm,val.real,val.imag))
    for shift,length in [(2,8),(4,8),(4,16),(4,32)]:
        tail_diagnostics(eta,shift,length)
print('NUMERICAL_DIAGNOSTICS_COMPLETE')
