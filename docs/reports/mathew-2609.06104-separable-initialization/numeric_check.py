#!/usr/bin/env python3
"""Exploratory floating-point random-product search, separate from exact evidence.

For each initial product state maximize over every time analytically: the power
is a*cos(2 sqrt(m)t)+b*sin(2 sqrt(m)t). Then use alternating largest-eigenvector
updates from random product seeds. These readings are not a global certificate.
"""
import numpy as np
from check import model


def dense(op,n):
    size=2**n
    out=np.zeros((size,size),dtype=complex)
    columns=np.arange(size)
    for (x,z),coefficient in op.items():
        signs=np.array([(-1)**((int(j)&z).bit_count()%2) for j in columns])
        out[columns^x,columns] += complex(coefficient)*signs
    return out


def random_vector(rng,d):
    v=rng.normal(size=d)+1j*rng.normal(size=d)
    return v/np.linalg.norm(v)


def expect(op,v,w):
    return float(np.einsum('i,a,iajb,j,b',v.conj(),w.conj(),op,v,w,optimize=True).real)


def check(m):
    rng=np.random.default_rng(260906104+m)
    nb,nq,*rest=model(m)
    hb,c=rest[-3:-1]
    hc=dense(c,nb+nq)
    # check.py numbers battery qubits first (least significant): kron is aux,battery.
    joint_hb=np.kron(np.eye(2**nq),dense(hb,nb))
    p=1j*(hc@joint_hb-joint_hb@hc)
    q=1j*hc@p/np.sqrt(m)
    db,da=2**nb,2**nq
    p=p.reshape(da,db,da,db).transpose(1,0,3,2)
    q=q.reshape(da,db,da,db).transpose(1,0,3,2)
    sample_best=0.0
    for _ in range(1024):
        v,w=random_vector(rng,db),random_vector(rng,da)
        a,b=expect(p,v,w),expect(q,v,w)
        sample_best=max(sample_best,float(np.hypot(a,b)))
    optimized_best=0.0
    for _ in range(32):
        v,w=random_vector(rng,db),random_vector(rng,da)
        for _ in range(20):
            a,b=expect(p,v,w),expect(q,v,w)
            angle=np.arctan2(b,a)
            op=np.cos(angle)*p+np.sin(angle)*q
            battery=np.einsum('a,iajb,b->ij',w.conj(),op,w,optimize=True)
            v=np.linalg.eigh(battery)[1][:,-1]
            auxiliary=np.einsum('i,iajb,j->ab',v.conj(),op,v,optimize=True)
            w=np.linalg.eigh(auxiliary)[1][:,-1]
        optimized_best=max(optimized_best,float(np.hypot(expect(p,v,w),expect(q,v,w))))
    bound=2*np.sqrt(2)*np.sqrt(m+2*(m//2))
    benchmark=4*np.sqrt(m)
    assert optimized_best <= bound+1e-9
    print(f'm={m}; seed={260906104+m}; random_products=1024; time_max_sample={sample_best:.12f}; '
          f'alternating_starts=32; sweeps=20; time_max_optimized={optimized_best:.12f}; '
          f'separable_bound={bound:.12f}; benchmark={benchmark:.12f}; '
          f'optimized_minus_bound={optimized_best-bound:.3e}')


if __name__=='__main__':
    for m in (1,2):check(m)
    print('EXPLORATORY_NUMERICS_OK')
