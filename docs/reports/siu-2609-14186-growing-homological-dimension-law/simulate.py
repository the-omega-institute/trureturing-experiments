#!/usr/bin/env python3
"""Exploratory floating-point Wishart simulation; not exact or proof evidence."""
import json, math, random
import numpy as np

SEED=260914186
SAMPLES=8000
rng=np.random.default_rng(SEED)
prng=random.Random(SEED)

def mp_integral(k,gamma):
    lo=(1-gamma**(-0.5))**2
    hi=(1+gamma**(-0.5))**2
    z,v=np.polynomial.legendre.leggauss(256)
    angle=(z+1)*np.pi/2
    t=(lo+hi)/2+(hi-lo)/2*np.cos(angle)
    jac=(hi-lo)/2*np.sin(angle)
    density=gamma/(2*np.pi*t)*np.sqrt((hi-t)*(t-lo))
    return float(np.dot(v,t**(-2/(k-2))*density*jac)*np.pi/2)

def moments(a):
    a=np.asarray(a)
    return {'mean':np.mean(a,axis=0).tolist(),'std':np.std(a,axis=0).tolist(),'samples':len(a)}

def point(k,D,w,J,i,kind):
    r=(k-2)/2
    S=float(np.sum(w[J]))
    wi=float(w[i])
    if kind=='removed':
        L=S**(-r)/k
        ell=L*np.expm1(-r*np.log1p(-wi/S))
    else:
        v=S**(-r)/k
        small=(S+wi)**(-r)/k
        ell=-v*np.expm1(-r*np.log1p(wi/S))
        L=-v if kind=='negative' else small
    return (D**r*L,D**(r+1)*ell)

print('EXPLORATORY_FLOATING_POINT_ONLY')
print(json.dumps({'seed':SEED,'samples_per_family':SAMPLES,'gamma':2,'k':[3,4],'D':[200,400]}))
for D in (200,400):
    X=rng.standard_normal((D,2*D))
    lam=np.linalg.eigvalsh(X@X.T/(2*D))[::-1]
    print(json.dumps({'D':D,'lambda_min':float(lam[-1]),'lambda_max':float(lam[0]),
                      'edge_targets':[(1-2**(-0.5))**2,(1+2**(-0.5))**2]}))
    for k in (3,4):
        r=(k-2)/2
        w=lam**(-1/r)
        I=mp_integral(k,2)
        wp=(1+2**(-0.5))**(-2/r)
        def target(s): return [1/(k*(s*I)**r),r/k*wp*(s*I)**(-r-1)]
        print(json.dumps({'D':D,'k':k,'empirical_I':float(np.mean(w)),'MP_I':I}))
        for theta in (1/3,1/2):
            q=int(theta*D)
            if k%2==0:
                for sup in (False,True):
                    m=q+2 if sup else D-q
                    A=[2**(m-i)*math.comb(D-i,m-i) for i in range(1,m+1)]
                    cumulative=[];total=0
                    for a in A:total+=a;cumulative.append(total)
                    import bisect
                    out=[]
                    for _ in range(SAMPLES):
                        i=bisect.bisect_right(cumulative,prng.randrange(total))+1
                        J=np.concatenate((np.arange(i),rng.choice(np.arange(i,D),m-i,replace=False)))
                        out.append(point(k,D,w,J,i-1,'removed'))
                    t=target(theta if sup else 1-theta)
                    print(json.dumps({'D':D,'k':k,'theta':theta,'q':q,'filtration':'super' if sup else 'sub',
                        'target':t,**moments(out)}))
            else:
                neg=[];pos=[]
                for kind,m,collector in [('negative',q+1,neg),('positive',D-q-1,pos)]:
                    for _ in range(SAMPLES):
                        J=rng.choice(np.arange(1,D),m,replace=False)
                        collector.append(point(k,D,w,J,0,kind))
                tneg=target(theta);tneg[0]*=-1
                tpos=target(1-theta)
                print(json.dumps({'D':D,'k':k,'theta':theta,'q':q,'filtration':'sub',
                    'negative_mass':(D-q-1)/D,'negative_target':tneg,'positive_target':tpos,
                    'negative':moments(neg),'positive':moments(pos)}))
                # Exact geometric transformation applied to numerical bar readings.
                for a in (neg,pos):
                    for j,(x,y) in enumerate(a):a[j]=(-x-y/D,y)
                tneg[0]*=-1;tpos[0]*=-1
                print(json.dumps({'D':D,'k':k,'theta':theta,'q':q,'filtration':'super',
                    'negative_family_target':tneg,'positive_family_target':tpos,
                    'negative_family':moments(neg),'positive_family':moments(pos)}))
print('SIMULATION_FINISHED')
