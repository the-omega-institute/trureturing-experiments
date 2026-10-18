"""The kernel identity behind the proof: for h = d, (M_d)_{ij} = 2 sum_{r,l in Z} K(i + d r, j + d l), where
K(x, y) = sqrt(2/(d a)) exp(-pi((x+y)^2 + D (x-y)^2)/(2 d a) + i pi c (x^2 - y^2)/(d a)) and D = a b - c^2.
Prints the largest entrywise error for d = 3, 5, 7 and four covariances. Requires numpy and check.py.
"""
import math,cmath,numpy as np
src=open('check.py').read().split('rng=np.random')[0]; ns={}; exec(src,ns); M=ns['M']
def K(x,y,a,b,c,d):
    D=a*b-c*c
    return math.sqrt(2/(d*a))*cmath.exp(-math.pi*((x+y)**2+D*(x-y)**2)/(2*d*a)+1j*math.pi*c*(x*x-y*y)/(d*a))
worst=0
for d in (3,5,7):
    s=(d-1)//2; I=list(range(-s,s+1))
    for (a,b,c) in [(1.0,1.0,0.0),(0.7,2.5,0.9),(2.0,1.3,-1.1),(0.4,4.0,0.6)]:
        if a*b-c*c<1: b=(1.2+c*c)/a
        Md=M(d,a,b,c,float(d))
        R=range(-30,31)
        for i in I:
            for j in I:
                v=2*sum(K(i+d*r,j+d*l,a,b,c,d) for r in R for l in R)
                worst=max(worst,abs(v-Md[I.index(i),I.index(j)]))
print('max |identity error|',worst)
