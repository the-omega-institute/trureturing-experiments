"""Common-parent Theorem 2.1, including its singular source lines.
The samples and responses are exactly XI 6.1.7--6.1.8; this is not a new lift.
"""
import numpy as np

def world(r,b,d):
    j=b+r*(1-r)/b;n=j+r;T=1+b;U=1+2*b
    q0=-d;q1=(24-d)/j;q2=(147/5-d)/T;q3=(24-d)/n
    g4=-q0+4*q1-6*q2+4*q3
    g5=-4*q0+15*q1-20*q2+10*q3
    g6=-10*q0+36*q1-45*q2+20*q3
    g7=-20*q0+70*q1-84*q2+35*q3
    return np.array([d+(U+n)*g5,d+(2*U+n)*g6,d+(T+2*j+4*r)*g7,d+(3*U+2*n)*g7,d+U*g4-18])

def secant(z,term):
    R,B,k,D,l,h=z
    wm=world(R-h,B-k*h,D-l*h);wp=world(R+h,B+k*h,D+l*h)
    diff=(wp-wm)/(2*h)
    return np.r_[(wp[4]+wm[4])/2,diff[[4,0,1,term]]]
def pair_worlds(z):
    R,B,k,D,l,h=z
    return [[R-h,B-k*h,D-l*h],[R+h,B+k*h,D+l*h]]
def valid(z):
    w=pair_worlds(z)
    return z[5]>1e-7 and all(0<r<1 and b>0 and np.all(np.isfinite([r,b,d])) for r,b,d in w)
