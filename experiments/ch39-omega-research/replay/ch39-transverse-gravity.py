"""Exact transverse metric identities used by the lapse-constraint test."""
from flint import fmpq_mpoly_ctx

ring = fmpq_mpoly_ctx.get(('E','q','P','Q','omega'))
E,q,P,Q,omega = ring.gens()
def mul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def trace(a):
    return a[0][0]+a[1][1]
def det(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]

# H, H^-1 and dH are represented after multiplication by E=exp(p).
h = [[E**2, E**2*q], [E**2*q,1+E**2*q**2]]
inverse = [[h[1][1],-h[0][1]],[-h[1][0],h[0][0]]]
deriv = [[E**2*P,E**2*(P*q+Q)],[E**2*(P*q+Q),-P+E**2*(P*q*q+2*q*Q)]]
assert det(h)==E**2
current = mul(inverse,deriv)
assert trace(current)==0
assert trace(mul(current,current))==2*E**4*(P**2+E**2*Q**2)

# At a principal axis of a rotating unimodular block the derivative is off-diagonal.
rotated = [[0,omega*(E**2-1)],[omega*(E**2-1),0]]
current = mul([[1,0],[0,E**2]],rotated)
assert trace(mul(current,current))==2*omega**2*E**2*(E**2-1)**2
print('unimodular shear norm and rotating-axis scalar-curvature identities verified over Q')
