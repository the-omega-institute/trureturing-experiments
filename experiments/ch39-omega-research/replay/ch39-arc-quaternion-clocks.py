"""Two clocked great-circle normals and the tracial tests at zero 35."""
from fractions import Fraction
from flint import acb,arb,ctx
ctx.prec=192
def q(f):
    v=Fraction(f)
    return arb(v.numerator)/v.denominator
def dot(x,y):
    return sum(a*b for a,b in zip(x,y))
def cross(x,y):
    return tuple(x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3))
P=(1,5,29)
v=cross(P,(13,23,13))
w=cross(P,(23,13,13))
angle_sq=Fraction(dot(v,w)**2,dot(v,v)*dot(w,w))
assert angle_sq==Fraction(8741547,12563572) and dot(v,w)>0
c=q(angle_sq).sqrt()
def theta(d):
    return (1-q(d)/1734).acos()
h=acb.zeta_zero(35).imag
p=[h*(theta(d)-theta(244)) for d in (724,804)]
def report(width):
    (s,x),(t,y)=[(z+arb(0,width.upper())).sin_cos() for z in p]
    change=(1-c)*s*t
    gram=(1-c*c)*(s*t)**2
    comm=1-2*gram
    five=(1-x)**2*gram**2/8
    for label,value in [('axis cosine',c),('mixed trace decrease',change),
                        ('three-anchor Gram determinant',gram),('commutator half-trace',comm),
                        ('five-anchor kappa Gram determinant',five)]:
        print(label,value.lower().str(12),value.upper().str(12))
    assert change>0 and gram>0 and five>0
report(q(0))
report(q('.01'))

short_a=(11,11,25)
short_b=(13,13,23)
assert Fraction(dot(short_a,short_b),867)==Fraction(287,289)
h2=Fraction(1,17**2)
s,x=p[0].sin_cos()
five=(1-x)**2*s**4*q(h2)**2/8
print('open shortest rotor, five-anchor determinant',five.str(15))
