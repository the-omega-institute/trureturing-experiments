import sys,random; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from wormald_hessian_core import *
import sympy as sp
random.seed(11)
for d in (3,5,7,9,11):
  A=circ(d,set(range(1,(d+1)//2)))
  for rep in range(3):
    for _ in range(40):
      i,j,k=random.sample(range(d),3)
      if A[i][j] and A[j][k] and A[k][i]: A[i][j]=A[j][k]=A[k][i]=0; A[j][i]=A[k][j]=A[i][k]=1
    K0,H=calH(A); r=Fr(d)/H
    S=sp.Matrix(d,d,lambda i,j:A[i][j]-A[j][i]); I=sp.eye(d)
    num=(d*I+S).det(); den=sp.expand((S+sp.I*d*I).det()/sp.I)
    print(d, sp.Rational(r.numerator,r.denominator)==sp.Rational(num,den)**2, sp.Rational(num,den))
