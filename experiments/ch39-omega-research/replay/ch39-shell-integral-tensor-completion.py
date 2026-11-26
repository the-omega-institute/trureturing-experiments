"""Integral mixed Green-Schwarz factorization and a full-rank mass on twelve charged zero modes."""
from fractions import Fraction as F
from flint import fmpz_mat

# Ordered keys are x^2, c2(A), c2(B). The determinant-neutral singlet restores rank.
mixed=(F(5,2),-3,-4)
det=(F(1,2),0,0)
Y=tuple(a-b for a,b in zip(mixed,det))
assert Y==(2,-3,-4)
rank=12-1+1
assert rank==12 and rank//12==1 and F(rank,24)==F(1,2)
print("integral Y coefficients for x^2,c2(A),c2(B),p1/2,c2(L)+c2(R)",*Y,-1,-1)
# Bianchi-relevant internal-flux component after splitting the bundle as
# A=L_-3 tensor A4, B=L_4 tensor B3. The two-dimensional base has t^2=0.
# Writing x=x4+12t gives: c2(A)=c2(A4)+9*x4*t; c2(B)=c2(B3)+8*x4*t.
# Signs follow c1(A4)=-x4, c1(B3)=x4.
axion=2*24-3*9-4*8
assert axion==-11
# 4D mixed anomaly for left weak: original six, minus 72 from the spectator.
assert (12-144)//2==6*axion
print("fiber integral of Y divided by x4",axion,"mixed weak anomalies",6*axion,-6*axion)
W=fmpz_mat([[j-i if i+j==11 else 0 for j in range(12)] for i in range(12)])
assert W.transpose()==-W and W.rank()==12
assert abs(int(W.det()))==(1*3*5*7*9*11)**2
print("Wronskian mass rank",W.rank(),"absolute Pfaffian",1*3*5*7*9*11)
