"""Flavor operators on the explicit proton-protecting Pati-Salam Higgs orbit."""
from itertools import product
from flint import fmpz_mat

I=fmpz_mat([[1,0,0],[0,1,0],[0,0,1]])
A=fmpz_mat([[-1,0,0],[0,0,0],[0,0,1]])
B=fmpz_mat([[0,1,1],[1,0,2],[1,2,0]])
v=fmpz_mat([[1],[1],[1]])
def flatten(x):
    return [x[i,j] for i,j in product(range(3),repeat=2)]
def matrix_from_cols(columns):
    return fmpz_mat([[x[i,0] for x in columns] for i in range(3)])
basis=[I]
ranks=[1]
for depth in range(1,4):
    words=[I]
    for _ in range(depth):
        words=[a*b for a in words for b in (A,B)]
    for w in words:
        if fmpz_mat([flatten(x) for x in basis+[w]]).rank()>len(basis):
            basis.append(w)
    ranks.append(len(basis))
assert ranks==[1,3,7,9]
N=matrix_from_cols([v,A*v,B*v])
assert abs(N.det())==1
directions=[I*v,A*v,B*v]
symmetric=[x*y.transpose()+y*x.transpose()
           for i,x in enumerate(directions) for y in directions[i:]]
assert fmpz_mat([flatten(x) for x in symmetric]).rank()==6
assert abs(sum((x*x.transpose() for x in directions),fmpz_mat(3,3)).det())==1

# Two independent gauge contractions select quarks and leptons. With a bidoublet
# vacuum in only one neutral component, H and its conjugate select the two R channels.
coefficient=fmpz_mat([[1,0,0,0],[0,1,0,0],[1,0,1,0],[0,1,0,1]])
assert abs(coefficient.det())==1
print("family matrix algebra ranks at operator-word lengths 0,1,2,3:",ranks)
print("determinant of (v,Av,Bv), dimension of symmetric mass space",N.det(),len(symmetric))
