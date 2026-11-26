"""Exact replay of the sphere comparison table using finite integer records only."""
from fractions import Fraction as Q
from itertools import product

records = [
    (4, -1, -2031855, -1, -2665916),
    (5, -1, 1867270, 1, -3971534),
    (8, 1, -3475440, -1, 1017009),
    (12, 1, 3628694, 1, -1821260),
    (15, 1, 2718931, -1, -1586390),
]
denominator = 2**22


def rotor(sign, numerator):
    t = Q(numerator,denominator)
    return sign*(1-t*t)/(1+t*t), sign*2*t/(1+t*t)


def mul(z,w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def dot(z,w):
    return z[0]*w[0]+z[1]*w[1]


def solve(A,b):
    a = [list(row)+[v] for row,v in zip(A,b)]
    for j in range(len(a)):
        k = next(i for i in range(j,len(a)) if a[i][j])
        a[k],a[j] = a[j],a[k]
        a[j] = [v/a[j][j] for v in a[j]]
        for i in range(len(a)):
            if i != j:
                a[i] = [v-a[i][j]*u for v,u in zip(a[i],a[j])]
    return [row[-1] for row in a]


def basis_and_weights():
    columns = []
    for _,su,tu,sv,tv in records:
        u,v = rotor(su,tu),rotor(sv,tv)
        assert dot(u,u)==dot(v,v)==1
        columns.append([(Q(1),Q(0)),u,v,mul(u,v)])
    equations = [[Q(1)]*5]+[[dot(c[i],c[j]) for c in columns]
                              for i,j in [(0,1),(0,2),(1,2),(0,3)]]
    weights = solve(equations,[Q(1)]+[Q(0)]*4)
    assert sum(weights)==1 and min(weights)>Q(1,10)
    assert all(sum(w*dot(c[i],c[j]) for w,c in zip(weights,columns))==int(i==j)
               for i,j in product(range(4),repeat=2))
    return columns,weights


def check(radius, digit, bound):
    cols,weights = basis_and_weights()
    # Include both signs of each frequency to obtain a complex orthonormal frame.
    atoms = [(w/2,[(x,s*y) for x,y in c]) for w,c in zip(weights,cols) for s in (1,-1)]
    points = [p for p in product(range(bound+1),repeat=3) if sum(x*x for x in p)==radius]
    assert digit**2+1 == radius
    def amplitudes(p):
        coeff = [(p[0],0),(p[1],0),(p[2],0),(digit,1)]
        return [tuple(map(sum,zip(*(mul(v,z) for v,z in zip(coeff,c))))) for _,c in atoms]
    vectors = list(map(amplitudes,points))
    for (p,v),(r,u) in product(zip(points,vectors),repeat=2):
        real = sum(w*dot(x,y) for (w,_),x,y in zip(atoms,v,u))
        imag = sum(w*(x[0]*y[1]-x[1]*y[0]) for (w,_),x,y in zip(atoms,v,u))
        assert imag==0 and real==radius+sum(i*j for i,j in zip(p,r))
    readings = [Q(radius+sum(i*j for i,j in zip(p,r)),2*radius)
                for p in points for r in points if p!=r]
    assert max(readings)==1-Q(1,2*radius)
    return points,readings


if __name__ == "__main__":
    p,readings = check(17,4,4)
    assert len(p)==9 and Q(14,17) in readings
    p,_ = check(26,5,4)
    assert len(p)==6
    print("Pure rational replay verified 81 + 36 ordered entries, including shortest 33/34 and 51/52")
