"""Sparse exact profiles of the twelve-return generator on the chiral ledger."""
from itertools import combinations
from math import gcd

zero = (0,0,0,0)
one = (1,0,0,0)
generator = (0,1,0,0)


def multiply(a,b):
    coeff = [0]*7
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            coeff[i+j] += x*y
    for i in range(6,3,-1):
        coeff[i-2] += coeff[i]
        coeff[i-4] -= coeff[i]
    return tuple(coeff[:4])


powers = [one]
for _ in range(12):
    powers.append(multiply(powers[-1],generator))
assert powers[-1]==one


def rank(columns):
    table = [list(row) for row in zip(*columns)]
    r = 0
    for j in range(len(columns)):
        pivot = next((i for i in range(r,len(table)) if table[i][j]!=zero),None)
        if pivot is None: continue
        table[r],table[pivot] = table[pivot],table[r]
        for i in range(r+1,len(table)):
            a,b = table[r][j],table[i][j]
            if b==zero: continue
            row = [tuple(u-v for u,v in zip(multiply(a,s),multiply(b,t)))
                   for s,t in zip(table[i],table[r])]
            content = gcd(*(u for entry in row for u in entry))
            table[i] = [tuple(u//content for u in entry) for entry in row] if content else row
        r += 1
    return r


def main():
    ys = [1,-4,2,-3,6,0,3]
    target = [powers[(4*((y+1)//3))%12] for y in ys]
    columns = [[powers[(k*y)%12] for y in ys] for k in range(12)]
    possible = []
    for k in range(1,5):
        for support in combinations(range(12),k):
            selected = [columns[i] for i in support]
            if rank(selected)==rank(selected+[target]): possible.append(support)
    assert possible==[(1,4,7,10)]
    # Fourier coefficients multiplied by four.
    coeff = [tuple(sum(row) for row in zip(*(powers[(v-3*i*j)%12]
              for j,v in enumerate((0,-1,2,1))))) for i in range(4)]
    assert all(c!=zero for c in coeff)
    for j in range(len(ys)):
        value = [sum(row) for row in zip(*(multiply(c,columns[k][j])
                    for c,k in zip(coeff,possible[0])))]
        assert tuple(value)==tuple(4*r for r in target[j])
    for n in (1,2,3,4,6):
        assert any((y-z)%n==0 and target[i]!=target[j]
                   for i,y in enumerate(ys) for j,z in enumerate(ys))
    changes = [{((4*((a+1)//3)-4*((b+1)//3))%12) for a in ys for b in ys
                if (a-b)%12==d} for d in range(12)]
    assert {d for d,row in enumerate(changes) if len(row)==1}=={0,4,8}
    assert 6%12==(-6)%12 and powers[8]!=powers[4]
    print('Minimum finite scalar order 12; unique four-power support 1,4,7,10.')
    print('On both orientations, the torsion-only data identify y=6 with y=-6.')


if __name__=='__main__': main()
