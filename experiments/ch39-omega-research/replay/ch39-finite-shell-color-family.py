"""Color and weak multiplicities in permutation records and an extra output register."""
from itertools import combinations, permutations

from sympy import I, Matrix, eye


def main():
    shell = list(permutations(range(4)))
    H = Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
    C = [H.T*H.extract(p,range(3)) for p in shell]
    X = Matrix.vstack(*[c.reshape(1,9) for c in C])
    assert X.T*X == 8*eye(9)
    assert X*X.T == Matrix([[sum(x==y for x,y in zip(p,q))-1 for q in shell] for p in shell])
    for p,c in zip(shell,C):
        first, second = [{i for i,v in enumerate(p) if v in pair} for pair in ((0,1),(0,2))]
        third = (first & second) | (set(range(4))-(first|second))
        for j,cut in enumerate((first,second,third)):
            signal = Matrix([1 if i in cut else -1 for i in range(4)])/2
            assert H.T*signal == c[:,j]
    # Color is the standard 3, weak is the nonconstant part of the
    # permutation representation on the three unordered two-pair partitions.
    cuts = [frozenset((frozenset((0,j)),frozenset(set(range(4))-{0,j}))) for j in (1,2,3)]
    partition_action = []
    for p in shell:
        image = [frozenset(frozenset(p[i] for i in side) for side in cut) for cut in cuts]
        partition_action.append(Matrix([[int(image[i]==cut) for cut in cuts] for i in range(3)]))
    color_char = [g.trace() for g in C]
    weak_char = [g.trace()-1 for g in partition_action]
    assert sum(x*x for x in color_char)==24 and sum(x*x for x in weak_char)==24
    assert sum(x*y for x,y in zip(color_char,weak_char))==0
    regular = [24]+[0]*23
    augmented = [x+g.trace() for x,g in zip(regular,partition_action)]
    multiplicity = lambda rep,chi: sum(x*y for x,y in zip(rep,chi))/24
    assert [multiplicity(regular,chi) for chi in (color_char,weak_char)]==[3,2]
    assert [multiplicity(augmented,chi) for chi in (color_char,weak_char)]==[3,3]
    for p,c,g in zip(shell,C,partition_action):
        if p[0]==0:
            assert c==g
    assert sum((a-(b+1))**2 for a,b in zip(color_char,weak_char))>0
    # A scalar convolution can act on the right without mixing colors.
    # Independent coefficients for each allowed label operation test local capacity.
    vals = (0,1,3,4)
    cost = [sum((vals[i]-vals[j])**2 for i,j in enumerate(p))//2 for p in shell]
    word = [sum(p[i]>p[j] for i,j in combinations(range(4),2)) for p in shell]
    for lengths, cutoffs, ranks in [
        (cost,(0,1,2,4,7,9,10,11),((1,1),(3,3),(3,3),(4,4),(8,6),(8,6),(8,6),(9,6))),
        (word,(0,1,2,3),((1,1),(4,4),(8,6),(9,6)))]:
        for bound,expected in zip(cutoffs,ranks):
            mats=[c for c,k in zip(C,lengths) if k<=bound]
            full=Matrix.hstack(*[c.reshape(9,1) for c in mats]).rank()
            symmetric=Matrix.hstack(*[(c+c.T).reshape(9,1) for c in mats]).rank()
            assert (full,symmetric)==expected
    reversal = C[shell.index((3,2,1,0))]
    for d in set(cost):
        radial = sum((c for c,k in zip(C,cost) if k==d),0*eye(3))
        assert radial*reversal==reversal*radial
    assert reversal.diagonal().tolist()[0]==[-1,-1,1]
    a,b,c = (C[shell.index(p)] for p in ((1,0,2,3),(0,2,1,3),(0,1,3,2)))
    up = 3*eye(3)+a+2*c
    down = eye(3)+b+I*a+c
    hu,hd = up*up.H,down*down.H
    assert (hu*hd-hd*hu).det()==-30720*I
    print("Nine orthogonal joint records. Color/weak multiplicities 3/2 -> 3/3 after the three-output register.")
    print("Fixing one label identifies family actions. Symmetric flavor capacity starts at chord cost 7, full capacity at 11.")


if __name__=="__main__":
    main()
