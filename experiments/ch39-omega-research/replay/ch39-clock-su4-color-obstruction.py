"""Conjugacy obstruction to replacing the central U(1)_A by minimal color SU(4)."""
# Fourth roots of unity are represented by exponents mod 4. Test the simplest
# cocharacter a=diag(1,1,1,-3); the proof in the notes treats arbitrary scaling.
def image(t,color):
    return [((4*t)*k)%4 for k in (1,1,-2)],sorted([(t+k)%4 for k in color]+[(-3*t)%4])

def chiral_singlet(t,color):
    weights=[(y*t-c)%4 for c in color for y in (-4,2)]
    weights += [(6*t)%4,0]
    # Imaginary part of chi minus its conjugate.
    return 2*sum((0,1,0,-1)[a] for a in weights)

left=(0,(1,1,2))
right=(1,(0,1,3))
assert image(*left)==image(*right)
assert sum(left[1])%4==sum(right[1])%4==0
assert chiral_singlet(*left)==-8 and chiral_singlet(*right)==0
print("common SU(3)_X, SU(4) image up to conjugation",image(*left))
print("three-family weak-singlet net character imaginary parts",-24,0)
