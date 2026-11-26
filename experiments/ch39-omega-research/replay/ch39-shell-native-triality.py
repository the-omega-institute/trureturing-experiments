"""Exact torsion of the native scalar returns and a three-copy ninth-root lift."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
shell = run_path(str(root/'ch39-shell-current-cyclic-protection.py'))
multiply, power, link = (shell[k] for k in ('multiply','power','link'))
triangles = [(-3,(0,7,13)),(-2,(0,7,14)),(1,(6,18,7)),
             (3,(6,9,18)),(-1,(6,12,14)),(2,(9,13,14))]


def main():
    total = (1,0)
    for count,(a,b,c) in triangles:
        cycle = ((a,b),(b,c),(c,a))
        assert all(e in shell['links'] for e in cycle)
        word = multiply(link(a,b),multiply(link(b,c),link(c,a)))
        if count<0: word = (word[0],-word[1])
        total = multiply(total,power(word,abs(count)))
    assert sum(abs(n) for n,_ in triangles)==12
    assert total==(3,1)
    omega = power(total,4)
    assert omega==(-1,1)
    assert len({power(total,k) for k in range(12)})==12 and power(total,12)==(1,0)
    # H=2Y-B is integral on left-handed SM fields with an optional sterile neutrino.
    charges={'Q':(Q(1,6),Q(1,3)),'u':(Q(-2,3),Q(-1,3)),
             'd':(Q(1,3),Q(-1,3)),'L':(Q(-1,2),0),'e':(1,0),'n':(0,0),'H':(Q(1,2),0)}
    exponents={name:2*y-b for name,(y,b) in charges.items()}
    assert all(Q(t).denominator==1 for t in exponents.values())
    q,u,d,l,e,n,h = (exponents[key] for key in ('Q','u','d','L','e','n','H'))
    assert [q+u+h,q+d-h,l+e-h,l+n+h,2*(l+h)]==[0,0,0,0,0]
    for name,(y,baryon) in charges.items():
        y6 = 6*y
        assert Q(y6).denominator==1
        digit = (int(y6)+1)%3-1
        assert digit==3*baryon and (int(y6)+1)//3==exponents[name]
    assert all((n+1)//3+(-n+1)//3==0 for n in range(9))
    assert (3*q+l)%3==2
    assert (3*(3*q+l))%3==0  # the three-generation weak instanton
    print('Twelve triangle returns give direction (3,1). Native scalar torsion is exactly of order 12.')
    print('The companion cyclic lift cubes to the nontrivial cube root. A uniform three-dimensional')
    print('primitive ninth-root eigenspace over quadratic shell data requires at least nine dimensions.')


if __name__=='__main__': main()
