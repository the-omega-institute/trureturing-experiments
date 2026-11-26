"""Sharp canonical dimensions for both lepton-number signs in triple-baryon local operators."""
from itertools import product
from pathlib import Path
import runpy

alg=runpy.run_path(str(Path(__file__).with_name("ch39-triple-baryon-operator.py")))
wedge,add,neutron=alg["wedge"],alg["sum_terms"],lambda j:alg["baryon"]("d",j)


def shifted(poly,bits):
    return {k<<bits:v for k,v in poly.items()}


def contract(left,right):
    return add([((-1)**j,wedge(left[j],right[1-j])) for j in range(2)])


if __name__=="__main__":
    n=[neutron(j) for j in range(2)]
    n2=[shifted(v,12) for v in n]
    # A single nonzero weak monomial L^1 H^2 suffices: the other component multiplies H^1.
    lh=[{1<<(24+j):1} for j in range(2)]
    poly=wedge(contract(n,n),contract(n2,lh))
    assert poly and all(k.bit_count()==10 for k in poly)
    assert 3*(-4)+6*2-3+3==0
    # q,u,l,e are signed multiplicities of Q, u^c+d^c, L, e^c+nu^c.
    # Colour divisibility needs residues modulo twelve to retain B modulo four.
    for q,u,l,e,h in product(range(12),range(12),range(4),range(4),range(2)):
        if (q-u)%3 or (q+l+h)%2 or (q+u+l+e)%2:
            continue
        B=(q-u)//3
        lep=l-e
        assert (-(q+u+l+e)+2*h-B+lep)%4==0
    print("dimension-sixteen component monomials, largest coefficient",len(poly),max(map(abs,poly.values())))
