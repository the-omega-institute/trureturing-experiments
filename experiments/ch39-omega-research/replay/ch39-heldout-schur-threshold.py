"""The ninth time node compared to the central 266.6444299 pivot threshold."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, arb_poly, ctx

root = Path(__file__).resolve().parent
old = run_path(str(root/'ch39-eight-time-pole-control.py'))
ctx.prec = 384

def real(x):
    return arb(x.numerator)/x.denominator

def enclose(cell):
    return real(cell.lo).union(real(cell.hi))

def time(j):
    if j==0: return arb(1)
    k=old['levels'][j]
    return arb(k).union(arb(k+1))/2**288

times = [time(i) for i in range(8)]
gram = arb_mat([[(i*times[i]+j*times[j])/(i+j) if i+j else arb(1)
                 for j in range(8)] for i in range(8)])
for n in range(1,9):
    assert arb_mat([[gram[i,j] for i in range(n)] for j in range(n)]).det()>0
inverse = gram.inv()
tables = [[list(map(enclose,part)) for part in old['quadrature'](ids,Q(0))]
          for ids in ([1,2,3,4,8,9,51,52],[0,1,2,3,5,6,15,16])]
criterion = 1/arb('266.6444299')**2

def profile(j,table):
    factors = [1+j*j*v for v in table[0]]
    return sum(w*prod(factors[:i]+factors[i+1:]) for i,w in enumerate(table[1])),prod(factors)

def form(j):
    terms = [j+i for i in range(8)]
    cofactors = [prod(terms[:i]+terms[i+1:]) for i in range(8)]
    a = [u*(i*times[i]) for i,u in enumerate(cofactors)]
    b = [u*j for u in cofactors]
    def dot(v,w):
        return sum(v[i]*w[h]*inverse[i,h] for i in range(8) for h in range(8))
    return prod(terms),dot(a,a),dot(b,b),dot(a,b)

def main():
    for j in range(8,29):
        g,a,b,c = form(arb(j))
        a,b,c = [x/g**2 for x in (a,b,c)]
        l,h = [u/v for u,v in (profile(arb(j),table) for table in tables)]
        def pivot(t): return 1-2*c-a/t-b*t
        assert a>0 and b>0 and l<h<(a/b).sqrt()
        if j<28:
            assert 0<pivot(l)<pivot(h)<criterion
            continue
        assert 0<pivot(l)<criterion<pivot(h)
        for x,low,high in ((l,'278.12384','278.12385'),(h,'262.12993','262.12994'),
                           (time(j),'269.76033','269.76034')):
            assert arb(low)<1/pivot(x).sqrt()<arb(high)
        d=1-2*c-criterion
        critical=2*a/(d+(d*d-4*a*b).sqrt())
        assert arb('0.62580730116807')<critical<arb('0.62580730116808')
        gap=critical-time(j)
        assert 1/arb(2)**19<gap<1/arb(2)**18
        print('J=28 critical time',critical.str(20),'independent RH inverse pivot',(1/pivot(time(j)).sqrt()).str(16))
    # Clear positive denominators and certify the entire half-line at once.
    j=arb_poly([29,1])
    g,a,b,c=form(j)
    for table in tables:
        n,w=profile(j,table)
        certificate=(g*g*(1-criterion)-2*c)*n*w-a*w*w-b*n*n
        assert all(v>0 for v in certificate.coeffs())
    print('Inverse pivot above the comparison for J=8..27, both sides at 28, below for every J>=29.')

if __name__=='__main__': main()
