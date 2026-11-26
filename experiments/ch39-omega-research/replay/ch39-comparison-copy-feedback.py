"""Visible transfer at a signed port of the ordered shell comparisons."""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from sympy import I, Matrix, Rational as Q, cancel, eye, sqrt, symbols, zeros


def copy_channel(size, x):
    pairs = [(i,j) for i,j in product(range(size), repeat=2) if i != j]
    incidence = zeros(size,len(pairs))
    for j, (a,b) in enumerate(pairs):
        incidence[a,j] = incidence[b,j] = 1
    gap2 = 2*(size-2)
    assert incidence*incidence.T == gap2*eye(size)+2*Matrix.ones(size)
    labels = range(size//2)
    # Each column is a real port component. The three receive synchronized clocks.
    port = zeros(len(pairs),3)
    for i, j in combinations(labels, 2):
        for offset, sign in ((0,1), (size//2,-1)):
            forward = pairs.index((i+offset,j+offset))
            backward = pairs.index((j+offset,i+offset))
            port[forward,:] = sign*Matrix([[3,0,0]])
            port[backward,:] = sign*Matrix([[2,1,2]])
    moved = incidence*port
    assert Matrix.ones(1,size)*moved == zeros(1,3)
    norm2 = port.norm()**2
    assert norm2 == 9*size*(size//2-1)
    assert moved.norm()**2/(gap2*norm2) == Q(5,6)
    adjacency = zeros(size+len(pairs))
    adjacency[:size,size:] = incidence
    adjacency[size:,:size] = incidence.T
    lifted = zeros(size,3).col_join(port)
    assert adjacency**3*lifted == gap2*adjacency*lifted
    assert (lifted.T*adjacency*lifted).trace() == 0
    assert (lifted.T*adjacency**2*lifted).trace()/norm2 == Q(5,6)*gap2
    return cancel(1/(6*x)+5*x/(6*(x**2-gap2)))


def verify():
    root = Path(__file__).parent
    run_path(str(root/'ch39-intrinsic-shell-comparisons.py'))['verify']()
    x = symbols('x', nonzero=True)
    four, six = [copy_channel(size,x) for size in (4,6)]
    assert cancel(six-four.subs(x,x/sqrt(2))/sqrt(2)) == 0
    aligned = cancel(six.subs(x,x/sqrt(Q(15,2)))/sqrt(Q(15,2)))
    assert cancel(aligned-four.subs(x,x/sqrt(15))/sqrt(15)) == 0
    jacobi = Matrix([[0,sqrt(50),0], [sqrt(50),0,sqrt(10)], [0,sqrt(10),0]])
    assert cancel(aligned-(x*eye(3)-jacobi).inv()[0,0]) == 0
    admittance = 12*aligned
    pade = (120+60*I*x-12*x*x-I*x**3)/(120-60*I*x-12*x*x+I*x**3)
    assert cancel((admittance+I)/(admittance-I)-pade) == 0
    # For scalar amplitudes the required ratio is a positive root of r^2-3r+1.
    for ratio in ((3+sqrt(5))/2, (3-sqrt(5))/2):
        assert ((ratio+1)**2/(2*(ratio**2+1))-Q(5,6)).simplify() == 0
    from flint import arb, ctx
    arithmetic = run_path(str(root/'ch39-shell-channel-prime-readings.py'))
    ctx.prec, ctx.cap = 384, 3
    rules = arithmetic['rules'](arithmetic['harmonic_measures'](16))
    targets = arithmetic['windows']['arithmetic_times']()
    for n, den, lo, hi in ((9,4,255,256), (4,9,277,278), (1,1,285,286),
                           (3,1,251,252), (1,3,264,265), (5,3,261,262), (3,5,296,298)):
        pivot = arithmetic['readings'](arb(n)/den, rules, targets)[2]
        assert lo < pivot < hi
    print('Squared copy frequencies: 4 and 8. Signed port spectral weights: 1/6, 5/12, 5/12.')
    print('Aligned tridiagonal coupling squares: 50 and 10; aggregate port factor: 12.')


if __name__ == '__main__':
    verify()
