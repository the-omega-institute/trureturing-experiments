"""One-copy orientation measurements on the shortest triangle of the radius-17 shell."""
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    left_right = run_path(str(Path(__file__).with_name('ch39-omega-projector-clock.py')))['multiply_by_vector']
    eye = sp.eye(4)
    right = left_right(sp.Matrix([0, 1, 4]), eye, right=True)

    def operator(vector):
        return left_right(vector, eye, right=False)*right

    axis = sp.ones(3, 1)
    tangents = [operator(3*sp.eye(3)[:, i]-axis) for i in range(3)]
    planes = [(eye+operator(2*axis+sp.eye(3)[:, i])/17)/2 for i in range(3)]
    states = [plane/2 for plane in planes]
    axial = operator(axis)
    assert axial**2 == 51*eye
    mean = sum(states, sp.zeros(4))/3
    assert mean == (eye+7*axial/51)/4
    a = sp.Symbol('a', real=True)
    effects = [(eye+a*t)/3 for t in tangents]
    assert sum(effects, sp.zeros(4)) == eye
    for i, state in enumerate(states):
        tangent = tangents[i]
        assert tangent.T == tangent and tangent**2 == 102*eye
        assert tangent*axial+axial*tangent == sp.zeros(4)
        assert state-mean == tangent/204
        for j, effect in enumerate(effects):
            k = sp.S.One if i == j else sp.Rational(33, 34)
            assert planes[j]*state*planes[j] == k*states[j]
            assert k/3/(sp.Rational(50, 51)) == sp.Rational(1, 3)+(int(i==j)-sp.Rational(1, 3))/100
            assert sp.simplify((state*effect).trace()-(1+a*(3*int(i==j)-1))/3) == 0
    q = (1-7/sp.sqrt(51))/2
    p = 1-q
    assert sp.simplify(p*q) == sp.Rational(1, 102)
    prob = sp.Symbol('P', positive=True)
    lo = sp.sqrt(q*(prob-q))/prob
    assert sp.simplify(lo.subs(prob, 2*q)**2) == sp.Rational(1, 4)
    assert sp.simplify(lo.subs(prob, 1)**2) == sp.Rational(1, 102)
    # Rational unsharp orientation effects are positive exactly when 102*a**2 <= 1.
    assert 102*sp.Rational(1, 11)**2 < 1 < 102*sp.Rational(1, 10)**2
    print('Projection: success 50/51, conditional orientation coefficient 1/100.')
    print('Complete readout: sharp coefficient 1/sqrt(102), correct-label rate', ((1+2/sp.sqrt(102))/3).evalf(14))
    print('Acceptance cutoff for conditional correct-label rate 2/3:', (2*q).evalf(14))


if __name__ == '__main__':
    verify()
