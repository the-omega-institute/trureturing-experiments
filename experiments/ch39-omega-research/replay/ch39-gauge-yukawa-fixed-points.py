"""All stationary squared couplings in the corrected gauge/Yukawa truncation."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb
from sympy import Matrix, Rational as Q, diag, eye, ones, zeros


def model():
    self_part = 10*eye(3)+ones(3)
    mix = Matrix(3,3,lambda i,j: 5*int(i==j or j==(i+1)%3))
    bridge = self_part.row_join(mix).col_join(mix.T.row_join(self_part))
    extra = Matrix([[Q(16,5),3],[1,Q(44,5)]])
    k = diag(*([Q(25,18)]*6+[Q(1,6),Q(1,2)]))
    mat = diag(bridge,extra)
    c5 = Matrix([Q(54,5)]*6+[Q(72,5),Q(108,5)])
    c10 = Matrix([Q(297,20)]*6+[0,0])
    d5 = Matrix([15]*6+[Q(12,5),Q(54,5)])
    d10 = Matrix([5]*6+[0,0])
    r = Q(33,8)
    assert k*c5==d5 and k*c10==r*d10
    assert k*mat==(k*mat).T and (k*mat).is_positive_definite
    one, self10, cross = run_path(str(Path(__file__).with_name("ch39-two-loop-product-beta.py")))["target_inverse_coefficients"]()
    assert str(one)=="94/3" and str(self10)=="22169/30" and r*int(cross)==594
    return mat,c5,c10,d5,d10,k


def verify():
    mat,c5,c10,d5,d10,k = model()
    # Normalize active Yukawa squares by z when w=0. The same linear systems
    # yield the w coefficients used to rule out jointly interacting gauge groups.
    axis = []
    for mask in product((0,1),repeat=8):
        active = [i for i,v in enumerate(mask) if v]
        v5,v10 = zeros(8,1),zeros(8,1)
        if active:
            inv = mat.extract(active,active).inv()
            for col,charge in [(v5,c5),(v10,c10)]:
                values = inv*charge.extract(active,[0])
                for i,value in zip(active,values):
                    col[i] = value
            assert all(v5[i]>0 for i in active)
            assert all(v10[i]>=0 for i in active)
        d = Q(8192,15)-d5.dot(v5)
        h = 594-d5.dot(v10)
        assert h>0 and d>470
        # Simultaneous zeros would require z >= (188/3)/(144 U)
        # and z <= (4/3)/(d U).
        assert Q(188,3)*d > Q(4,3)*144
        assert Q(94,3) > 144*Q(2,3)/d  # w is marginally UV attractive at w=0
        assert all((mat*v5-c5)[i]<0 for i,v in enumerate(mask) if not v)
        axis.append(d)
    assert len(axis)==256 and min(axis)==Q(102963932,217005) and max(axis)==Q(8192,15)
    alpha = lambda b: 8*arb.pi()*int(b.q)/(3*int(b.p))
    assert arb(".01533")<alpha(max(axis))<arb(".01535")
    assert arb(".01765")<alpha(min(axis))<arb(".01767")
    print("Origin and 256 non-Gaussian fixed points; all have g10=0.")
    print("Endpoints of alpha5 on the non-Gaussian axis",alpha(max(axis)),alpha(min(axis)))


if __name__ == "__main__":
    verify()
