"""Polar transport current and volume-independent weighted resolvent bounds."""
from fractions import Fraction as F
import sympy as sp


def polar_current():
    eye = sp.eye(2)
    v = sp.Matrix([12, 5])/13
    w = sp.Matrix([24, 7])/25
    dv = sp.Matrix([0, 5*sp.I])/13
    dw = sp.Matrix([0, -7*sp.I])/25
    P, Q = v*v.H, w*w.H
    dP, dQ = dv*v.H+v*dv.H, dw*w.H+w*dw.H
    T = P*Q+(eye-P)*(eye-Q)
    dT = dP*(2*Q-eye)+(2*P-eye)*dQ
    B = T.H*T
    dB = dT.H*T+T.H*dT
    b, db = B[0, 0], dB[0, 0]
    assert B == b*eye and dB == db*eye and B == eye-(P-Q)**2
    O = T/sp.sqrt(b)
    dO = dT/sp.sqrt(b) - T*db/(2*b**sp.Rational(3, 2))
    Y = Q*B.inv()*T.H
    actual = sp.im((Q*O.H*dO).trace())
    assert actual != 0 and actual == sp.im((Y*dT).trace())
    assert actual == sp.im((((2*Q-eye)*Y)*dP + (Y*(2*P-eye))*dQ).trace())
    assert O.H*O == eye and O*Q*O.H == P


def higgs_generator():
    phi = sp.eye(5)[:, 0]
    tangent = sp.Matrix([3*sp.I, 4, 0, 0, 0])/5
    phase = (phi.H*tangent)[0]
    normal = tangent-phase*phi
    line = phi*phi.H
    generator = normal*phi.H-phi*normal.H + phase*(line-(sp.eye(5)-line)/4)
    assert generator*phi == tangent and generator.H == -generator and generator.trace() == 0
    assert 8*2*2**45 + 60*4*896 < 2**50


def constants():
    # Verify the coefficient using positivity in the variable delta^2 <= 1.
    unit_step = F(1, 2**18)
    assert 10*((1+unit_step)**4-1) < F(1, 4096)
    assert F(7,2048) < F(1,256) and F(1,4)+F(2,256) < F(1,3)
    assert 1/(1-F(1,3)**2) < 2
    p_loss = 7*32*4
    inverse_loss = 4*F(2,3)*(2*p_loss)
    t_loss = 10*p_loss
    y_loss = 26*p_loss+26*inverse_loss+4*t_loss
    c_loss = 104*p_loss+5*y_loss
    product_loss = 2*128*260*2 + 4*c_loss
    assert 2*7*240*2048*product_loss < 2**45
    assert (44+7+1)*32 < 2048
    assert F(3*32**2, 256**2) < F(1, 16)  # sqrt(3)*32*a < delta/2 in the large-coefficient example
    print("Exact polar-current identity and 2^45 locality prefactor verified.")


if __name__ == "__main__":
    polar_current()
    higgs_generator()
    constants()
