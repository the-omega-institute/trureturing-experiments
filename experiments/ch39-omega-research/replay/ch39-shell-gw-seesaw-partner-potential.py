"""Scalar partners for the seesaw effective potential in a constant light background."""
from flint import arb
from sympy import Matrix, Rational, Symbol, exp, log, series, simplify, sqrt

q = arb


def algebra():
    u = Symbol("u", positive=True)
    h = (1+2*u+sqrt(1+4*u))/2
    f = lambda t: t**2*(log(t)-Rational(3,2))
    assert series(f(1+2*u)-f(h),u,0,4).removeO().expand() == -2*u**2+4*u**3
    mass = Matrix([[1+u,sqrt(u)],[sqrt(u),u]])
    assert simplify(mass.det()-u**2)==0 and simplify(mass.trace()-1-2*u)==0
    r = Symbol("r", positive=True)
    soft = Matrix([[1+u,sqrt(u)],[sqrt(u),r+u]])
    assert simplify(soft.det()-(u**2+r*u+r))==0
    # Derivatives of the heavy and light scalar eigenvalues from the quadratic
    # characteristic polynomial at u=0 (dimensionless mass-squared units).
    upper = 1+1/(1-r)
    lower = 1-1/(1-r)
    net = -(upper-2)+r*(log(r)-1)*lower
    assert simplify(net+r*(1+r*log(r)/(1-r)))==0


def numbers():
    Q = q("4.352e13").union(q("4.473e13"))
    weak = q(125)**2/2
    vev = q("246.22")
    mfloor = q(".002421").sqrt()*q("1e-9")
    yfloor = 2*Q.lower()*mfloor/vev**2
    width = 2*weak*16*arb.pi()**2/Q.lower()**2
    print("full allowed width for a single singlet scalar portal",width.upper().str(8),
          "unavoidable |one-loop quartic matching| above", (yfloor**2/(16*arb.pi()**2)).lower().str(8))
    soft = q("7.41e10")
    # c(r)=1+r log(r)/(1-r) decreases on (0,1) by log(r)<=r-1.
    rmax = q(".01")
    factor = 1+rmax*rmax.log()/(1-rmax)
    rem = yfloor*soft**2*factor/(16*arb.pi()**2)
    max_soft = (weak*16*arb.pi()**2/(yfloor*factor)).sqrt()
    print("soft light scalar for m_soft<=M/10: mass limit GeV",max_soft.upper().str(8),
          "cost if m_soft>=7.41e10 GeV",rem.lower().str(8))
    assert width<q("1.31e-21") and rem>q("2.3e18") and max_soft<q("4290")


if __name__=="__main__":
    algebra()
    numbers()
