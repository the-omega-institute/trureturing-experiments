"""CKM obstruction after bounded endpoint distortion, for all common left singular bases."""
from flint import arb

condition_squared=arb("1.809")
mu=(condition_squared-1)/(condition_squared+1)
up=arb(".00316")/arb(".8573")
down=arb(".000271")/arb(".01392")


def direction(r,e):
    den=1-2*mu-(1+mu)*(r+e)
    assert den>0
    return ((2*(1+mu)).sqrt()*mu*(r+e)+(1-2*mu)*e)/den


def ceilings(total):
    s=max((direction(up,e)+direction(down,total-e)).upper() for e in (arb(0),total))
    return s,mu+arb(2).sqrt()*s


for budget in (arb(0),arb(".0035")):
    a,b=ceilings(budget)
    assert a<arb(".035") and b<arb(".98")
    print("normalized noncommuting correction, same-mode Vcb and distinct-mode Vtb",
          budget.str(5),a.str(12),b.str(12))
print("sum of relative noncyclic kernel norms and relative commutator norms must exceed",
      (arb(".0035")/condition_squared).str(12),
      (arb(".0035")*3/(2*condition_squared)).str(12))
