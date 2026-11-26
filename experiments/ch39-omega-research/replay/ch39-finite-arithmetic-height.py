"""A finite Dirichlet readout and the winning binary height, with a uniform tail bound."""
from math import factorial
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

source = run_path(str(Path(__file__).with_name("ch39-hardy-height-selection.py")))
weight, hardy = source["weight"], source["hardy"]


def finite_value(g):
    s = acb(arb(1)/2,g)
    n, order = 4096, 32
    log_n = arb(n).log()
    divisor = list(range(n+1))
    for prime in range(2,65):
        if divisor[prime] == prime:
            for j in range(prime*prime,n+1,prime):
                divisor[j] = min(divisor[j],prime)
    powers = [acb(0), acb(1)]
    for i in range(2,n+1):
        p = divisor[i]
        powers.append((-s*arb(i).log()).exp() if p == i else powers[p]*powers[i//p])
    prefactor = (-s*log_n).exp()
    result = sum(powers[1:n])+prefactor*(n/(s-1)+arb(1)/2)
    rising = s
    for k in range(1,order+1):
        result += prefactor*arb.bernoulli(2*k)*rising/(factorial(2*k)*n**(2*k-1))
        rising *= (s+2*k-1)*(s+2*k)
    return result


def verify():
    ctx.prec = 240
    root = acb.zeta_zero(10875).imag
    best = weight(root)
    strength = arb(1)/256
    n, order = 4096, 32
    factor = (arb(11644)**2+(2*order)**2).sqrt()/(2*arb.pi()*n)
    epsilon = arb("3.65e-22")
    assert 4*arb(n).sqrt()*factor**(2*order)/(2*order-arb(1)/2) < epsilon
    assert abs(finite_value(root)) < epsilon
    # Independently certify a uniform margin outside the sharp neighborhood used for the grid.
    ctx.prec, ctx.cap = 96, 2
    near = arb(1)/1024
    stack = [(arb(10603),(root-near).upper()),((root+near).lower(),arb("10621.28"))]
    slack = arb("1e-10")
    assert source["modulus"](-arb(59)/108)/(arb("10621.28")**2+arb(1)/4) < best*(1-slack)
    while stack:
        lo,hi = stack.pop()
        cell = lo.union(hi)
        gain = weight(cell)/best
        if gain < 1-slack or gain < (1-slack)*(1+strength*(abs(hardy(cell))-epsilon)):
            continue
        assert hi-lo > arb("1e-9")
        middle = (lo+hi)/2
        stack.extend([(lo,middle),(middle,hi)])
    ctx.prec = 240
    # Any continuous maximizer beats the exact root evaluated with the finite response.
    radius = 2*strength*epsilon/(5*strength-arb(".015")*(1+strength*epsilon))
    assert radius < arb("6.30e-22")
    grid = 2**64
    center = int((root*grid).floor().unique_fmpz())
    points = range(center-10,center+12)
    scored = {m: weight(arb(m)/grid)/(1+strength*abs(finite_value(arb(m)/grid))) for m in points}
    winner = max(scored,key=lambda m:scored[m].mid())
    assert winner == 195923358475422178854846
    assert all(scored[winner]>v for m,v in scored.items() if m!=winner)
    # If a different grid height matches this score, the linear Z and W bounds force it nearby.
    lower = scored[winner].lower()/best
    assert lower > 1-slack
    max_distance = (1-lower+lower*strength*epsilon)/(5*lower*strength-arb(".015"))
    assert 0 < max_distance < 10/arb(grid)
    gap = (scored[winner]-max(v for m,v in scored.items() if m!=winner))/best
    print("64-bit height numerator",winner,"error", (arb(winner)/grid-root).str(12),
          "relative winning margin among checked neighbors",gap.str(10))
    print("Continuous maximizers within",radius.upper().str(9),"of zero 10875")


if __name__ == "__main__":
    verify()
