"""Uniformly synthesize the dyadic phase from the shell Cayley coordinate."""
from math import prod
from pathlib import Path
from flint import acb, arb, arb_poly, arb_series, ctx

# Eight positive nodes and weights for the omitted odd poles, in units of 2^-160.
QUADRATURE = [
    (11771794845886840877817047808581836621856353573, 2367108840908316590904553320080098212464164616),
    (104483973755039697500284119678813109992930587364, 2076330151016903062504838391868253072275105630),
    (280432621623710640291032317772834592786268076093, 1670610740665150692238107339770558247300854588),
    (517288706669619191626708970588920070793086111394, 1280548434893349424188690121072362479125751999),
    (781262943480170371739760321465178767824439289218, 942473373550056135208125399803406135725973008),
    (1033229060180711011527667751851876954677615059479, 656824973795791777633765175044832682205389738),
    (1247453261736273691561056674369814810925823043178, 479124805012370686006606792398140455293385571),
    (1461397655663372280677333877996356860619001378570, 493900362394829198136487199282374656399645492),
]
ANCHOR = 11674
SPAN = 121**2 - ANCHOR
SLOPE, SCALE, ROTATION = (arb(n)/2**50 for n in (248413601658037, 1433540284805665, 689006131071067))


def cayley(tangent):
    return acb(1, tangent) / acb(1, -tangent)


def dyadic_response(height, nodes, weights):
    u = SLOPE*(height-11107)
    q = (u**2-ANCHOR)/SPAN
    terms = [(1, k*k-u*u) for k in range(1, 120, 2)]
    terms.extend((v, 1-q*x) for x, v in zip(nodes, weights))
    numerator, denominator = arb(0), arb(1)
    for weight, pole_distance in terms:
        numerator = numerator*pole_distance + weight*denominator
        denominator *= pole_distance
    scaled = SCALE*u*numerator
    return acb(1, 1)*cayley(ROTATION)*acb(denominator, scaled)/acb(denominator, -scaled)


def verify():
    ctx.prec, ctx.cap = 448, 20
    nodes, weights = zip(*[(arb(x)/2**160, arb(v)/2**160) for x,v in QUADRATURE])
    assert all(0 < x < 1 and v > 0 for x,v in zip(nodes, weights))
    assert all(x<y for x,y in zip(nodes,nodes[1:]))
    s = arb_series([ANCHOR, SPAN]).sqrt()

    def digamma(series):
        return series.lgamma().derivative()/series.derivative()

    moments = (digamma((121+s)/2)-digamma((121-s)/2))/(4*s)
    residual = max(abs(moments[k]-sum(v*x**k for x,v in zip(nodes,weights))) for k in range(16))
    assert residual < arb(2)**-152
    polynomial = arb_poly([1])
    for x in nodes:
        polynomial *= arb_poly([-x,1])
    norm = sum(c*moments[k] for k,c in enumerate((polynomial**2).coeffs()))
    radius = arb("118.49")
    q_min, q_max = -arb(ANCHOR)/SPAN, (radius**2-ANCHOR)/SPAN
    floor = 1-q_max
    max_q = max(abs(q_min), q_max)
    # Newton form of the degree-15 Hermite interpolant at duplicated nodes.
    hermite = sum((2*max_q/floor)**k/floor for k in range(16))
    remainder = max(prod(abs(q)/(1-q*x) for x in nodes)**2
                    * norm/(1-max(arb(0),q)) for q in (q_min,q_max))
    log_two, pi = arb(2).log(), arb.pi()
    assert 537*SLOPE < radius
    offset_error = abs(acb(1,1)*cayley(ROTATION)/arb(2).sqrt()-acb(0,11107*log_two).exp())
    # Multiplying a real tangent by k>0 changes its Cayley phase by <= |k-1|/sqrt(k).
    ratio = SCALE*pi/4
    bound = (offset_error + 537*abs(pi*SLOPE-log_two) + abs(ratio-1)/ratio.sqrt()
             + 2*SCALE*radius*(residual*hermite + remainder))
    assert bound < arb("5.2e-7")
    assert 1074*log_two-2*arb("1e-6").asin() > 2*118*pi
    for odd in (1, 117):
        at_pole = dyadic_response(11107+odd/SLOPE,nodes,weights)
        assert abs(at_pole+acb(1,1)*cayley(ROTATION)) < arb("1e-90")
    print("Normalized dyadic phase error <", bound.upper().str(8), "continuous order 119 to 136")

    # Reuse the independently certified 48-bit shell coordinate for all possible labels.
    # The dyadic column is deliberately ignored by the ranking.
    table = Path(__file__).with_name("ch39-omega-rational-register.tsv")
    scores = {}
    for line in table.read_text().splitlines()[1:]:
        n, k, _, _ = map(int, line.split())
        height = arb(2)**47/k
        assert 10570 < height < 11644
        d = dyadic_response(height,nodes,weights)
        scores[n] = abs((d-1)*(d-2)*(9+4*d)/13)/(height**2+arb(1)/4)
    assert list(scores) == list(range(10855,12093))
    gap = scores[10875]-max(s for n,s in scores.items() if n!=10875)
    assert gap > arb("2e-11")
    print("Single 48-bit coordinate, 1238 labels, winning weight margin", gap.str(12))


if __name__ == "__main__":
    verify()
