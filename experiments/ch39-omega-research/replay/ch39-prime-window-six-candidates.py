"""Which of the six center solutions carries most dyadic-window mass?"""
from itertools import combinations
from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    h, log_two = arb(2).sqrt(), arb(2).log()
    heights = [z.imag for z in acb.zeta_zeros(11485,6)]
    cosines = [(g*log_two).cos() for g in heights]
    intercept = [((3*h-4*v)/(g*g+arb(1)/4))**2 for v,g in zip(cosines,heights,strict=True)]
    slopes = [a*v for a,v in zip(intercept,cosines,strict=True)]
    winners = [1,2,0,3,4]
    breaks = [arb(0)]
    for i,j in zip(winners,winners[1:]):
        crossing = (intercept[i]-intercept[j])/(slopes[j]-slopes[i])
        assert breaks[-1] < crossing < 1
        breaks.append(crossing)
    breaks.append(arb(1))
    for row,i in enumerate(winners):
        def score(j,v):
            return intercept[j]+v*slopes[j]
        for j in range(6):
            if j==i:
                continue
            lo, hi = breaks[row:row+2]
            if row>0 and j==winners[row-1]:
                assert slopes[i] > slopes[j] and score(i,hi)>score(j,hi)
            elif row+1<len(winners) and j==winners[row+1]:
                assert slopes[i] < slopes[j] and score(i,lo)>score(j,lo)
            else:
                assert score(i,lo)>score(j,lo) and score(i,hi)>score(j,hi)
    kbreaks = [v/(h*(1+(1-v*v).sqrt())) for v in breaks]
    print('Small-branch switches', [x.str(12) for x in kbreaks[1:-1]])

    k = arb(4)/9
    v = 2*h*k/(1+2*k*k)
    weights = [(a+v*b).sqrt() for a,b in zip(intercept,slopes,strict=True)]
    order = [3,4,0,2,5,1]
    assert all(weights[i]>weights[j] for i,j in zip(order,order[1:]))
    tolerance = (weights[3]-weights[4])/(weights[3]+weights[4])
    assert arb('.01280') < tolerance < arb('.01281')
    assert arb('.001467') < tolerance
    fraction = weights[3]/sum(weights)
    assert arb('.19429') < fraction < arb('.19430')
    # The envelope bound applies to each positive index, with multiplicity already enumerated.
    a,b = 1+2*k*k,2*h*k
    envelope = lambda z: (3*h-4*z)*(a+b*z).sqrt()
    peak = (3*h*b-8*a)/(12*b)
    assert -1<peak<1 and envelope(peak)>max(envelope(arb(-1)),envelope(arb(1)))
    first,second = [z.imag for z in acb.zeta_zeros(1,2)]
    first_weight = envelope((first*log_two).cos())/(first**2+arb(1)/4)
    rest_upper = envelope(peak)/(second**2+arb(1)/4)
    assert first_weight*(1-arb('.001467')) > rest_upper*(1+arb('.001467'))
    curve = [1+(1+(g/312)**2).sqrt() for g in heights]
    separation = min(curve[3]-curve[2],curve[4]-curve[3])/2
    assert arb('36.71034083') < curve[3] < arb('36.71034084')
    assert arb('.00074281') < separation < arb('.00074282')

    for n in range(163,366):
        k = arb(n)/(624-n)
        v = 2*h*k/(1+2*k*k)
        if n<=204 or n>=317:
            assert breaks[3]<v<breaks[4]
        else:
            assert v>breaks[4] or v<breaks[3]
    print('Band winner 11488, fraction',fraction.str(12),'relative tolerance',tolerance.str(12))
    print('Global winner index 1. Curvature',curve[3].str(14),'identification half-gap',separation.str(12))


if __name__=='__main__':
    verify()
