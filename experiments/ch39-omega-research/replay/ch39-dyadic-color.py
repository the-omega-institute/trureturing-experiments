"""Finite rational audit of color harmonics. Use --audit-targets to regenerate the input intervals."""
from fractions import Fraction as F
from functools import cache

COUNTS = ((1,1,23,2,29),(1,1,21,4,30),(1,1,19,8,31),(1,1,34,10,49))
GRID = 1 << 160
# Each record holds lower endpoints in units of 1e-20 and four rows of integer seeds.
DATA = {
    11485: (
        (86093000604017441251, 56229610094183775512, 28917806861315791956, 12123455673967549264, 4263146914078590199),
        (
            (816603599293186, 276370964377441, 240541098997117, 174495538482501, 132540657254920),
            (820917190953328, 301567186430891, 241175992162826, 210633271361956, 132723610461955),
            (833616567742356, 337024332133444, 243465476383734, 211898886747865, 132250515433734),
            (915776971617562, 593106008463193, 244347225985709, 206042339559303, 132209712862343),
        ),
    ),
    11486: (
        (86091008104231048071, 56225073612866928525, 28913178517705443695, 12120387167345733910, 4261633633173088807),
        (
            (816843346683177, 277324411979280, 240528678498074, 174754207089621, 132552403475176),
            (821152651824876, 302265175829279, 241148127667752, 210889142432619, 132737771289525),
            (833853039354988, 337458593162200, 243460827951342, 212006263956522, 132265028406513),
            (915926085029330, 593529972880925, 244382247558695, 206089725122321, 132225010989606),
        ),
    ),
    11487: (
        (86087965849254175642, 56218147967101398144, 28906113583839731459, 12115703910014308564, 4259324351375442948),
        (
            (817209270231203, 278731094099764, 240511992245878, 175146180504257, 132570440100481),
            (821512047926644, 303314809756626, 241106021274820, 211281446824852, 132759415059263),
            (834214047089691, 338118484009160, 243453854221023, 212170320906775, 132287193699226),
            (916152872235862, 594178342639437, 244435723968217, 206162075894727, 132248368195032),
        ),
    ),
    11488: (
        (86086656644271259302, 56215167915675110175, 28903073952559878675, 12113689222793713488, 4258331048720306609),
        (
            (817366693269925, 279319722352630, 240505587980825, 175313822536942, 132578238399983),
            (821666667087731, 303760782357689, 241088043489125, 211450922162770, 132768737364873),
            (834369384817671, 338401306441221, 243450897214873, 212240963402636, 132296734567937),
            (916250136244186, 594457750153838, 244458738356908, 206193211356692, 132258419447367),
        ),
    ),
    11489: (
        (86083890877099721881, 56208873053041298169, 28896653938931206775, 12109434484303478898, 4256233580473661512),
        (
            (817699165503031, 280533084228463, 240493459038464, 175666052536167, 132594781764813),
            (821993223468918, 304692127755749, 241050301273067, 211810339143395, 132788446966528),
            (834697507846439, 338996543986089, 243444733297297, 212390287671634, 132316894574171),
            (916454951776971, 595048784930685, 244507360135401, 206258986782359, 132279652683217),
        ),
    ),
    11490: (
        (86082584126899365759, 56205899211846958251, 28893621304636275214, 12107424895552238503, 4255243026136893918),
        (
            (817856207132673, 281093049406959, 240488347092453, 175831607667362, 132602629658821),
            (822147472678864, 305127244202458, 241032565432193, 211980848266727, 132797766527432),
            (834852519936962, 339276735454448, 243441858664650, 212460882639951, 132326421691107),
            (916551408544185, 595328398353531, 244530333844106, 206290063935629, 132289684536431),
        ),
    ),
}


def decrement(j):
    x = F(1,3*2**26)
    cos = (1-x*x)/(1+x*x)
    prev,cur = F(1),cos
    for _ in range(1,j):
        prev,cur = cur,2*cos*cur-prev
    return 1-cur


@cache
def response(j, tick, bits):
    step = 1-decrement(j)*F(tick,2**bits)
    assert 0 < step <= 1
    num,den = step.as_integer_ratio()
    lo, hi = GRID*num//den,-(-GRID*num//den)
    for _ in range(56):
        lo,hi = lo*lo//GRID,-(-hi*hi//GRID)
    return lo,hi


def error_bounds(bits,model):
    counts = COUNTS[model]
    low_error,high_error = F(0),F(0)
    for cells,seeds in DATA.values():
        shift = 52-bits
        ticks = [(k+2**(shift-1))//2**shift if shift else k for k in seeds[model]]
        assert all(0 <= t <= 2**bits for t in ticks)
        for j,target in enumerate(cells,1):
            intervals = [response(j,t,bits) for t in ticks]
            lo,hi = [F(sum(c*v for c,v in zip(counts,col)),GRID*sum(counts))
                     for col in zip(*intervals)]
            diff_lo,diff_hi = lo-F(target+1,10**20),hi-F(target,10**20)
            low_error = max(low_error,diff_lo,-diff_hi)
            high_error = max(high_error,abs(diff_lo),abs(diff_hi))
    return low_error,high_error


def audit_targets():
    from pathlib import Path
    from runpy import run_path
    from flint import acb,arb,ctx
    earlier = run_path(str(Path(__file__).with_name("ch39-fifth-shell-harmonic.py")))
    two,q = earlier["two"],earlier["q"]
    ctx.prec = 420
    scale = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/scale)**2/(8*arb.pi())
    for zero,(cells,_) in DATA.items():
        targets = two["harmonics"](acb.zeta_zero(zero).imag,power,190,5)
        assert all(c < v*10**20 and v*10**20 < c+1 for c,v in zip(cells,targets))
    print("All 30 rational target intervals regenerated.")


def verify():
    tol = F(1,10**12)
    for model,(counts,expected) in enumerate(zip(COUNTS,(42,42,42,41))):
        for bits in range(1,53):
            low,high = error_bounds(bits,model)
            if high < tol:
                assert bits == expected
                print(sum(counts),"first rounding precision under 1e-12:",bits)
                break
            assert low > tol
        else:
            raise AssertionError("No certified rounding precision")
    for size,expected in ((56,44),(57,44),(60,42),(95,45),(280,46),(285,46)):
        def variation(bits):
            r = 2**bits % size
            return F(r*(size-r),size*2**bits)
        assert variation(expected) < tol and all(variation(b) >= tol for b in range(expected))


if __name__ == "__main__":
    import sys
    if "--audit-targets" in sys.argv:
        audit_targets()
    verify()
