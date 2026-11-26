"""Integer sign proof for the transition between the seventh offsets 17 and 18."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
lib = run_path(str(root/'ch39-last-offset-closure.py'))
old = lib['old']
B, node, mul, rounded, root_cell = (old[k] for k in ('B','node','mul','rounded','root_cell'))
p, modes, q2, target = old['base']['physical_inputs']()
mults = old['base']['MODES']
poles = [B(1)/node(i) for i in range(1,5)]
top = poles[2]
H = B(125**2)*top/q2
data = old['base']['data']('ch39-last-offset-witnesses.json')


def response(table, y):
    return sum((rounded(w/(B(1)-B(y)*v)) for v,w in zip(*table)), B(0))


def root_points(rows):
    assert len(rows)==3
    for j,row in enumerate(rows):
        for raw,lam,n in zip(row,modes,mults,strict=True):
            y = Q(raw,10**5 if j else 10**12)
            assert (poles[j-1].hi if j else H.hi) < y < poles[j].lo
            yield y,lam,n


def check_seventeen():
    ids = [1,2,3,4,8,9,52,53]
    a = lib['interpolate'](ids,17)
    d = lib['interpolate']([0,1,2,3,5,6,16,17],17)
    vs,ws = a
    assert min(ws[0].lo,ws[2].lo,d[1][2].lo)>Q('.05')
    assert ((B(1)-p)/(B(1)-H*vs[0])).hi < 1
    squares = [1,4,9,16,25,36,289]
    def factor(y):
        return B(y)*mul(*(B(y+b) for b in squares))/mul(*(B(1)-B(y)*v for v in vs))
    cap = Q('1.93e-26')
    at = Q(320)
    assert factor(at).hi<0
    assert ((response(d,at)-response(a,at))/factor(at)).hi < cap
    slope = mul(*(B(1)+B(b)*vs[0] for b in squares))/(vs[0]*mul(*(vs[0]-v for v in vs[1:])))
    assert len(data['seventeen'])==4
    for k, rows in enumerate(data['seventeen']):
        lower = k*cap/4
        upper = (k+1)*cap/4
        scale = vs[0]*old['b41']*(ws[0]+slope*B(lower))
        exposure = mul(*((top/v)**25 for v in poles[:2]))
        for y,lam,n in root_points(rows):
            remainder = B(lower,upper/(1-y*node(5).hi))
            reading = response(a,y)+factor(y)*remainder
            assert (scale*(B(y)-H)-lam+lam*(B(1)-p)*reading).lo>0
            exposure = mul(exposure,(top/B(y))**n)
        ratio = rounded(exposure**11/target)
        assert ratio.lo>Q('1.002')
        print('offset 17, part', k, 'lower ratio',float(ratio.lo))


def check_counterexamples():
    c = B(Q('.000834926'))
    cosines = old['cosines']
    loss = 1-root_cell((B(1)+cosines[41])/B(2)).lo
    sine = [B(0)]+[root_cell((B(1)-v)/B(2)) for v in cosines[1:]]
    increment = min((sine[j]-sine[j-1]).lo for j in range(1,42))
    assert [r['j'] for r in data['witnesses']]==list(range(18,232))
    high = B(650)
    H650 = B(125**2)*high/q2
    assert top.hi<650<poles[3].lo and (c*high).hi<8
    assert ((B(1)-p)/(B(1)-H650*node(1))).hi<1
    upper = [[2060,1980,1332,2036,1590,2016],
             [339,336,299,338,315,337],[569,566,526,568,546,568]]
    ratio_max = 0
    for record in data['witnesses']:
        j = record['j']
        k,l = record['tail_pairs']
        assert 4<k<k+1<l<l+1<=len(old['prefix'])
        v = record['floor256']
        cell = B(Q(v,2**256),Q(v+1,2**256))
        table = lib['interpolate']([1,2,3,4,k,k+1,l,l+1],j,last_cell=cell,error=Q(0))
        vs,ws = table
        assert min(w.lo for w in ws)>0
        edges = [root_cell(w*z) for w,z in zip(ws,vs)]
        first,other = edges[0],max(s.hi for s in edges[1:])
        length = sum(s.lo for s in edges)
        assert first.lo>other and first.lo**2>c.hi and first.hi**2*old['b41'].hi<c.lo
        assert length>2*first.hi+other*loss and length>2*other+first.hi*loss
        assert first.lo*increment>other*sine[1].hi
        exposure = mul(*((top/t)**25 for t in poles[:2]))
        for y,lam,n in root_points(record['roots']):
            f = c*(B(y)-H)-lam+lam*(B(1)-p)*response(table,y)
            assert f.hi<0 and (q2/high).lo*y>1000**2
            exposure = mul(exposure,(top/B(y))**n)
        ratio = rounded(exposure**11/target)
        assert ratio.hi<1
        ratio_max = max(ratio_max,ratio.hi)
        exposure = mul(*((high/t)**25 for t in poles[:3]))
        for gap,row in enumerate(upper):
            for v,lam,n in zip(row,modes,mults,strict=True):
                y=Q(v,10**6 if gap==0 else 1)
                assert (poles[gap-1].hi if gap else 0)<y<poles[gap].lo
                assert (c*(B(y)-H650)-lam+lam*(B(1)-p)*response(table,y)).lo>0
                exposure = mul(exposure,(high/B(y))**n)
        assert (exposure**11/target).lo>1
        for lam in modes:
            assert (c*(high-H650)-lam+lam*(B(1)-p)*response(table,Q(650))).hi<0
    print('Offsets 18 through 231 admit 144 open doublets. Maximum boundary upper ratio',float(ratio_max))


def check_source():
    from flint import arb, ctx
    origin = run_path(str(root/'ch39-extra-offset-cutoffs-origin.py'))
    mass_numerator = old['base']['data']('ch39-natural-mass-cell.json')
    ctx.prec, ctx.cap = 560, 3
    mass = arb(mass_numerator).union(arb(mass_numerator+1))/2**288
    for row in data['witnesses']:
        value = origin['level'](row['j'],mass)*2**256
        assert int(value.floor().unique_fmpz())==row['floor256']
    print('All 214 independent RH source cells checked with Arb.')


if __name__=='__main__':
    import sys
    if '--origin' in sys.argv: check_source()
    check_seventeen()
    check_counterexamples()
