"""Finite observed CKM and four quark ratios for the actual common Maxwell zero-mode metric."""
from itertools import product
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat

finite=runpy.run_path(str(Path(__file__).with_name("ch39-finite-orthogonal-ckm.py")))
fit,common=finite["fit"],finite["common"]
adj=lambda m:m.conjugate().transpose()
ncell=4096
qlo,qhi=arb(".855"),arb(".925")
step=(qhi-qlo)/ncell
F=acb_mat([[fit["omega"]**(-i*j)/arb(3).sqrt() for j in range(3)] for i in range(3)])
model=fit["model"](*common["angles"])
source=[sum((h*adj(h) for h in orbit),acb_mat(19,19)) for orbit in model[1]]
# Coefficients pick explicit stiffnesses. Every cell center and its full metric tolerance
# are subsequently checked with complex balls and Hermitian spectral estimates.
raw="""
0.000201012686358921 1.10824561760408e-19 -9.17415842466489e-20 9.9365549215183e-20 -6.35909663872808e-20 7.23170821320333e-20
0.0145464049637832 2.57161130130669e-05 -6.23069171715494e-07 1.77260368840631e-08 -4.52301005357976e-09 1.1819167931004e-10
0.158423001350456 0.000929050440711829 -6.43594563846525e-06 2.24781726084396e-07 -2.28616387387507e-08 6.80017771050487e-10
0.286890108517912 0.000284465244780425 -1.83270118688847e-06 6.52637694565678e-08 -5.78724324710107e-09 1.77327464216269e-10
0.976042481516899 0.0012677829664349 -1.3502475368214e-05 2.85544487250914e-06 -5.21183005703277e-07 2.83458209843767e-08
0.315534656964864 0.00804911529316249 -5.27043435670397e-05 1.85383267673769e-06 -1.61332160989286e-07 4.9864044022957e-09
0.315538135691804 0.00804913318262595 -5.27044609182482e-05 1.85383680181939e-06 -1.61332521089033e-07 4.98641551710341e-09
1.04269567641344 0.0159721682446085 8.36157844264776e-05 -1.95808946331342e-06 1.2006284629719e-07 -3.50374981729837e-09
0.000333819579941842 2.81534979347476e-19 -2.08868972651491e-19 2.52424890311104e-19 -1.8066326105755e-19 1.83711876691447e-19
0.554979071177063 -0.00145414414170724 0.00011092666903616 -4.4441779927234e-06 6.96017277188914e-07 -2.87489371613656e-08
0.0326887382749005 0.00271230627141198 6.66191574313833e-05 1.08869276150436e-06 4.76861482177438e-07 8.02456608562223e-09
0.000144126511839488 6.64508750991452e-20 -4.06760681993017e-20 5.95800028005614e-20 -6.40798226808684e-20 4.33616277471006e-20
1e-06 1.43472570228597e-23 3.34148777627268e-22 1.28637826413419e-23 -1.36918832548765e-23 9.36211023991806e-24
0.10460184737121 -0.00825832776554343 0.000165467211389654 -1.39401061173167e-06 5.73198891244361e-07 6.44561041561393e-09
2.17414759058164 0.0313212080429928 0.000170902802362495 4.17679978296074e-07 -4.98001675280074e-07 1.43041937501694e-08
0.0625905981178776 0.00374993756810751 1.71698976193494e-05 4.34616157931623e-07 7.20307552702768e-08 2.39100324655244e-09
0.810257231146456 0.0290814508566019 4.72853459556699e-05 9.22353788654902e-07 -2.31087479619165e-07 -5.42570849640235e-09
3.40955773542016e-06 9.71861196757184e-22 -9.68671433448179e-23 8.71372916280749e-22 -1.04837638282077e-21 6.34174995780888e-22
0.315270664736033 0.0116824940206874 3.27530527541291e-05 -8.2099355029972e-07 -1.46735024598667e-07 -7.18544145386969e-09
0.343747748934028 -0.000872109356485189 7.9669052962888e-05 -3.30778974990091e-06 6.67432167219251e-07 -2.03928161766723e-08
"""
coefficients=[list(map(arb,row.split())) for row in raw.strip().splitlines()]


def metric_interval():
    osc=arb(".04822")
    limits=[(arb(lo)*(-osc/2).exp(),arb(hi)*(osc/2).exp())
            for lo,hi in (("7.70","7.93"),("4.73","4.88"),("6.08","6.33"))]
    def verify(box,d):
        a,b,c=[lo.union(hi) for lo,hi in box]
        assert a>c>b
        # Mix (a,b,c) and (c,a,b) so the last two coordinates coincide.
        t=(c-b)/(a+c-2*b)
        q=(((1-t)*b+t*a)/((1-t)*a+t*c)).sqrt()
        if qlo<q<qhi:
            return 1
        assert d<27
        j=d%3
        lo,hi=box[j]
        mid=(lo+hi)/2
        count=0
        for part in ((lo,mid),(mid,hi)):
            new=list(box)
            new[j]=part
            count+=verify(new,d+1)
        return count
    return verify(limits,0)


def weight_row(q):
    x=(q-arb(".89"))/arb(".035")
    terms=[arb(1),x]
    for _ in range(4):
        terms.append(2*x*terms[-1]-terms[-2])
    w=[sum(a*r for a,r in zip(row,terms)) for row in coefficients]
    assert len(w)==20 and min(w)>0
    return w


def kernel(w):
    wall,J,block=fit["W"],fit["J"],fit["block"]
    m=model[0]+sum((val*v for val,v in zip(w[:4],source)),acb_mat(19,19))
    k=block(m,wall,wall)-block(m,wall,J)*block(m,J,J).solve(block(m,J,wall))
    for end,a in product((1,2),range(3)):
        finite["add_edge"](k,fit["inner"][a],fit["outer"][(a+end)%3],w[4+3*(end-1)+a])
    finite["fix_phase"](k,(7,6,8),arb.pi()-arb(".002"))
    finite["fix_phase"](k,(11,12,18),arb(".03606860"))
    degree,margin=finite["flux"]["degrees"](k)
    assert degree==(3,[3]) and margin>arb(".0002")
    t=-block(k,fit["ii"],fit["oo"])
    return adj(F)*(t*t)*F


def center_test(q,w):
    diag=acb_mat([[q,0,0],[0,1,0],[0,0,1]])
    u,d=[diag*kernel(part)*diag for part in (w[:10],w[10:])]
    vals=finite["observables"](u,d)
    targets=common["up"]+list(reversed(common["down"]))
    assert all(abs(a/b-1)<arb(".0001") for a,b in zip(vals[:4],targets))
    assert vals[0]/vals[1]<arb(".052") and vals[2]/vals[3]<arb(".052")
    for value,ref in zip(vals[4:8],(".22501",".04183",".003732","1.147")):
        assert abs(value-arb(ref))<arb(".000002")
    return vals


def extend():
    # Parameter q changes the two Gram factors by the same positive diagonal congruence.
    dilation=(1+step/(2*qlo))**2
    assert dilation*arb("1.0001")<arb("1.00013")
    assert dilation*arb(".052")<arb(".06")
    # Parallel-transport each singular vector. A rank-one unit projection has
    # off-diagonal entries of absolute value <=1/2 in each orthonormal basis.
    slope=(arb(2).sqrt()/qlo)*(1+arb(".06"))/(1-arb(".06"))
    dv=slope*step/2
    assert dv<arb(".0000161")
    s12,s23,s13,delta=map(arb,(".22501",".04183",".003732","1.147"))
    err=arb(".000002")
    angles=[s+arb(0,err) for s in (s12,s23,s13)]
    a,b,c=angles
    c13=(1-c*c).sqrt()
    abs_center=[a*c13,b*c13,c,(1-a*a).sqrt()*(1-b*b).sqrt()-a*b*c]
    abs_low=[v.lower()-dv for v in abs_center]
    assert min(abs_low)>0
    new_c=c+arb(0,dv)
    new_c13=(1-new_c*new_c).sqrt()
    sa=(a*c13+arb(0,dv))/new_c13
    sb=(b*c13+arb(0,dv))/new_c13
    assert max(abs(sa-s12).upper(),abs(sb-s23).upper(),abs(new_c-s13).upper())<arb(".000019")
    def phase_offset(x,y,z):
        return (x*y*z/((1-x*x)*(1-y*y)).sqrt()).asin()
    phi=err+phase_offset(a,b,c)+dv*sum(1/v for v in abs_low)+phase_offset(sa,sb,new_c)
    assert phi<arb(".006")
    full_delta=delta+arb(0,phi.upper())
    j=sa*sb*new_c*((1-sa*sa)*(1-sb*sb)).sqrt()*(1-new_c**2)*full_delta.sin()
    assert arb(".0000308")<j<arb(".0000316")
    return dv,phi,j


if __name__=="__main__":
    boxes=metric_interval()
    for i in range(ncell):
        q=qlo+(arb(i)+arb(1)/2)*step
        center_test(q,weight_row(q))
    dv,phi,j=extend()
    print("moment boxes, q cells, CKM entry error, phase error, signed J",
          boxes,ncell,dv.upper().str(8),phi.upper().str(8),j.str(8))
