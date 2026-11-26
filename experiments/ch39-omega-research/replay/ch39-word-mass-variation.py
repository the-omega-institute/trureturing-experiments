"""达成 7.2% 末阶份额的最小字质量总变差与最大逐字质量下界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


TV_SHIFTS = {
    ('direct', 1): (
        (0, 44367185), (2, 87996412), (4, 85786559), (6, 82115913),
        (8, 77003070), (10, 70473929), (12, 62561567), (14, 53306063),
        (16, 42754300), (18, 30959730), (20, 17982011), (22, 3886752),
        (67, 5404948), (69, 19280015), (71, 32052856), (73, 43669891),
        (75, 54085999), (77, 63266253), (79, 71187816), (81, 77842572),
        (83, 83241098), (85, 87418974), (87, 90447291), (89, 92451149),
        (91, 93644201), (93, 94397935), (95, 95393499), (97, 97994819),
        (99, 105321028), (101, 126126436), (103, 195367048), (105, 626725133),
    ),
    ('direct', -1): (
        (1, 133923532), (3, 128590850), (5, 117952762), (7, 102063682),
        (9, 81004889), (11, 54884111), (13, 23834970), (74, 28056876),
        (76, 59574574), (78, 85921439), (80, 106883173), (82, 122248810),
        (84, 131801283), (86, 135301783), (88, 132463152), (90, 122902753),
        (92, 106054221), (94, 80989721), (96, 46026143),
    ),
    ('diff', 1): (
        (64, 4347572), (66, 39738099), (68, 69527815), (70, 94220429),
        (72, 114237480), (74, 129931049), (76, 141591316), (78, 149451460),
        (80, 153689810), (82, 154428724), (84, 151729049), (86, 145577853),
        (88, 135864857), (90, 122338265), (92, 104519747), (94, 81530973),
        (96, 51706820), (98, 11617381),
    ),
    ('diff', -1): (
        (57, 15580460), (59, 30140590), (61, 42763672), (63, 53672573),
        (65, 63056747), (67, 71078320), (69, 77877107), (71, 83574785),
        (73, 88278473), (75, 92083952), (77, 95078731), (79, 97345257),
        (81, 98964621), (83, 100021356), (85, 100610340), (87, 100847656),
        (89, 100889146), (91, 100964579), (93, 101445806), (95, 102995961),
        (97, 106936504), (99, 116304388), (101, 139686428), (103, 211600529),
        (105, 643517342),
    ),
}

FLOOR_SHIFTS = {
    ('direct', 1): (4633511331, (
        (0, 287779962), (2, 569864181), (4, 552804690), (6, 524468597),
        (8, 484999972), (10, 434599487), (12, 373523385), (14, 302082172),
        (16, 220639075), (18, 129608899), (20, 29451897), (69, 66735650),
        (71, 162999094), (73, 249978932), (75, 327252245), (77, 394453093),
        (79, 451277685), (81, 497490487), (83, 532933165), (85, 557537971),
        (87, 571349286), (89, 574560725), (91, 567583543), (93, 551182271),
        (95, 526768059), (97, 497106867), (99, 468305402), (101, 456751206),
        (103, 524100999), (105, 1173897519),
    )),
    ('direct', -1): (4670709047, (
        (1, 345914990), (3, 334410291), (5, 311459675), (7, 277180410),
        (9, 231747652), (11, 175393538), (13, 108405924), (15, 31125953),
        (72, 28515739), (74, 107990563), (76, 176915816), (78, 234838069),
        (80, 281331100), (82, 315988519), (84, 338409634), (86, 348175740),
        (88, 344809766), (90, 327705099), (92, 295993275), (94, 248280003),
        (96, 182066905), (98, 92314399),
    )),
    ('diff', 1): (4668763233, (
        (64, 81875933), (66, 159683454), (68, 225225731), (70, 279617017),
        (72, 323793091), (74, 358541944), (76, 384520818), (78, 402268116),
        (80, 412210394), (82, 414663666), (84, 409827329), (86, 397767209),
        (88, 378380821), (90, 351330859), (92, 315916754), (94, 270813947),
        (96, 213498444), (98, 138812823), (100, 34682883),
    )),
    ('diff', -1): (4630745553, (
        (59, 101958087), (61, 207603889), (63, 298623711), (65, 376595274),
        (67, 442857960), (69, 498555159), (71, 544667276), (73, 582038663),
        (75, 611400064), (77, 633387798), (79, 648560765), (81, 657416450),
        (83, 660407428), (85, 657960562), (87, 650502741), (89, 638500419),
        (91, 622528201), (93, 603401192), (95, 582458700), (97, 562249393),
        (99, 548461283), (101, 556708783), (103, 646992061), (105, 1317361355),
    )),
}


def verify(root):
    geometry = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 224
    unit = geometry["unit"]((108469,109))
    theta = unit.arg()
    direct = [[(q*j*theta).cos() for j in range(107)] for q in range(1,107)]
    difference = [[b-a for a,b in zip(old,new)] for old,new in zip(direct,direct[1:])]
    columns = {"direct": direct[:-1], "diff": difference}
    sizes = [1]+[2]*106
    reference = [arb(n)/213 for n in sizes]
    threshold = arb(9)/125
    slope = (threshold/(1-threshold)).sqrt()
    for name, rows in columns.items():
        A = arb_mat(rows)
        last = rows[-1]
        equal_last = sum(a*b for a,b in zip(last, reference))
        zero_distance = abs(equal_last)/(max(v.upper() for v in last)-min(v.lower() for v in last))
        assert zero_distance > arb(".00186" if name == "direct" else ".00200")
        tv_lowers, tv_uppers, floor_lowers, floor_uppers = [], [], [], []
        for sign in (1,-1):
            pair = (name, sign)
            change = [arb(0)]*107
            for j, amount in TV_SHIFTS[pair]:
                change[j] = arb(amount)/10**12
            change[-1] = -sum(change[:-1])
            base, extras = FLOOR_SHIFTS[pair]
            mass = [arb(n*base) for n in sizes]
            for j, amount in extras:
                mass[j] += amount
            total = sum(mass)
            floor_change = [v/total-u for v,u in zip(mass,reference)]
            cuts = []
            for mode, delta in (("TV",change),("floor",floor_change)):
                factor = arb("1.0000002") if mode == "TV" else arb("1.00001")
                weights = [u+factor*d for u,d in zip(reference,delta)]
                word_floor = min(w/n for w,n in zip(weights,sizes))
                assert word_floor > 0 and (sum(weights)-1).contains(0)
                moments = (A*arb_mat([[w] for w in weights])).entries()
                prior = sum(v*v for v in moments[:-1]).sqrt()
                assert prior > 0 and sign*moments[-1] > slope*prior
                # 唯一用作对偶的方向是见证在前104维的精确单位向量。
                score = [(sign*rows[-1][j]-slope*sum(v*row[j] for v,row in
                          zip(moments[:-1],rows[:-1]))/prior) for j in range(107)]
                mx, mn = max(v.upper() for v in score), min(v.lower() for v in score)
                at_uniform = sum(v*u for v,u in zip(score,reference))
                assert at_uniform < 0 and mx > 0
                if mode == "TV":
                    lower = -at_uniform/(mx-mn)
                    tv = sum(abs(factor*d) for d in delta)/2
                    tv_lowers.append(lower)
                    tv_uppers.append(tv)
                    cuts += [lower,tv,word_floor]
                else:
                    upper = mx/(213*(mx-at_uniform))
                    floor_lowers.append(word_floor)
                    floor_uppers.append(upper)
                    cuts += [word_floor,upper]
            print(name,sign,"总变差下上界、见证逐字最小质量、最优逐字质量下上界",
                  *(x.str(14) for x in cuts))
        if name == "direct":
            bounds = (".0018004786",".0018004792",".004670708",".004670710")
        else:
            bounds = (".001856046",".001856050",".004668762",".004668765")
        t0,t1,f0,f1=map(arb,bounds)
        assert min(tv_lowers) > t0 and min(tv_uppers) < t1
        assert max(floor_lowers) > f0 and max(floor_uppers) < f1


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
