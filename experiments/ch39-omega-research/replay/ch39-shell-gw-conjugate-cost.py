"""Bound the norm of a single conjugate Higgs profile realizing the quark data."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

c = run_path(str(Path(__file__).with_name("ch39-shell-gw-conjugate-universality.py")))
m = c["c"]["a"]
F = c["output_basis"]()
adj, hs = m["adj"],m["hs"]
rows, count = c["rows"],c["count"]
G = adj(F)*F
up = [arb(t) for t in (".00000615",".00311",".8616")]
down = [arb(t) for t in (".0000135",".000268",".01401")]
u2,d2 = (sum(t**2 for t in sector) for sector in (up,down))


def fit_cost(units,scale_u=arb(1),scale_d=arb(1)):
    U,A,B = [m["four"]["exact_frame"](np.array(z).reshape(3,3)) for z in units]
    ck,du,dd = m["orb"]["reference"]()
    u = U*du*A*up[-1]/scale_u
    d = U*ck*dd*B*down[-1]/scale_d
    target = acb_mat([[z] for z in u.entries()+d.conjugate().entries()])
    paired = acb_mat(rows(F)+rows(F.conjugate()))
    v = adj(paired)*(paired*adj(paired)).solve(target)
    return hs(v,v).real.sqrt()


# Three orbit representatives used for the upper bound.
units = [
    [
        (-0.7242346374056949+0.4186116597868701j), (-0.42785485455237354+0.2909247071004517j), (-0.15243543939598822+0.09651395403233244j),
        (0.5193114878685622+0.10066579789290536j), (-0.8416591131454901-0.10558232578769518j), (-0.022187792937301942+0.012328312427774327j),
        (0.11182349363903417-0.08902780602166945j), (0.09529527103217857-0.06069727583490398j), (-0.7802207471423737+0.5983809771199677j),
    ],
    [
        (-0.3859907468126761-0.7682102103153244j), (0.19192658748805916+0.025055824129104096j), (0.42523505824786684-0.20633892524576958j),
        (-0.4125495183432653-0.21572960408801697j), (0.06193559885915357-0.059893979057427015j), (-0.4104703139903662+0.7793294855412763j),
        (0.1971265403843873+0.07258771282453605j), (0.7625891347196707-0.6112044626526425j), (0.013536418088186685+0.023996517441216107j),
    ],
    [
        (-0.26977275187095506+0.022990914531496147j), (0.9459899607782739+0.1597385342522211j), (0.01632933162829668-0.07755016381150752j),
        (0.07245934860991592-0.8979053659671852j), (0.11064507736945561-0.23156098092043076j), (-0.11883161929376189+0.32944168192783696j),
        (0.11696366406708784+0.31863643553636717j), (0.08516439539268704+0.08044127987614078j), (0.14307839397710095+0.9222771398353791j),
    ],
]


split_units = [
    [
        (-0.12091850536588458+0.042996546910222734j), (0.01618872729070074+0.05810297534392888j), (-0.9658824651163185-0.21670958668136983j),
        (-0.08989846257425499-0.4246531610861978j), (-0.88024698582933+0.15649839635816926j), (-0.03644816350034939+0.10456092223947738j),
        (-0.7838481143540957-0.42508834945660423j), (0.22452043250668036-0.38291779878681875j), (0.042497512763456015+0.07772081412360492j),
    ],
    [
        (0.11023849090126693+0.05609601352464558j), (0.08330810068743368-0.32172219892603315j), (0.9011162209756685-0.24948918967174463j),
        (0.05727848727196709-0.07837973081760453j), (-0.9094530864188013-0.21161593663408673j), (0.09785429258097653+0.3303242484015327j),
        (-0.653890124546564+0.7400704088826275j), (-0.1234965236126308-0.04916273689928055j), (0.004567204198983624-0.08387055215464151j),
    ],
    [
        (-0.00016737067064144078-0.08432611088880579j), (-0.22854305584596568+0.4985347146901123j), (-0.048075636476268224-0.8305474231254297j),
        (-0.017328311862952438-0.11256220762954747j), (-0.8092807966647747-0.18022283904650968j), (-0.5245771323848927+0.15630940994067422j),
        (-0.6910281696555831-0.7088008371703897j), (0.07220852584717083-0.08124936371450579j), (0.09089930179355678-0.0017968287172219638j),
    ],
]


if __name__ == "__main__":
    lam = arb(".000131336")
    assert count(rows(G-62*G.conjugate()),lam)==54
    lower = ((u2-62*d2)/lam).sqrt()
    upper = fit_cost(units)
    assert lower > arb("74.56") and upper < arb("76.25")
    def pencil_bound(t,limit):
        t,limit = arb(t),arb(limit)
        assert count(rows(G-t*G.conjugate()),limit) == 54
        return t,limit
    # Apply the same pencil after conjugating the scalar coefficients to reverse the roles.
    t,limit = pencil_bound("1.235",".000787")
    ratio_lo = ((u2-(4*arb.pi())**2*limit)/(d2*t)).sqrt()
    t,limit = pencil_bound("1.2",".000921")
    ratio_hi = ((t*u2+(4*arb.pi())**2*limit)/d2).sqrt()
    assert ratio_lo > arb("50.48") and ratio_hi < arb("72.65")
    split = fit_cost(split_units,arb(1),1/arb("61.5"))
    assert split < arb("2.437")
    print("ratio 61.5 effective up strength",split.str(9),"down",(split/arb("61.5")).str(9))
    print("equal-coupling norm cost between",lower.lower().str(9),upper.upper().str(9),
          "coupling ratio for up strength <= 4π in",ratio_lo.lower().str(8),ratio_hi.upper().str(8))
