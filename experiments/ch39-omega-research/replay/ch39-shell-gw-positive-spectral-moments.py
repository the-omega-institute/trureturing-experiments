"""Exclude simultaneous mass-norm fits from every positive fourth/sixth moment mixture."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb

folder=Path(__file__).resolve().parent
s=run_path(str(folder/"ch39-shell-gw-sextic-spectral-mode.py"))
kin=run_path(str(folder/"ch39-shell-gw-native-kinetic-cost.py"))
m,n=kin["m"],kin["n"]
hs,adj,count,rows=kin["hs"],kin["adj"],kin["count"],kin["rows"]
p=m["c"]["p"]
G=m["G"]
t=n["relative_phase"](arb(676))
f=kin["arc"]["conjugate_output"](t)
ratio=m["d2"]/m["u2"]
constraint=adj(f)*f-ratio*G


def site_frame():
    basis=[s["eye"](3)/arb(3).sqrt()]+s["w"]["gen"]
    return acb_mat([[p[i//9]*basis[j%9][(i%9)//3,i%3] if i//9==j//9 else 0
                     for j in range(54)] for i in range(54)])


if __name__=="__main__":
    W=site_frame()
    Q=adj(W)*s["quadratic"](s["a"],s["b"],6)*W
    K=kin["scalar_kinetic"]()
    assert count(rows(Q+100*constraint),arb("176"))==0
    assert count(rows(K+10*constraint),arb("1.09"))==0
    singlet=acb_mat([[1/arb(6).sqrt() if i%9==0 else 0] for i in range(54)])
    q0=hs(singlet,Q*singlet).real
    k0=hs(singlet,K*singlet).real
    assert q0 < arb("161.398") and k0 < arb(".01377")
    # At unit HS mass, write E=V²-I/18. If ||E||_HS>=.05, the onsite
    # fourth-power surplus is at least .03. Otherwise the anticommutator
    # coefficient is at least (2/3-12*.05)*kinetic.
    split=arb(".05")
    slack=min(12*split**2,(arb(2)/3-12*split)*arb("1.09"))
    assert slack-2*k0/3 > arb(".0208")
    print("for fixed mass: excess in the sixth moment exceeds",
          (176-q0).lower().str(8),"times mass + .0208 times mass squared;",
          "no upper bound on the common Yukawa strength used")
