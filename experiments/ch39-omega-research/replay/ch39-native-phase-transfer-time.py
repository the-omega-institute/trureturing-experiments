"""Charged native hopping time for a dyadic phase on the full matching band."""
from collections import Counter
from fractions import Fraction
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


def verify():
    ctx.prec = 192
    shell = run_path(str(Path(__file__).with_name("ch39-rotor-shell-returns.py")))
    _, adj, _ = shell["charged_shell"]()
    edges = [(i,j,a) for i,row in enumerate(adj) for j,a,q in row if q==1]
    assert len(edges)==6 and all(a==Fraction(27,52136) for _,_,a in edges)
    out = {i:j for i,j,_ in edges}
    sources = set(out)-set(out.values())
    chains = [(i,out[i],out[out[i]]) for i in sources]
    assert len(chains)==3 and len(set(v for chain in chains for v in chain))==9
    assert {e for chain in chains for e in zip(chain,chain[1:])} == set(out.items())

    a = arb(27)/52136
    speed = arb(2).sqrt()*a
    pulse = arb.pi()/speed
    # Three-vertex perfect transfer checked from exp(-i H pi/(sqrt(2) a))=I-H^2/a^2.
    w = acb(3,4)/5
    hop = acb_mat([[0,a/w,0],[a*w,0,a/w],[0,a*w,0]])
    endpoint = (acb_mat([[1,0,0],[0,1,0],[0,0,1]])-hop*hop/a**2)[2,0]
    assert abs(endpoint+w**2) < arb("1e-50")
    log_two, pi, eps = arb(2).log(), arb.pi(), arb("1e-6")

    def angle(g):
        return 2*(1/(2*g)).atan()

    def time_bound(lo,hi):
        angle_span = angle(lo)-angle(hi)
        diff = (hi-lo)*log_two
        assert diff < pi
        return (2*(diff/2).sin()-2*eps)/(angle_span*speed)

    continuous = time_bound(arb(11644)-arb(1)/32,arb(11644))
    z1,z2 = acb.zeta_zeros(12053,2)
    discrete = time_bound(z1.imag,z2.imag)
    assert continuous > 128303659431 and discrete > 123395966066

    # Finite Fourier clock: remove a known integer carrier before phase estimation.
    alpha_min, alpha_max = 2*angle(arb(11644)), 2*angle(arb(10570))
    width = arb("1e-6")
    high = 1/(2*((alpha_min-width)/4).tan())
    low = 1/(2*((alpha_max+width)/4).tan())
    carrier = 43000000
    derivative = max((high**2+arb(1)/4)*log_two/2-carrier,
                     carrier-(low**2+arb(1)/4)*log_two/2)
    second = log_two*high*(high**2+arb(1)/4)/2 + derivative**2
    local = (second+derivative*width/3)/(1-width**2/12)
    far = 2/(1-width.cos())+derivative/(width/2).tan()
    clock = 2**34
    error = max(local,far)*(1-(pi/(clock+1)).cos())
    total_time = (2*(clock-1)+carrier)*pulse
    assert error < arb("3.764e-7") and total_time < arb("1.476e14")
    # Direct finite-register check of the first Fourier moment at off-grid input.
    for m in (7,16):
        count, alpha = 2*m, arb("0.1728")
        amps = [(arb(2)/(m+1)).sqrt()*(pi*(j+1)/(m+1)).sin() for j in range(m)]
        first = acb(0)
        for k in range(count):
            delta = 2*pi*k/count-alpha
            response = sum(v*acb(0,-j*delta).exp() for j,v in enumerate(amps))
            first += abs(response)**2*acb(0,delta).exp()/count
        assert abs(first-(pi/(m+1)).cos()) < arb("1e-50")
    print("Continuous and two-zero necessary native times", continuous.str(14),discrete.str(14))
    print("Finite-clock amplitude error",error.str(9),"native hopping time",total_time.str(12))


if __name__ == "__main__":
    verify()
