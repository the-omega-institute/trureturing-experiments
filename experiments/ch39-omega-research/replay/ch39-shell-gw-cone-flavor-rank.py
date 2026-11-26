"""The nine light fields transmit a full flavor matrix; the ground field alone approaches rank one."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

cone = run_path(str(Path(__file__).with_name("ch39-shell-gw-adjoint-cone-thresholds.py")))
first, shell, full = (cone[key] for key in ("first","shell","full"))
q, eye = arb,shell["yuk"]["identity"]
w = q("1e-18")


def normalized_trial(K):
    guess = np.linalg.eigh(full["model"]["array"](K))[1][:,0]
    v = acb_mat([[acb(float(a.real),float(a.imag))] for a in guess])
    return v/shell["hs"](v,v).real.sqrt()


def fermion(sign):
    kernel = shell["yuk"]["model"]["action_kernel"]
    full_center = acb_mat(kernel(sign,q(0),q(0))[1])
    box = acb_mat(kernel(sign,q(0,w),q(0,w))[1])
    delta = sum(abs(z).abs_upper()**2 for z in (box-full_center).entries()).sqrt()
    assert delta < q("2.89e-18")
    # At the center the full kernel is K tensor I_3.
    center = acb_mat([[full_center[3*i,3*j] for j in range(12)] for i in range(12)])
    scale,gap = q(8),q(".33")
    for cut in (-gap,gap):
        assert full["count"](full["rows"](center),cut) == 6-sign
    assert max(sum(abs(z) for z in row).upper() for row in full["rows"](center)) < scale
    x,b = center/scale,gap/scale
    unit = eye(12)
    while 1-b > q("1e-24"):
        x = x*(3*unit-x*x)/2
        b = b*(3-b*b)/2
    offset = 0 if sign>0 else 6
    block = acb_mat([[(unit[i,j]-sign*x[i,j])/2 for j in range(offset,offset+6)]
                     for i in range(offset,offset+6)])
    gap0 = q(".005")
    assert full["count"](full["rows"](block),gap0+(1-b)/2) == 1
    v = normalized_trial(block)
    rem = block*v
    norm = shell["hs"](rem,rem).real.sqrt()+(1-b)/2
    # Index and opposite-chirality exclusion give the exact zero.
    assert 2*norm/gap0 < q("1e-11")
    move = delta/(2*(gap-delta))
    # Each unit vector in the new kernel has center-complement component <= move/gap0.
    return v,2*move/gap0


def main():
    single,_,phases,_ = shell["transported_spectra"]()
    v = normalized_trial(single)
    e = shell["hs"](v,single*v).real
    rem = single*v-e*v
    assert e < q(".014") and full["count"](full["rows"](single),q("3.7"))==1
    assert 2*shell["hs"](rem,rem).real.sqrt()/(q("3.7")-e) < q("1e-11")
    l,err_l = fermion(1)
    r,err_r = fermion(-1)
    eta = abs(sum(phases[i]*v[i,0]*l[i,0]*r[i,0] for i in range(6)))
    assert eta-q("3e-11") > q(".40812589") and eta+q("3e-11") < q(".40812591")
    s = cone["velocity"]()
    # Complement at distance >3 and nine-band effective error <= (8+16/3)t².
    D = 8+q(16)/3
    isolation = s/2-D*w
    field_error = 2*(4*w/3+D*w/isolation)
    matrix_error = err_l+err_r+field_error
    assert matrix_error < q("5.5e-15")
    ratio = matrix_error/(q(".40812589")-matrix_error)
    assert ratio < q("1.35e-14")
    # Davis-Kahan projection across the 3.69 gap, using its midpoint as the sign cut.
    band_move = 4*w/(2*(q("3.69")/2-4*w))
    assert err_l+err_r+2*band_move < q("5e-15")
    mixing = q(".00311")/(4*q.pi())-matrix_error
    assert mixing > q(".00024748")
    print("Ground second/first singular ratio <",ratio.upper().str(9),
          "; nine light singular values in (.40812588,.40812592);",
          "minimum outside-ground norm at strength <= 4pi",mixing.lower().str(9))


if __name__ == "__main__":
    main()
