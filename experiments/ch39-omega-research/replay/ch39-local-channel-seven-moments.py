"""Count channels from consecutive exact moments for two identical independent shells."""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

local = run_path(str(Path(__file__).with_name("ch39-local-product-contact.py")))


def convolve(left,right,k):
    return sum(comb(k,j)*left[j]*right[k-j] for j in range(k+1))


def verify():
    laplacian,port = local["lifted_laplacian"]()
    single = [(laplacian**j)[port,port] for j in range(8)]
    pair = [convolve(single,single,j) for j in range(8)]
    bare = [comb(2*j,j)+int(j==7) for j in range(8)]
    joint = [convolve(bare,pair,j) for j in range(8)]
    recovered_pair, recovered = [1], [fmpq(1)]
    for j in range(1,8):
        recovered_pair.append(joint[j]-sum(comb(j,i)*bare[i]*recovered_pair[j-i]
                                           for i in range(1,j+1)))
        known = sum(comb(j,i)*recovered[i]*recovered[j-i] for i in range(1,j))
        recovered.append((recovered_pair[j]-fmpq(known))/2)
    assert recovered == single
    polynomial = [-175,199,-29,1]
    assert sum(a*b*single[1+i+j] for i,a in enumerate(polynomial)
               for j,b in enumerate(polynomial))==0
    zero_mass = sum(a*single[i] for i,a in enumerate(polynomial))/fmpq(polynomial[0])
    assert zero_mass == fmpq(1,25)

    # Positive part of a competing 26-state port, before lifting its soft mode.
    moments = [fmpq(25,26),*single[1:]]
    gram = fmpq_mat([[moments[i+j] for j in range(4)] for i in range(4)])
    shifted = fmpq_mat([[moments[i+j+1]+int(i==j==3) for j in range(4)]
                        for i in range(4)])
    for matrix in (gram,shifted):
        for size in range(1,5):
            assert fmpq_mat([[matrix[i,j] for j in range(size)] for i in range(size)]).det()>0
    operator = gram.inv()*shifted
    positive = [(gram*operator**k)[0,0] for k in range(8)]
    assert positive[:7] == moments[:7] and positive[7] == moments[7]+1
    print("clock-site four moments 1 through 7:",*[int(x) for x in joint[1:]])
    print("recovered single-shell zero weight",zero_mass,
          "-> total channels 625; 676-state competitor matches orders 1..6.")


if __name__=="__main__":
    verify()
