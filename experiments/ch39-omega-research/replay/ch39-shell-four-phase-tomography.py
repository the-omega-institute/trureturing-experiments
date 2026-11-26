"""Four rational clock configurations recover a common arbitrary complex edge operator."""
from pathlib import Path
from runpy import run_path

from flint import nmod_mat

previous=run_path(str(Path(__file__).with_name("ch39-shell-nine-phase-probes.py")))


def main():
    edges=previous["orientation"]()
    prime=65537
    assert all(prime%k for k in range(2,257))
    i_mod=256
    assert (1+i_mod**2)%prime==0
    phase=(3+4*i_mod)*pow(5,-1,prime)%prime
    codes=[106*j+j*j%53 for j in range(1,49)]
    # Complexify the Hermitian observation map using all ordered entries.
    # Each e,f contributes at the four ordered pairs of endpoints.
    size=len(edges)**2
    A=nmod_mat(size,size,prime)
    for k in range(1,5):
        wave=[pow(phase,k*code,prime) for code in codes]
        offset=(k-1)*24**2
        for e,(u,v) in enumerate(edges):
            ref=pow(wave[e],-1,prime)
            for f,(a,b) in enumerate(edges):
                col=e*48+f
                A[offset+24*u+a,col]=1
                A[offset+24*u+b,col]=-wave[f]
                A[offset+24*v+a,col]=-ref
                A[offset+24*v+b,col]=ref*wave[f]
    assert A.rank()==size==2304
    assert 3*24**2<48**2 and 23**2+3*24**2<48**2
    print("Exact reconstruction with k=1,2,3,4; Gaussian-rational tomography determinant nonzero modulo 65537.")
    print("Minimum four arbitrary backgrounds; minimum five when the flat background is required.")


if __name__=="__main__":
    main()
