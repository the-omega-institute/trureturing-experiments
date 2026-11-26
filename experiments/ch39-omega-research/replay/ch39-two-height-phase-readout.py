"""Bottleneck and robust quarter-turn resolution for the 98-probe two-height design."""
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
model = run_path(str(root/"ch39-clock-schur-mediators.py"))
B, Q, root_cell = [model[k] for k in ("B","Q","root_cell")]


def main():
    ids, z, p, remainder, k, gamma, theta = model["calibration"]()
    b = [B(1),B(1)]+[B(k)*root_cell(x/remainder) for x in p[2:]]
    star = {i: theta*root_cell(z[0]*v)*y for i,v,y in zip(ids[1:],z[1:],b[1:])}
    weakest = min(star,key=lambda i:star[i].lo)
    value = star[weakest]
    assert all(value.hi < s.lo for i,s in star.items() if i!=weakest)
    for j in range(2,len(ids)):
        for i in range(1,j):
            mag = theta*root_cell(z[i]*z[j])*b[i]*b[j]
            if i>=2:mag = mag*B(gamma)
            assert mag.hi < min(star[ids[i]].lo,star[ids[j]].lo)
    bits = 1
    while Q(2)**(1-bits) >= value.lo:
        bits += 1
    assert value.hi < Q(2)**(2-bits)
    assert weakest==2315 and bits==28
    assert Q("1.13001e-8")<value.lo<=value.hi<Q("1.13003e-8")
    # A native phase kick C^30 is read with the identity and C^n.
    word = model["base"]["word"]
    def cosine(n):
        a,b = word["power"](n)[0]
        return B(Q(a,word["DEN"]),Q(b,word["DEN"]))
    kick = B(1)-cosine(30)
    threshold = Q(1,2**54)
    least = None
    for n in range(31):
        c = cosine(n)
        absolute = c if c.lo >= 0 else B(-c.hi,-c.lo)
        assert absolute.lo >= 0
        square = value*value*kick*(B(1)-absolute)
        if square.lo > threshold:
            least = n
            break
        assert square.hi <= threshold
    assert least==19
    print("Native C^30 kick at 28 bits: least second probe exponent after identity:", least)
    print("50 sites, 98 real two-site measurements. Bottleneck edge",ids[0],weakest)
    print("Amplitude",float(value.lo),float(value.hi),"least uniform quarter-turn precision",bits)


if __name__=="__main__":
    main()
