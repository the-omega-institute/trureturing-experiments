"""Choose one additional offset to distinguish the two fixed integer weight tables."""
import json
from fractions import Fraction as Q
from pathlib import Path


def main():
    root=Path(__file__).resolve().parent
    data=json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())
    den=1<<data["depth"]
    cases=[data["cases"][tick]["rows"] for tick in ("2","8")]
    # Floors and ceilings on each input height allow every actual ordinate in the cell.
    def bounds(rows,j):
        t=(j*den)**2
        lo=sum(w*h*h//(h*h+t) for _,h,w in rows)
        hi=sum(-((-w*(h+1)**2)//((h+1)**2+t)) for _,h,w in rows)
        return lo,hi
    at=[bounds(rows,2651) for rows in cases]
    gap=at[1][0]-at[0][1]
    assert gap > Q(".0063270893")*den
    last=9000
    # Positive averages obey F(j) <= sum(w*h^2)/j^2. The absolute
    # difference between two nonnegative numbers is at most their maximum.
    assert max(sum(w*(h+1)**2 for _,h,w in rows) for rows in cases) < gap*den**2*last**2
    for j in range(last):
        if j==2651:
            continue
        a,b=[bounds(rows,j) for rows in cases]
        assert max(b[1]-a[0],a[1]-b[0])<gap
    print("Unique maximum of |F2-F8| over nonnegative integer offsets at 2651.")
    print("Fixed witness values",*[float(Q(a+b,2*den)) for a,b in at])


if __name__=="__main__":
    main()

