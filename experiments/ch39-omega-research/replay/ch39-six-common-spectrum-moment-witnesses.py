"""同一补空间固定谱正矩阵通过六壳十二个外包矩条件的四组见证。"""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, acb_mat


WITNESSES = (
    ((1, 1, 1, 1, 16, 16),
     (-4511, -53432, -6890, -13891, -19619, 259, 8684, 38821, 9912, 10024, 8917, -1230),
     (-19717, -22567, -32480, -10797, 31997, 1940, -27805, 8520, -19052, -15977, -6827, 5417)),
    ((1, 16, 16, 16, 16, 16),
     (5735, 47219, -16359, 6651, -327, 412, 22242, -21538, 4360, -12805, -12163, -380),
     (-4079, 49353, 21576, 811, 10639, 17551, 14352, 45446, 17790, -7923, -9152, -7848)),
    ((16, 16, 1, 1, 1, 1),
     (9762, -31197, -23363, 9443, 7153, -14950, 26330, 48808, 2998, 1566, -1632, -3593),
     (5836, 58171, 4256, -134, -9907, -2981, -32719, 133, -15077, 4463, 7013, 3036)),
    ((16, 16, 1, 1, 1, 16),
     (-126, -21516, 7138, 5794, -29224, 6533, -11841, 18766, -4446, -26922, 6363, -3808),
     (-9901, -19465, 19944, 606, 2684, 10658, -43260, 31925, 36734, -8705, -43858, -5355)),
)

def certify(root, levels, blocks, model, frames):
    sieve = run_path(str(root/"ch39-six-cyclic-line-sieve.py"))
    trivial = sieve["squared"](blocks[0])
    trace = sum(sieve["squared"](m) for m in blocks[1:])
    gram = sum((sieve["adj"](m)*m for m in blocks[1:]), acb_mat(2,2))
    lam = (1+(1-4*gram.det().real/trace**2).sqrt())/2
    assert arb(".99828") < lam < arb(".99830")
    for (labels, real, imag), floor in zip(WITNESSES, (".004181",".001024",".004160",".010093"), strict=True):
        q = acb_mat([[acb(x,y)] for x,y in zip(real,imag,strict=True)])
        norm2 = sum(x*x+y*y for x,y in zip(real,imag,strict=True))
        t = arb(sum(x*x+y*y for x,y in zip(real[:6],imag[:6],strict=True)))/norm2
        assert 1-lam < t < lam
        intervals = sieve["bounds"](model,levels,trivial,trace,[i-1 for i in labels])
        gap = arb(10)
        for shell, rows, (u,v) in zip(model,frames,intervals,strict=True):
            other = [i for i in range(19) if i not in shell["triv"]]
            amp = acb_mat([rows[i] for i in other])*q
            prob = [abs(x)**2/norm2 for x in amp.entries()]
            weights = [-shell["w"][i] for i in other]
            u += arb(0,(arb(".04")/(shell["ck"]*trace**2)).upper())
            v += arb(0,(arb(".045")/(shell["cw"]*trace**4)).upper())
            for power, box in ((2,u),(4,v)):
                moment = sum(w*p**power for w,p in zip(weights,prob,strict=True))
                margin = min((moment-box.lower()).lower(),(box.upper()-moment).lower())
                gap = min(gap,margin)
        assert gap > arb(floor)
        print(labels, t.str(10), gap.str(12), flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment = run_path(str(root/"ch39-six-trivial-shared-endpoints.py"))["experiment"]
    experiment(root,*source["verify"](root,return_operator="spectral"),check=certify)
