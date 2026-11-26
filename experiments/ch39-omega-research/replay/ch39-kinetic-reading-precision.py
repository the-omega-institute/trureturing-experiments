"""Optimal pairs of finite gauge and kinetic measurement errors near the third pole."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
gauge=run_path(str(root/"ch39-kinetic-reading-gauge.py"))
cross,B=gauge["cross"],gauge["B"]
MARK=Q(15,8)
WINDOWS={
    3:(".00082895008463",".00082895008469"),
    4:(".00082775992148",".00082775992155"),
    5:(".00082716497872",".00082716497879"),
    6:(".00082686754202",".00082686754210"),
}


def generate(model):
    gen=run_path(str(root/"ch39-kinetic-reading-origin.py"))["candidate_rows"]
    bare_q,bare_t=model.q2,model.target
    wide,targets=gauge["input_boxes"](model)
    model.q2=wide
    chamber=gen(model,[(Q(1,1152),None,-1,620),(Q(1,2048),None,1,640)])
    rows={}
    for k,(lo,hi) in WINDOWS.items():
        model.q2,model.target=bare_q,bare_t
        model.q2,targets=gauge["input_boxes"](model,k)
        rows[str(k)]=gen(model,[(lo,None,-1),(hi,None,1)])
    return {"chamber":chamber,"windows":rows}


def main():
    model=cross["ReadingBounds"]()
    baseq,baset=model.q2,model.target
    records=cross["base"]["data"]("ch39-kinetic-reading-precision-roots.json")
    model.q2,targets=gauge["input_boxes"](model)
    assert [(Q(r["c"]),r["Y"],r["sign"]) for r in records["chamber"]]==[(Q(1,1152),620,-1),(Q(1,2048),640,1)]
    for row,target in zip(records["chamber"],targets,strict=True):
        model.target=target
        model.verify_sample(row)
    # Between these two Y bounds only the third pure pole can be crossed.
    max_higgs=Q(125**2*640,model.q2.lo*1152)
    for lam in model.modes:
        def secular(y,c):
            return B(c*y)-lam+lam*(B(1)-model.p)*model.response(y)
        assert secular(Q(1,10**5),Q(1,1152)).hi<0
        assert Q(125**2*640,model.q2.lo)<Q(1,10**5) and model.q2.lo*Q(1,10**5)/640>1000**2
        assert secular(Q(620),Q(1,2048)).lo>max_higgs
        assert secular(Q(640),Q(1,1152)).hi<0
    assert model.poles[1].hi<620<model.poles[2].lo<=model.poles[2].hi<640<model.poles[3].lo
    assert model.reading(Q(1,1152)).hi<MARK-Q(1,256) and model.reading(Q(1,2048)).lo>MARK+Q(1,256)
    for k,(a,b) in WINDOWS.items():
        model.q2,model.target=baseq,baset
        model.q2,opts=gauge["input_boxes"](model,k)
        model.target=opts[1]
        rows=records["windows"][str(k)]
        assert [(Q(r["c"]),r["mass"],r["sign"]) for r in rows]==[(Q(a),None,-1),(Q(b),None,1)]
        for row in rows:model.verify_sample(row)
        critical=B(model.reading(Q(b)).lo,model.reading(Q(a)).hi)
        if k==3:
            assert critical.hi<MARK
        else:
            bits={4:12,5:10,6:9}[k]
            assert MARK+Q(1,2**bits)<critical.lo<=critical.hi<MARK+Q(2,2**bits)
        print(k,"lowest possible critical value bracket",float(critical.lo),float(critical.hi))
    # Even perfectly specified gauges admit a population above threshold with 8-bit kinetic tolerance.
    assert MARK<Q("1.8780365378188")<Q("1.8780365378190")<MARK+Q(1,256)
    print("Frontier (additional gauge bits, kinetic error exponent): (4,12), (5,10), (6,9).")


if __name__=="__main__":
    import json,sys
    if "--write" in sys.argv:
        path=root/"ch39-kinetic-reading-precision-roots.json"
        path.write_text(json.dumps(generate(cross["ReadingBounds"]()))+"\n")
    main()
