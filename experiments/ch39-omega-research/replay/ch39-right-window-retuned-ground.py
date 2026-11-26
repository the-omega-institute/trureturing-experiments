"""Recover the zero-index ground for each top group after the W_R window is included."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
window = run_path(str(root/"ch39-first-charged-vector-window.py"))
study, q, retuned = window["study"], window["q"], window["retuned"]
adj, empirical = window["adj"], study["empirical"]
log_sets = (("24.32546407460937","26.71295730013009","26.81727886982861"),
            ("24.32519852473211","26.71283158255943","26.81727833146807"),
            ("24.32404818846499","26.71228692560664","26.81727599938087"))


def predictions(logs, external_time):
    h3,h8,cal = map(q,logs)
    g,u,w,_,_,active = retuned["current"]["neutral"]["parameters"](logs)
    L = (2*(cal-external_time)).exp()
    c,z,rho = 18*u**2,1/(2*w**2*g**2),1/g**2
    masses = (h3-cal).exp(),(h8-cal).exp()
    nodes = tuple(map(q,("1/1024","1/50","1/5","1","1/2")))
    readings = [1+(s/(L*g**2)+4/(1/(2*w*w)+active(-s/L)))/c for s in nodes]
    widths = run_path(str(root/"ch39-current-width-normalization.py"))
    currents = retuned["current"]["light_current"](logs)[2]
    return (L,z,rho,rho,c,masses),readings,widths["reference_width"](currents)


def verify():
    word = empirical["endpoint"]
    bands = study["tables"]["table_bands"]()
    nominal = study["nominal"]
    gap = window["graph_gap"]()
    _,reference,reference_width = predictions(retuned["logs"],nominal[2])
    for index,logs in enumerate(log_sets):
        h3,h8,cal = map(q,logs)
        model,readings,width = predictions(logs,nominal[2])
        onset,_,end,_,first,second = window["window_data"](model,gap)
        inputs = []
        for upper in (False,True):
            inputs.append([q(abs(Fraction(w,word[0]))) for w in word[1:]]+
                          [q(row[index if i==2 else (-1 if upper else 0)][0]).sqrt()
                           for i,row in enumerate(bands)])
        sm,trip,both = adj["flow"](False,False),adj["flow"](True,False),adj["flow"](True,True)
        first_flow = {**sm,"NEUTRAL_YUKAWA_BOUNDS":(first[3],first[2])}
        second_flow = {**trip,"NEUTRAL_YUKAWA_BOUNDS":(second[3],second[2])}
        segments = ((q(0),onset,160,sm),(onset,h3,8,first_flow),(h3,end,10,second_flow),
                    (end,h8,30,trip),(h8,cal,10,both),(cal,q(33),64,both))
        offsets = [a+b for a,b in zip(first[0],second[0],strict=True)]
        rates = [max(a,b) for a,b in zip(first[1],second[1],strict=True)]
        empirical["triplet"]["verify_ground"](inputs[0],segments,4,cal,q("1.5e-6"),q(".136"),
            upper_values=inputs[1],logarithm=(onset,end,offsets),stationary_readout=True,
            maximum_log_rates=rates)
        assert all(v/r < q(".997") for v,r in zip(readings,reference))
        assert width/reference_width < q(".999")
        print("top group",index,"mass and scale changes (%)",
              *((100*((v-r).exp()-1)).str(10) for v,r in zip((h3,h8,cal),nominal)))
        print("ratio changes (%)",*( (100*(v/r-1)).str(9) for v,r in zip(readings,reference)))
        print("lepton width change (%)",(100*(width/reference_width-1)).str(10))


if __name__ == "__main__":
    verify()
