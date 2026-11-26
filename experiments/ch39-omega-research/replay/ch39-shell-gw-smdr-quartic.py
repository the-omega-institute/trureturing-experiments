"""Validate the low-energy quartic ordering for the pinned SMDR beta functions.

The bundled upstream C and license are from davidgrobertson/SMDR at
5c4dd92c918cadb7f91ad3795694f8398b5c287a. Decimal coefficients and the
Zeta3 macro are interpreted as the indicated exact decimals. The independent
upstream C application was also run on the same MZ input points.
"""
import ast
from functools import lru_cache
import hashlib
from pathlib import Path
import re

from flint import arb, arb_series, ctx

ctx.prec = 144
folder = Path(__file__).resolve().parent
names = ("g3","g","gp","yt","yb","ytau","k","yc","ys","yu","yd","ymu","ye")
loop = 16*arb.pi()**2


@lru_cache(None)
def evaluator(order):
    raw = (folder/"ch39-smdr"/"betas.c").read_bytes()
    assert hashlib.sha256(raw).hexdigest()=="a066fecd4a5e5bbbd53a8d60737d5e01379266e1ae148b31b666d3527d5539f8"
    source = re.sub(r"/\*.*?\*/","",raw.decode(),flags=re.S)
    pattern = r"\b(Y2S|Y4S|chi4S|beta(?:"+"|".join(names)+r")[1-5])\s*=\s*([^;]+);"
    assignments = re.findall(pattern,source)
    scope = {"Zeta3":arb("1.2020569031595942854"),"__builtins__":{}}
    coefficients = {}
    def literal(m):
        token = m.group()
        if token not in coefficients:
            name = "c"+str(len(coefficients))
            coefficients[token] = name
            scope[name] = arb(token)
        return coefficients[token]
    lines = ["def beta(state):","    "+",".join(names)+" = state"]
    for v,exponents in {"g3":(2,4,6,8),"g":(2,4,6,8),"gp":(2,4,6,8),
                        "yt":(2,4,6,8),"yb":(2,4,6,8),"ytau":(2,4,6,8),"k":(2,3,4)}.items():
        for exponent in exponents:
            lines.append(f"    {v}{exponent} = {v}**{exponent}")
    terms = {n:[] for n in names}
    for symbol,body in assignments:
        if symbol.startswith("beta"):
            level = int(symbol[-1])
            if level>order:
                continue
            terms[symbol[4:-1]].append(f"factor{level}*{symbol}")
            scope["factor"+str(level)] = loop**(1-level)
        body = " ".join(body.split())
        expr = ast.parse(body, mode="eval")
        assert all(isinstance(n,(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Load,ast.Name,ast.Constant,
                                  ast.Add,ast.Sub,ast.Mult,ast.Div,ast.UAdd,ast.USub)) for n in ast.walk(expr))
        body = re.sub(r"(?<![.\w])(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?",literal,body)
        lines.append(f"    {symbol} = {body}")
    assert sorted(len(v) for v in terms.values())==sorted(min(order,n) for n in (5,4,4,5,3,3,4,3,3,3,3,2,2))
    lines.append("    return ["+",".join("+".join(terms[n]) for n in names)+"]")
    exec("\n".join(lines),scope)
    return scope["beta"]


def init(top,quartic):
    return [arb(x) for x in ("1.2123",".65096")]+[arb(".461228")*(arb(3)/5).sqrt()]+[
        arb(x) for x in (top,".01630",".0099378")]+[arb(quartic)/4]+[
        arb(x) for x in (".00356",".000293",".00000704",".0000154",".0005850",".00000277713")]


def coefficients(beta,y,order):
    series = [[v] for v in y]
    for k in range(order):
        if k==0:
            vals = beta(y)
        else:
            args = [arb_series(v,prec=k+1) for v in series]
            vals = [p[k]/(k+1) for p in beta(args)]
        for cs,v in zip(series,vals):
            cs.append(v)
    return series


def advance(beta,y,h):
    pad = arb(0,".022")
    tube = [v+pad for v in y]
    bounds = beta(tube)
    assert all(abs(h*z)<arb(".0219") for z in bounds)
    tube = [v+arb(0).union(h)*z for v,z in zip(y,bounds)]
    assert beta(tube)[6]<0
    # Evaluate the cubic Taylor polynomial and bound its fourth coefficient
    # everywhere in the verified trajectory tube.
    initial = coefficients(beta,y,3)
    last = coefficients(beta,tube,4)
    result = []
    for row,tail in zip(initial,last):
        v = tail[4]
        for a in reversed(row):
            v = a+h*v
        result.append(v)
    return result


def trajectory(top,quartic,order,scales):
    beta = evaluator(order)
    y = init(top,quartic)
    t = arb(0)
    out = []
    h = arb(1)/1536
    for mu in scales:
        target = (arb(mu)/arb("91.1876")).log()/loop
        steps = int(float(((target-t)/h).lower()))
        for _ in range(steps):
            y = advance(beta,y,h)
        remainder = target-t-steps*h
        assert 0<remainder<h
        y = advance(beta,y,remainder)
        out.append(y[6])
        t = target
    return out


if __name__=="__main__":
    ctx.cap = 5
    central = trajectory(".967",".55853",5,("1000","7e10","7.4e10","1.44e11","1.48e11","1.70e11","1.80e11"))
    print("SMDR central quartics",[z.str(9) for z in central])
    assert central[3]>0 and central[4]<0
    for top,k,expected in ((".963",".56010",1),(".971",".55696",-1)):
        ends = trajectory(top,k,5,("7e10","1.80e11"))
        assert all(expected*z>0 for z in ends)
        print("top, four-times quartic, below and above both shell thresholds",top,k,[z.str(9) for z in ends])
