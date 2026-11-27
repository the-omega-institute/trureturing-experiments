"""Static census of the arithmetic used by the ch39 replay scripts and of their dependence on Riemann zeros.

The scripts are parsed with ast and never executed.  The population is the frozen list of the
1706 ch39-*.py files tracked in replay/ at commit fe26883c21fcb459d80abf1d5fc14f8fe61f47e6, stored with their
SHA-256 in ch39-arithmetic-tier-census.json; scripts added later are outside the census.  Run with
--write (inside that checkout) to regenerate the stored lists.

Direct tier of one script (first rule that fires, weakest first):
  4 float     floating evidence that is not display-only and not steering-only:
              import of scipy/cmath/statistics/mpmath (decimal is not counted: the archive only
              parses and compares decimal strings with it, which is exact); math transcendental or float
              functions and constants (sqrt, exp, log*, trig, pi, e, inf, hypot, isclose, ...);
              numpy float operations (np.linalg, np.exp, np.pi, np.random, dtype=float/complex,
              astype(float/complex), ...); sympy N/evalf/Float/lambdify/nroots/nsolve;
              float or complex literals; float()/complex() calls; int-literal/int-literal
              division; ** with a float literal exponent; random()/uniform()/gauss() draws.
  3 steer     only floating evidence that sits inside key=..., int(...), round(...) or range(...)
              (search grids, orderings, rounding of candidates) or in a tol=... argument
              (target tolerance of a certified flint root finder).
  2 ball      a flint ball type is referenced (arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ...).
  1 symbolic  sympy functions that leave Q: sqrt, root, pi, E, exp, log, trig, RootOf, series, oo, ...
  0 rational  none of the above: integers, Fraction, flint fmpz/fmpq/fmpz_mat/fmpq_mat/...,
              sympy Rational/Matrix/Poly/symbols, Gaussian rationals via sympy I.
Display-only floating evidence (inside print(...), an f-string, .format/.write/.dump/.dumps,
or an assert message) and a numpy import without the listed float operations are recorded
but do not lower the tier.
Closure tier: the weakest direct tier over the run_path/string-reference closure, where a script
depends on every population script whose file name appears in one of its string constants (plus,
for scripts that build names as "ch39-"+name+".py", every constant name that completes to a
population script).
Zero dependence of the closure: H if some script in it calls zeta_zero or zeta_zeros (python-flint,
the only source of zero heights in the archive; mpmath is never imported); N if no heights but
zeta_nzeros or riemann_siegel_z is called; D if neither but a script in it reads a data file that
is also read by some script whose closure calls one of these functions; 0 otherwise.
"""
import ast
import hashlib
import json
import random
import re
import sys
from collections import Counter
from pathlib import Path

R = Path(__file__).resolve().parent
STORE = R / "ch39-arithmetic-tier-census.json"
COMMIT = "fe26883c21fcb459d80abf1d5fc14f8fe61f47e6"
SOURCES = {}


def git_blob(name):
    """The bytes of replay/<name> at COMMIT, or None when no checkout containing COMMIT is available."""
    import subprocess
    try:
        out = subprocess.run(["git", "show", "%s:./%s" % (COMMIT, name)], cwd=R, capture_output=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout


def load_sources(digests):
    """Read each census script as it was at COMMIT: from the working tree when its SHA-256 still
    matches, otherwise from the git object store of a checkout that contains COMMIT."""
    missing = []
    for k, h in digests.items():
        data = (R / k).read_bytes() if (R / k).exists() else b""
        if hashlib.sha256(data).hexdigest() != h:
            data = git_blob(k)
            if data is None or hashlib.sha256(data).hexdigest() != h:
                missing.append(k)
                continue
        SOURCES[k] = data.decode("utf-8")
    if missing:
        raise SystemExit("census: %d script(s) differ from the pinned commit %s and cannot be read from git "
                         "(first: %s). Run inside a checkout of the trureturing-experiments repository whose "
                         "history contains that commit." % (len(missing), COMMIT, missing[0]))


SELF = Path(__file__).name

BALL = {"arb", "acb", "arb_mat", "acb_mat", "arb_poly", "acb_poly", "arb_series", "acb_series"}
EXACT_FLINT = {"fmpz", "fmpq", "fmpz_mat", "fmpq_mat", "fmpz_poly", "fmpq_poly", "nmod", "nmod_mat", "nmod_poly", "fmpz_mod"}
FLOAT_MOD = {"scipy", "cmath", "statistics", "mpmath"}
MATH_FLOAT = {"sqrt", "exp", "log", "log2", "log10", "log1p", "expm1", "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
              "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "pi", "e", "tau", "hypot", "pow", "fsum", "erf", "erfc",
              "gamma", "lgamma", "fabs", "inf", "nan", "dist", "cbrt", "exp2", "isclose", "ldexp", "frexp", "fmod", "modf",
              "copysign", "radians", "degrees"}
SYM_IRR = {"sqrt", "root", "cbrt", "pi", "E", "exp", "log", "cos", "sin", "tan", "acos", "asin", "atan", "atan2", "cosh",
           "sinh", "tanh", "acosh", "asinh", "RootOf", "CRootOf", "rootof", "real_root", "expand_trig", "chebyshevt",
           "nsimplify", "oo", "integrate", "series", "zeta", "gamma"}
SYM_FLOAT = {"Float", "N", "evalf", "lambdify", "nroots", "nsolve"}
NP_FLOAT = {"linalg", "exp", "sqrt", "pi", "sin", "cos", "tan", "log", "log2", "log10", "random", "linspace", "polynomial",
            "inf", "nan", "mean", "std", "isfinite", "float64", "float32", "complex128", "complex64", "arctan2", "arccos",
            "arcsin", "arctan", "angle", "fft", "rint", "round", "floor", "ceil", "real", "imag", "conj", "sign", "hypot",
            "sinh", "cosh", "tanh", "e", "divide", "true_divide", "reciprocal", "power", "float_", "cdouble", "double", "finfo"}
ZERO_HEIGHT = {"zeta_zero", "zeta_zeros"}
ZERO_OTHER = {"zeta_nzeros", "riemann_siegel_z"}
TIER_NAMES = ["rational", "symbolic", "ball", "steer", "float"]
# Seeded samples whose sources were read by hand; verdicts are in the companion note.
REVIEWED_EXACT = ["ch39-oriented-height-reconstruction.py", "ch39-native-phase-alignment.py", "ch39-clock-baryon-anomaly-free.py",
                  "ch39-rh-polarity-completion.py", "ch39-shell-bulk-parent-couplings.py", "ch39-clock-schur-mediators.py",
                  "ch39-two-preparation-pole.py", "ch39-four-dimensional-plaquette.py", "ch39-higgs-baryon-triality-measure.py",
                  "ch39-rh-shell-two-thresholds.py", "ch39-finite-shell-first-triplet.py", "ch39-rh-shell-euclidean-distance.py"]
# Hand verdicts on the float-tier sample (what the floating arithmetic is used for).
FLOAT_VERDICTS = {"ch39-wall-maxwell-moments.py": "float decides",
                  "ch39-common-source-tower.py": "float proposes, ball certifies",
                  "ch39-shell-finite-order-target.py": "float proposes, ball certifies",
                  "ch39-shell-gw-infinitesimal-mixing.py": "float proposes, ball certifies",
                  "ch39-majorana-vacuum-energy.py": "float proposes, ball certifies",
                  "ch39-shell-gw-conjugate-cost.py": "float constants as model input",
                  "ch39-rh-reference-kinetic-continuum.py": "tolerance or step count",
                  "ch39-shell-gw-nonlinear-crossing.py": "tolerance or step count",
                  "ch39-shell-baseband-explicit.py": "exact imaginary unit 1j",
                  "ch39-hinge-occupation-response.py": "float literal as ball comparison bound"}
# Closure-rational scripts with zero class 0 that nevertheless take zero-derived data from hard-coded
# constants or from a table no zero-calling script reads (known false negatives of the D rule).
HIDDEN_ZERO_DATA = ["ch39-dyadic-shortlist.py", "ch39-finite-rh-pivot-mass.py", "ch39-rational-local-flavor.py", "ch39-rational-gauge-flow.py"]
REVIEWED_FLOAT = ["ch39-common-source-tower.py", "ch39-shell-gw-nonlinear-crossing.py", "ch39-rh-reference-kinetic-continuum.py",
                  "ch39-wall-maxwell-moments.py", "ch39-shell-finite-order-target.py", "ch39-shell-gw-infinitesimal-mixing.py",
                  "ch39-hinge-occupation-response.py", "ch39-majorana-vacuum-energy.py", "ch39-shell-baseband-explicit.py",
                  "ch39-shell-gw-conjugate-cost.py"]


def analyze(path, names):
    tree = ast.parse(SOURCES[path.name] if path.name in SOURCES else path.read_text(encoding="utf-8"))
    parent = {c: p for p in ast.walk(tree) for c in ast.iter_child_nodes(p)}

    def role(n):
        c = n
        while c in parent:
            p = parent[c]
            if isinstance(p, ast.Assert) and c is p.msg:
                return "display"
            if isinstance(p, ast.keyword) and p.arg in ("key", "tol"):
                return "steer"
            if isinstance(p, ast.Call) and isinstance(p.func, ast.Name) and p.func.id in ("int", "round", "range") and c in p.args:
                return "steer"
            if isinstance(p, ast.JoinedStr):
                return "display"
            if isinstance(p, ast.Call):
                f = p.func
                if isinstance(f, ast.Name) and f.id == "print":
                    return "display"
                if isinstance(f, ast.Attribute) and f.attr in ("format", "write", "dump", "dumps"):
                    return "display"
            c = p
        return "float"

    ev = {k: set() for k in ("float", "steer", "display", "ball", "symbolic", "exact")}
    zh = zo = 0
    alias = {}
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                alias[a.asname or a.name.split(".")[0]] = a.name
                if a.name.split(".")[0] in FLOAT_MOD:
                    ev["float"].add("import " + a.name)
        elif isinstance(n, ast.ImportFrom) and n.module:
            top = n.module.split(".")[0]
            for a in n.names:
                alias[a.asname or a.name] = n.module + "." + a.name
                if top in FLOAT_MOD:
                    ev["float"].add("from %s import %s" % (n.module, a.name))
                if top == "math" and a.name in MATH_FLOAT:
                    ev["float"].add("math." + a.name)
                if top == "flint" and a.name in BALL:
                    ev["ball"].add("flint." + a.name)
                if top == "flint" and a.name in EXACT_FLINT:
                    ev["exact"].add("flint." + a.name)
                if top == "sympy" and a.name in SYM_FLOAT:
                    ev["float"].add("sympy." + a.name)
                elif top == "sympy" and a.name in SYM_IRR:
                    ev["symbolic"].add("sympy." + a.name)

    def module(expr):
        if isinstance(expr, ast.Name):
            return alias.get(expr.id, "")
        if isinstance(expr, ast.Attribute):
            b = module(expr.value)
            return b + "." + expr.attr if b else ""
        return ""

    npnames = {k for k, v in alias.items() if v.split(".")[0] == "numpy"}
    for n in ast.walk(tree):
        if isinstance(n, ast.Attribute):
            top = module(n.value).split(".")[0]
            if top == "math" and n.attr in MATH_FLOAT:
                ev["float"].add("math." + n.attr)
            if top == "flint" and n.attr in BALL:
                ev["ball"].add("flint." + n.attr)
            if top == "flint" and n.attr in EXACT_FLINT:
                ev["exact"].add("flint." + n.attr)
            if top == "sympy":
                if n.attr in SYM_FLOAT:
                    ev["float"].add("sympy." + n.attr)
                elif n.attr in SYM_IRR:
                    ev["symbolic"].add("sympy." + n.attr)
            if top == "numpy" and isinstance(n.value, ast.Name) and n.value.id in npnames and n.attr in NP_FLOAT:
                ev["float"].add("numpy." + n.attr)
            if n.attr == "evalf":
                ev["float"].add(".evalf")
            if n.attr in ZERO_HEIGHT:
                zh += 1
            if n.attr in ZERO_OTHER:
                zo += 1
        elif isinstance(n, ast.Name):
            full = alias.get(n.id, "")
            if full.startswith("flint.") and full.split(".")[1] in BALL:
                ev["ball"].add(full)
        elif isinstance(n, ast.Constant) and type(n.value) in (float, complex):
            ev[role(n)].add(type(n.value).__name__ + " literal")
        elif isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in ("float", "complex") and n.func.id not in alias:
            ev[role(n)].add(n.func.id + "()")
        elif isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr in ("random", "uniform", "gauss", "normalvariate", "expovariate") and not n.args:
            ev["float"].add("random float draw")
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div) and all(isinstance(x, ast.Constant) and type(x.value) is int for x in (n.left, n.right)):
            ev[role(n)].add("int/int literal division")
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Pow) and isinstance(n.right, ast.Constant) and type(n.right.value) is float:
            ev[role(n)].add("float power")
        elif isinstance(n, ast.keyword) and n.arg == "dtype" and isinstance(n.value, ast.Name) and n.value.id in ("float", "complex") and npnames:
            ev["float"].add("numpy dtype=" + n.value.id)
        elif isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "astype" and n.args and isinstance(n.args[0], ast.Name) and n.args[0].id in ("float", "complex"):
            ev["float"].add("astype(" + n.args[0].id + ")")
    if npnames:
        ev["display"].add("numpy imported")

    consts = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
    deps = {m for s in consts for m in re.findall(r"ch39-[A-Za-z0-9_-]+\.py", s) if m in names}
    dynamic = "ch39-" in consts or any(isinstance(n, ast.JoinedStr) and any(isinstance(v, ast.Constant) and "ch39-" in v.value for v in n.values) for n in ast.walk(tree))
    if dynamic:
        deps |= {"ch39-" + s + ".py" for s in consts if "ch39-" + s + ".py" in names}
    deps.discard(path.name)
    data = {m for s in consts for m in re.findall(r"[A-Za-z0-9_一-鿿-]+\.(?:json|txt|tsv)", s) if (R / m).exists()}
    data.discard(STORE.name)
    if ev["float"]:
        tier = 4
    elif ev["steer"]:
        tier = 3
    elif ev["ball"]:
        tier = 2
    elif ev["symbolic"]:
        tier = 1
    else:
        tier = 0
    return {"tier": tier, "deps": sorted(deps), "data": sorted(data), "zh": zh, "zo": zo,
            **{k: sorted(v) for k, v in ev.items()}}


def census(population):
    names = set(population)
    raw = {k: analyze(R / k, names) for k in population}
    closure = {}
    for k in population:
        seen, stack = {k}, [k]
        while stack:
            for y in raw[stack.pop()]["deps"]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        closure[k] = seen
    ctier = {k: max(raw[x]["tier"] for x in closure[k]) for k in population}
    calls = {k: any(raw[x]["zh"] or raw[x]["zo"] for x in closure[k]) for k in population}
    zero_data = sorted({d for k in population if calls[k] for d in raw[k]["data"]})
    zero = {}
    for k in population:
        if any(raw[x]["zh"] for x in closure[k]):
            zero[k] = "H"
        elif any(raw[x]["zo"] for x in closure[k]):
            zero[k] = "N"
        elif any(d in zero_data for x in closure[k] for d in raw[x]["data"]):
            zero[k] = "D"
        else:
            zero[k] = "0"
    return raw, closure, ctier, zero, zero_data


def main():
    write = "--write" in sys.argv
    if write:
        import subprocess
        out = subprocess.run(["git", "ls-tree", "-r", "--name-only", COMMIT, "--", "./"], cwd=R, capture_output=True, text=True, check=True).stdout
        # ls-tree prints paths relative to the working directory R; keep the top level of replay/ only.
        population = sorted(x for x in out.split("\n") if x.startswith("ch39-") and "/" not in x and x.endswith(".py"))
        digests = {k: hashlib.sha256(git_blob(k)).hexdigest() for k in population}
    else:
        stored = json.loads(STORE.read_text(encoding="utf-8"))
        digests = stored["population"]
        population = sorted(digests)
    load_sources(digests)
    assert len(population) == 1706 and SELF not in population
    raw, closure, ctier, zero, zero_data = census(population)
    direct = Counter(raw[k]["tier"] for k in population)
    ctotal = Counter(ctier.values())
    ztotal = Counter(zero.values())
    cross = Counter((ctier[k], zero[k]) for k in population)
    print("direct tiers", {TIER_NAMES[t]: direct[t] for t in range(5)})
    print("closure tiers", {TIER_NAMES[t]: ctotal[t] for t in range(5)})
    print("zero dependence of closure", dict(sorted(ztotal.items())))
    for t in range(5):
        print("  %-8s" % TIER_NAMES[t], {z: cross[t, z] for z in "HNDO".replace("O", "0")})
    assert [direct[t] for t in range(5)] == [328, 34, 1114, 123, 107]
    assert [ctotal[t] for t in range(5)] == [231, 38, 423, 482, 532]
    assert (ztotal["H"], ztotal["N"], ztotal["D"], ztotal["0"]) == (1283, 0, 112, 311)
    assert sum(1 for k in population if raw[k]["zh"]) == 537
    assert len(zero_data) == 38
    assert sum(1 for k in population if raw[k]["deps"]) == 1456
    assert {(t, z): cross[t, z] for t in range(5) for z in "HND0" if cross[t, z]} == {
        (0, "D"): 79, (0, "0"): 152, (1, "0"): 38, (2, "H"): 289, (2, "D"): 33, (2, "0"): 101,
        (3, "H"): 472, (3, "0"): 10, (4, "H"): 522, (4, "0"): 10}
    exact = sorted(k for k in population if ctier[k] == 0)
    exact_zero = [k for k in exact if zero[k] != "0"]
    assert len(exact_zero) == 79 and all(zero[k] == "D" for k in exact_zero)
    # Seeded spot checks; the verdicts from reading the sources are recorded in the companion note.
    sample_exact = random.Random(20261126).sample(exact, 12)
    sample_float = random.Random(20261127).sample(sorted(k for k in population if raw[k]["tier"] == 4), 10)
    assert sample_exact == REVIEWED_EXACT and sample_float == REVIEWED_FLOAT
    assert sum(zero[k] == "D" for k in sample_exact) == 7 and sum(zero[k] == "H" for k in sample_float) == 8
    assert sorted(FLOAT_VERDICTS) == sorted(sample_float)
    assert sorted(Counter(FLOAT_VERDICTS.values()).values()) == [1, 1, 1, 1, 2, 4]
    assert Counter(FLOAT_VERDICTS.values())["float decides"] == 1
    assert all(ctier[k] == 0 and zero[k] == "0" for k in HIDDEN_ZERO_DATA)
    # A known false negative of the float rule: a division of integer expressions, harmless here.
    assert ctier["ch39-rational-band-lift.py"] == 0 and "(20*36-183)/537 == 1" in SOURCES["ch39-rational-band-lift.py"]
    print("exact-tier sample", sample_exact)
    print("float-tier sample", sample_float)
    lists = {"rules": __doc__,
             "tier_names": TIER_NAMES,
             "counts": {"direct": [direct[t] for t in range(5)], "closure": [ctotal[t] for t in range(5)],
                        "zero": {z: ztotal[z] for z in "HND0"},
                        "closure_by_zero": {TIER_NAMES[t]: {z: cross[t, z] for z in "HND0"} for t in range(5)}},
             "zero_related_data_files": zero_data,
             "sample_exact": sample_exact, "sample_float": sample_float,
             "scripts": {k: {"direct": TIER_NAMES[raw[k]["tier"]], "closure": TIER_NAMES[ctier[k]], "zero": zero[k],
                             "closure_size": len(closure[k]), "float": raw[k]["float"], "steer": raw[k]["steer"],
                             "ball": raw[k]["ball"], "symbolic": raw[k]["symbolic"]} for k in population},
             "population": digests}
    if write:
        STORE.write_text(json.dumps(lists, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    else:
        stored = json.loads(STORE.read_text(encoding="utf-8"))
        assert stored["scripts"] == lists["scripts"] and stored["counts"] == lists["counts"]
        assert stored["sample_exact"] == sample_exact and stored["sample_float"] == sample_float
    print("closure-rational scripts:", ctotal[0], "of which read zero-related data:", len(exact_zero))
    classification(lists["scripts"])


FIRST_ROW, ROUNDS_ROW, LAST_ROW = "213 个原生两面字要移动多少概率", "沿对角时间的交叉偏振负比率支", "截止 56 时，保持三坐标置换对称"
BOOK_REPO = R.parents[2] / ".." / "trureturing-books"
BOOK_COMMIT = "cbdd318582de6c2f71e48d9d072fb8e07d94905a"
BOOK_CH15 = "books/omega-unification/zh/ch15-unified-verdict-and-open-roadmap.md"


def book_ch15():
    """Chapter 15 at the pinned book commit, or None when that commit is not available."""
    import subprocess
    try:
        out = subprocess.run(["git", "show", "%s:%s" % (BOOK_COMMIT, BOOK_CH15)], cwd=BOOK_REPO, capture_output=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.decode("utf-8")
R4_REASONS = {"两量恰好相等", "全体零点量词", "选择问题", "残差远小于包围", "缺结构原理"}


def table_questions():
    """Questions of the README summary table from FIRST_ROW to LAST_ROW, as recorded at COMMIT."""
    import subprocess
    try:
        text = subprocess.run(["git", "show", "%s:./../README.md" % COMMIT], cwd=R, capture_output=True, check=True).stdout.decode("utf-8")
    except (OSError, subprocess.CalledProcessError):
        text = (R.parent / "README.md").read_text(encoding="utf-8")
    rows = [line[2:].split(" | ")[0] for line in text.split("\n") if line.startswith("| ") and " | " in line]
    i = next(k for k, q in enumerate(rows) if q.startswith(FIRST_ROW))
    j = next(k for k, q in enumerate(rows) if q.startswith(LAST_ROW))
    return rows[i:j + 1], text


def classification(scripts):
    """Part two: the written R1-R4 classification (a judgement, not a certificate) is checked for form and counted."""
    from fractions import Fraction
    from math import isqrt
    entries = json.loads((R / "ch39-rational-map-classification.json").read_text(encoding="utf-8"))
    questions, readme = table_questions()
    assert len(questions) == 25
    for x in entries:
        assert x["class"] in ("R1", "R2", "R3", "R4") and x["claim"] and x["rational_form"] and x["certificate"]
        assert all(k in scripts for k in x["scripts"]), x["scripts"]
        ps = sorted({int(v) for k in x["scripts"] for v in re.findall(r"ctx\.prec *= *(\d+)", SOURCES[k])})
        expected = ("无脚本" if not x["scripts"] else "脚本自身未设" if not ps else
                    "%d 位" % ps[0] if len(ps) == 1 else "%d–%d 位" % (ps[0], ps[-1]))
        assert x["working_precision"] == expected, (x["claim"], x["working_precision"], expected)
        if x["class"] == "R2":
            assert x["r2_kind"] in ("全称", "存在见证", "环面子区域", "以全称带界为前提的推导"), x["claim"]
        if x["class"] == "R3":
            assert x.get("bridge"), x["claim"]
        if x["class"] == "R4":
            assert x.get("why_not_resolution") and x.get("r4_reason") and set(x["r4_reason"]) <= R4_REASONS, x["claim"]
        if x.get("zero_use") == "有限枚零点的数值位置" and x["scripts"]:
            assert any(scripts[k]["zero"] == "H" for k in x["scripts"]), x["claim"]
        if x["source"] == "README":
            assert sum(q.startswith(x["row"]) for q in questions) == 1, x["row"]
            assert x["note"] and all((R / n).exists() for n in x["note"]), x["note"]
        assert not re.search(r"README:\d", json.dumps(x, ensure_ascii=False))
    by = lambda xs: [sum(1 for x in xs if x["class"] == c) for c in ("R1", "R2", "R3", "R4")]
    rows = [x for x in entries if x["source"] == "README"]
    probs = [x for x in entries if x["source"].startswith("ch15:")]
    extra = [x for x in entries if x["source"].startswith("extra:")]
    assert len(rows) + len(probs) + len(extra) == len(entries) == 91 and len(extra) == 1 and extra[0]["class"] == "R3"
    assert {next(q for q in questions if q.startswith(x["row"])) for x in rows} == set(questions)
    k = next(i for i, q in enumerate(questions) if q.startswith(ROUNDS_ROW))
    rounds = [x for x in rows if any(q.startswith(x["row"]) for q in questions[k:])]
    assert len(questions) - k == 21 and len(rounds) == 57
    assert by(rows) == [19, 13, 26, 5] and by(rounds) == [17, 13, 22, 5] and by(probs) == [4, 1, 6, 16]
    assert by(entries) == [23, 14, 33, 21]
    assert Counter(x["r2_kind"] for x in entries if x["class"] == "R2") == {"全称": 8, "存在见证": 3, "以全称带界为前提的推导": 2, "环面子区域": 1}
    per = {q: by([x for x in probs if x["source"] == "ch15:prob-ou-q%d:gap" % q]) for q in range(1, 5)}
    assert per == {1: [0, 0, 0, 4], 2: [0, 0, 1, 3], 3: [1, 0, 1, 4], 4: [1, 0, 0, 4]}
    gaps = [x for x in probs if x["source"].endswith(":gap")]
    assert all(x.get("missing_sentence") for x in gaps) and not any(x.get("missing_sentence") for x in extra)
    assert Counter(r for x in entries if x["class"] == "R4" for r in x["r4_reason"]) == {"缺结构原理": 16, "选择问题": 5, "全体零点量词": 2}
    book = book_ch15()
    if book is not None:
        assert all(m in book for x in gaps for m in x["missing_sentence"]), "ch15 quote"
        assert "消除对低能输入偏移" in book
        print("ch15 quotes checked against book commit", BOOK_COMMIT[:12])
    else:
        print("ch15 quotes not checked: pinned book commit not available next to this repository")
    rh = [x for x in entries if x.get("zero_use") == "用到对全体零点的分布断言"]
    assert len(rh) == 2 and all(x["class"] == "R4" for x in rh)
    precs = sorted({int(v) for x in entries if x["class"] == "R3" for k in x["scripts"]
                    for v in re.findall(r"ctx\.prec *= *(\d+)", SOURCES[k])})
    assert precs[0] == 112 and precs[-1] == 576, precs
    # The circle point u=(108469+109*sqrt(867) i)/|.| of the word problems has irrational cosine.
    n = 108469**2 + 867*109**2
    assert n == 11775824788 and isqrt(n)**2 != n
    # 4.21603 and 29.984 are only checked to appear in README; their ratio is computed exactly.
    assert 7 < Fraction("29.984")/Fraction("4.21603") < Fraction(72, 10)
    for text in ("约 120 位", "\\ge1/5", "4.21603", "29.984", "6.59\\times10^4", "1.8\\times10^{-120}", "0.0398046", "-17/282",
                 "29\\le n_{\\rm aux}\\le31", "102 枚", "6.6\\times10^{19}", "另 23 枚未判定", "另 8 枚", "#12610", "最坏名次"):
        assert text in readme, text
    assert "8.685<9.6" in (R / "十九点壳重叠核与测度输运项的显式指数局域性.md").read_text(encoding="utf-8")
    assert "相对半径约 $10^{-12}$" in (R / "两点截止在匹配牛顿项时不能驻定真空读数而四点测度对共同截止标度一阶驻定.md").read_text(encoding="utf-8")
    assert "0.2976" in (R / "每腿至多两次邻边输运时第九项即解除三零点共同系数的质量比障碍.md").read_text(encoding="utf-8")
    assert "Q(3141,1000)" in SOURCES["ch39-shell-bulk-parent-couplings.py"] and "Q(22,7)" in SOURCES["ch39-shell-bulk-parent-couplings.py"]
    assert "1.72781" in (R / "六项算术时间强迫的第八至第十二项二进制码.md").read_text(encoding="utf-8")
    # Quoted anchors of the RH/zero-count section must be present in the cited sources.
    for name, text in (("ch39-scale-weight-contrast.py", "0.28 log(t)"),
                       ("Ω复算说明.md", "t[log(t/(2π))−1]/(2π)+7/8"),
                       ("Ω复算说明.md", "Brent、Platt、Trudgian"),
                       ("Ω复算说明.md", "0.00093642<W<0.00771815"),
                       ("Ω复算说明.md", "由 RH 下完成 ξ 函数的乘积展开"),
                       ("Ω复算说明.md", "在 RH 下，经典的 ξ 函数对数导数恒等式给出全部正零点的权重和"),
                       ("Ω复算说明.md", "由 RH 的成对零点表示得到"),
                       ("Ω复算说明.md", "黎曼–冯·曼戈尔特零点计数公式"),
                       ("ch39-双时钟的三代手征环面.md", "由 RH 和完成函数的对称 Hadamard 乘积"),
                       ("ch39-双时钟的三代手征环面.md", "经典 RH 给 $\\rho(1-\\rho)=1/4+\\gamma^2$"),
                       ("三窗味概率下CP亏额与Vcb偏移之比的双侧界.md", "RH 下的双钟联合等分布"),
                       ("两步质量读数下容许秩含二三四时五荷一致的三代费米能级窗为空而全环面共同两步分界只剩秩15.md", "零点相位在环面上稠密"),
                       ("ch39-零点与标准模型规范群积分.md", "Landau 零点显式公式"),
                       ("ch39-双时钟的三代手征环面.md", "Guinand–Weil 显式公式")):
        assert text in (R / name).read_text(encoding="utf-8"), (name, text)
    print("classification entries", len(entries), "table rows:", by(rows), "of which from the reflection row on:", by(rounds),
          "ch15:", by(probs), "per-problem gaps:", per, "(R1,R2,R3,R4)")


if __name__ == "__main__":
    main()
