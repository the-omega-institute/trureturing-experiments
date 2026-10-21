import importlib.util, sys
spec = importlib.util.spec_from_file_location("chk", "check.py"); chk = importlib.util.module_from_spec(spec)
sys.argv = ["x"]
src = open("check.py").read().split("colors = [1, 2, 3, 4]\nnv, es = wheel()")[0]
exec(src, chk.__dict__)
def star(r):
    nv, es = chk.wheel()
    out = []
    for j in range(r):
        m = {0: 0}
        for a in range(1, 8): m[a] = 7 * j + a
        for (u, v, c) in es: out.append((m[u], m[v], c))
    return 7 * r + 1, out
for r in (1, 2, 3):
    nv, es = star(r)
    print("star r=%d: n=%d m=%d greedy_sizes=%s" % (r, nv, len(es), sorted(chk.all_runs(nv, es, [1, 2, 3, 4]))))
