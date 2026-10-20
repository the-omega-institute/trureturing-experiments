#!/usr/bin/env python3
"""Kumar-Peng-Tran, arXiv:2211.04379v1, Conjecture 4.4: E(m, 2) is the third most inefficient for odd m, where
E(m, s) = max_r ||r||_s and ||r||_s is the least number of x's in a (+, *)-expression evaluating to r mod m at x = s.
Exact layered closure: L_1 = {s}; L_n = (union over i of L_i + L_{n-i} and L_i * L_{n-i}) minus earlier layers. A shortest
expression with n > 1 leaves has shortest children, so ||r|| = the n with r in L_n. For every odd m <= 101, counts the unit
bases s with E(m, s) > E(m, 2); prints full rankings for m <= 17 and the minimal failure. usage: check.py"""
import sys, math

def layers(m, s):
    best = {s % m: 1}; L = [None, {s % m}]; n = 1
    while len(best) < m:
        n += 1; new = set()
        for i in range(1, n):
            for a in L[i]:
                for b in L[n - i]:
                    for r in ((a + b) % m, (a * b) % m):
                        if r not in best: new.add(r)
        for r in new: best[r] = n
        L.append(new)
        if n > 4 * m: return None, L  # unreachable residues
    return best, L

def E(m, s):
    best, _ = layers(m, s)
    return None if best is None else max(best.values())

fails = []
for m in range(3, 102, 2):
    units = [s for s in range(m) if math.gcd(s, m) == 1]
    Em = {s: E(m, s) for s in units}
    above = sorted(s for s in units if Em[s] > Em[2])
    if len(above) > 2: fails.append(m)
    if m <= 17:
        rank = sorted(((v, s) for s, v in Em.items()), reverse=True)
        print(f"m={m}: E(m,2)={Em[2]}; unit bases strictly above 2: {[(s, Em[s]) for s in above]}; top: {rank[:5]}")
print("odd m <= 101 with at least three unit bases strictly above base 2:", fails)
m = 17; best2, L2 = layers(m, 2); best3, L3 = layers(m, 3)
print("m=17 base 2 layers:", [sorted(x) for x in L2[1:] if x])
print("m=17 base 3 layers:", [sorted(x) for x in L3[1:] if x])
ok = fails and fails[0] == 17 and E(17, 1) == 9 and E(17, 16) == 7 and E(17, 3) == 7 and E(17, 2) == 6
print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
