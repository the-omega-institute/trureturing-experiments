"""OEIS A343502: numbers k such that tau(tau(k)) and tau(k+1) are both prime. Dequatre's conjecture: all
members except 2 and 3 are composite.

Usage: python3 check.py N   (sieve tau up to N + 1)
"""
import sys

DATA = [2,3,4,6,8,10,15,16,22,36,46,58,82,100,106,120,166,168,178,196,210,226,256,262,270,280,312,330,346,358,378,382,408,456]

def main():
    N = int(sys.argv[1])
    tau = [0] * (N + 2)
    for d in range(1, N + 2):
        for m in range(d, N + 2, d):
            tau[m] += 1
    isp = [False, False] + [True] * N
    for i in range(2, int((N + 1) ** 0.5) + 1):
        if isp[i]:
            for m in range(i * i, N + 2, i): isp[m] = False
    def prime(x): return x < len(isp) and isp[x]
    member = [k for k in range(1, N + 1) if prime(tau[tau[k]]) and prime(tau[k + 1])]
    data = member[:len(DATA)] == DATA
    primes_in = [k for k in member if prime(k)]
    odd = [k for k in member if k % 2 == 1]
    # odd members have k + 1 = 2^(t-1) with t prime
    odd_shape = all(((k + 1) & k) == 0 and prime((k + 1).bit_length()) for k in odd)
    print(f"members <= {N}: {len(member)}; data {data}; prime members {primes_in}; odd members {odd}; odd k + 1 = 2^(t-1), t prime: {odd_shape}")
    ok = data and primes_in == [2, 3] and odd_shape
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
