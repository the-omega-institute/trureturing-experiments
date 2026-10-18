/* Glen–Simpson–Smyth arXiv:1701.00928 Open Problem 1 over a k-letter alphabet:
   min over k-ary Lyndon words of length n of the number of distinct Lyndon factors, compared with the
   binary Christoffel minimum (Theorem 8: L = d1 + ... + dn + 3, computed here directly).
   usage: op1k K NMIN NMAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, K;
static int a[64];

static int is_lyndon(const int *s, int len) {
  int i = 0, j = 1;
  while (j < len) {
    if (s[i] < s[j]) { i = 0; j++; }
    else if (s[i] == s[j]) { i++; j++; }
    else return 0;
  }
  return i == 0;
}

#define HSZ 8192
static uint64_t hkeys[HSZ];
static int hused[HSZ];
static int hstamp = 1;
static int insert(uint64_t k) {
  uint64_t h = (k * 0x9E3779B97F4A7C15ULL) >> 51;
  while (hused[h] == hstamp) { if (hkeys[h] == k) return 0; h = (h + 1) & (HSZ - 1); }
  hused[h] = hstamp; hkeys[h] = k; return 1;
}
/* key: length (6 bits) and the factor packed in base 4 (2 bits per letter, K <= 4, length <= 28) */
static int count_distinct_lyndon(const int *w, int len) {
  int cnt = 0; hstamp++;
  for (int i = 0; i < len; i++) {
    uint64_t bits = 0;
    for (int j = i; j < len; j++) {
      bits = (bits << 2) | (uint64_t)w[j];
      int l = j - i + 1;
      if (is_lyndon(w + i, l)) cnt += insert(((uint64_t)l << 58) | bits);
    }
  }
  return cnt;
}

static int best, bestw[64], beststar;
static long long nl, nbest, nfull;
static int usesall(const int *w) { int seen = 0; for (int i = 0; i < n; i++) seen |= 1 << w[i]; return seen == (1 << K) - 1; }
static void visit(void) {
  nl++;
  int c = count_distinct_lyndon(a + 1, n);
  if (c < best) { best = c; nbest = 1; memcpy(bestw, a + 1, sizeof(int) * n); }
  else if (c == best) nbest++;
  if (usesall(a + 1)) { nfull++; if (c < beststar) beststar = c; }
}
static void gen(int t, int p) {
  if (t > n) { if (p == n) visit(); return; }
  a[t] = a[t - p]; gen(t + 1, p);
  for (int j = a[t - p] + 1; j < K; j++) { a[t] = j; gen(t + 1, t); }
}
static int gcd(int x, int y) { while (y) { int r = x % y; x = y; y = r; } return x; }

int main(int argc, char **argv) {
  K = atoi(argv[1]); int nmin = atoi(argv[2]), nmax = atoi(argv[3]);
  if (K < 2 || K > 4) return 2;
  for (n = nmin; n <= nmax; n++) {
    best = 1 << 30; beststar = 1 << 30; nbest = 0; nl = 0; nfull = 0;
    a[0] = 0; gen(1, 1);
    int cbest = 1 << 30;
    for (int p = 1; p < n; p++) {
      if (gcd(p, n) != 1) continue;
      int w[64];
      for (int i = 0; i < n; i++) w[i] = ((i + 1) * p) / n - (i * p) / n;
      int c = count_distinct_lyndon(w, n);
      if (c < cbest) cbest = c;
    }
    if (n == 1) cbest = 1;
    printf("K=%d n=%d lyndon_words=%lld min_distinct=%d attained_by=%lld min_using_all_%d_letters=%d christoffel_min=%d %s example=",
           K, n, nl, best, nbest, K, beststar == (1 << 30) ? -1 : beststar, cbest, best < cbest ? "COUNTEREXAMPLE" : "OK");
    for (int i = 0; i < n; i++) putchar('a' + bestw[i]);
    putchar('\n'); fflush(stdout);
  }
  return 0;
}
