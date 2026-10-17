/* Glen–Simpson–Smyth arXiv:1701.00928, Open Problem 1: is min over binary Lyndon words of length n of
   #distinct Lyndon factors attained by a Christoffel word of length n?
   Enumerates binary Lyndon words of length n (FKM algorithm), counts distinct Lyndon factors exactly
   (hash set on (length, bits)), and compares with the lower Christoffel words a..b of slope p/n, gcd(p,n)=1.
   usage: op1 NMIN NMAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static int a[64];

static int is_lyndon(const int *s, int len) {
  /* s is Lyndon iff strictly smaller than all its proper rotations; Duval-style check */
  int i = 0, j = 1;
  while (j < len) {
    if (s[i] < s[j]) { i = 0; j++; }
    else if (s[i] == s[j]) { i++; j++; }
    else return 0;
  }
  return i == 0;
}

#define HSZ 4096
static uint64_t hkeys[HSZ];
static int hused[HSZ];
static int hstamp = 1;

static int insert(uint64_t k) {
  uint64_t h = (k * 0x9E3779B97F4A7C15ULL) >> 52;
  while (hused[h] == hstamp) {
    if (hkeys[h] == k) return 0;
    h = (h + 1) & (HSZ - 1);
  }
  hused[h] = hstamp; hkeys[h] = k; return 1;
}

static int count_distinct_lyndon(const int *w, int len) {
  int cnt = 0; hstamp++;
  for (int i = 0; i < len; i++) {
    uint64_t bits = 0;
    for (int j = i; j < len; j++) {
      bits = (bits << 1) | (uint64_t)w[j];
      int l = j - i + 1;
      if (is_lyndon(w + i, l)) {
        uint64_t key = ((uint64_t)l << 40) | bits;
        cnt += insert(key);
      }
    }
  }
  return cnt;
}

static int best;
static long long nbest, nlyndon;
static int bestw[64];

static void visit(void) {
  nlyndon++;
  int c = count_distinct_lyndon(a + 1, n);
  if (c < best) { best = c; nbest = 1; memcpy(bestw, a + 1, sizeof(int) * n); }
  else if (c == best) nbest++;
}

static void gen(int t, int p) {
  if (t > n) { if (n % p == 0 && p == n) visit(); return; }
  a[t] = a[t - p]; gen(t + 1, p);
  for (int j = a[t - p] + 1; j < 2; j++) { a[t] = j; gen(t + 1, t); }
}

static int gcd(int x, int y) { while (y) { int r = x % y; x = y; y = r; } return x; }

int main(int argc, char **argv) {
  int nmin = atoi(argv[1]), nmax = atoi(argv[2]);
  for (n = nmin; n <= nmax; n++) {
    best = 1 << 30; nbest = 0; nlyndon = 0;
    a[0] = 0; gen(1, 1);
    int cbest = 1 << 30, cp = -1;
    for (int p = 1; p < n; p++) {
      if (gcd(p, n) != 1) continue;
      int w[64];
      for (int i = 0; i < n; i++) w[i] = ((i + 1) * p) / n - (i * p) / n;
      if (!is_lyndon(w, n)) { printf("not lyndon p=%d n=%d\n", p, n); return 1; }
      int c = count_distinct_lyndon(w, n);
      if (c < cbest) { cbest = c; cp = p; }
    }
    printf("n=%d lyndon_words=%lld min_distinct=%d attained_by=%lld christoffel_min=%d (p=%d) %s example=",
           n, nlyndon, best, nbest, cbest, cp, cbest == best ? "OK" : "COUNTEREXAMPLE");
    for (int i = 0; i < n; i++) putchar(bestw[i] ? 'b' : 'a');
    putchar('\n'); fflush(stdout);
  }
  return 0;
}
