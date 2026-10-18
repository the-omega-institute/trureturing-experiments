/* Binary Lyndon words of length n attaining the minimum number of distinct Lyndon factors:
   print each with Christoffel flag and its standard factorization (u, v), v the longest proper Lyndon suffix,
   together with L(u), L(v), f(|u|), f(|v|) and |Lyn(u) ∪ Lyn(v)|.  usage: minimizers NMAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int is_lyndon(const int *s, int len) { int i = 0, j = 1; while (j < len) { if (s[i] < s[j]) { i = 0; j++; } else if (s[i] == s[j]) { i++; j++; } else return 0; } return i == 0; }
#define HSZ 8192
static uint64_t hk[HSZ]; static int hu[HSZ]; static int hs = 1;
static int ins(uint64_t k) { uint64_t h = (k * 0x9E3779B97F4A7C15ULL) >> 51; while (hu[h] == hs) { if (hk[h] == k) return 0; h = (h + 1) & (HSZ - 1); } hu[h] = hs; hk[h] = k; return 1; }
static int addfactors(const int *w, int len) { int c = 0; for (int i = 0; i < len; i++) { uint64_t b = 0; for (int j = i; j < len; j++) { b = (b << 1) | (uint64_t)w[j]; int l = j - i + 1; if (is_lyndon(w + i, l)) c += ins(((uint64_t)l << 40) | b); } } return c; }
static int L(const int *w, int len) { hs++; return addfactors(w, len); }
static int f[64]; static int n; static int a[64];
static int store[200000][40]; static int ns;
static int best;
static void visit(void) { int c = L(a + 1, n); if (c < best) { best = c; ns = 0; } if (c == best && ns < 200000) memcpy(store[ns++], a + 1, sizeof(int) * n); }
static void gen(int t, int p) { if (t > n) { if (p == n) visit(); return; } a[t] = a[t - p]; gen(t + 1, p); for (int j = a[t - p] + 1; j < 2; j++) { a[t] = j; gen(t + 1, t); } }
static int christoffel(const int *w, int len) { int p = 0; for (int i = 0; i < len; i++) p += w[i]; for (int i = 0; i < len; i++) if (w[i] != ((i + 1) * p) / len - (i * p) / len) return 0; return 1; }
int main(int argc, char **argv) {
  int nmax = atoi(argv[1]); f[1] = 1;
  for (n = 2; n <= nmax; n++) {
    best = 1 << 30; ns = 0; a[0] = 0; gen(1, 1); f[n] = best;
    printf("n=%d f=%d minimizers=%d\n", n, best, ns);
    for (int s = 0; s < ns; s++) {
      int *w = store[s]; int k;
      for (k = 1; k < n; k++) if (is_lyndon(w + k, n - k)) break; /* v = w[k..], longest proper Lyndon suffix */
      int lu = L(w, k), lv = L(w + k, n - k);
      hs++; int un = addfactors(w, k); un += addfactors(w + k, n - k);
      printf("  ");
      for (int i = 0; i < n; i++) putchar(w[i] ? 'b' : 'a');
      printf(" %s u=", christoffel(w, n) ? "C" : "-");
      for (int i = 0; i < k; i++) putchar(w[i] ? 'b' : 'a');
      printf(" v="); for (int i = k; i < n; i++) putchar(w[i] ? 'b' : 'a');
      printf(" L(u)=%d f(|u|)=%d L(v)=%d f(|v|)=%d |Lyn(u)uLyn(v)|=%d\n", lu, f[k], lv, f[n - k], un);
    }
    fflush(stdout);
  }
  return 0;
}
