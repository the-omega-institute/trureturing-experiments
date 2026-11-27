// 四个角的二进闭盒，所有决定覆盖的包围均以定点整数向外取整。
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using Z = std::int64_t;
constexpr Z Q = 1 << 24;
constexpr int period = 1 << 16;
using Pair = std::array<Z, 2>;
struct Row {
  std::array<Pair, 5> coeff;
  Z diag, off, error;
};
struct Span {
  int lo, hi;
};
using Box = std::array<Span, 4>;
std::array<Row, 38> rows;
std::array<Pair, 19> weights;
std::array<Pair, 3> gain;
std::array<Pair, 2> goal;
std::array<Z, period / 4 + 1> quarter;
Z cap;
Z up(Z a, Z b) {
  assert(a >= 0);
  return (a + b - 1) / b;
}
Z sine(int n) {
  n %= period;
  if (n < 0)
    n += period;
  if (n > period / 2)
    return -sine(n - period / 2);
  if (n > period / 4)
    n = period / 2 - n;
  return quarter[n];
}
Z cosine(int n) { return sine(n + period / 4); }
Z root_floor(Z n) {
  Z r = static_cast<Z>(std::sqrt(static_cast<double>(n)));
  while ((r + 1) * (r + 1) <= n)
    ++r;
  while (r * r > n)
    --r;
  return r;
}

bool excluded(const Box &box, int &split) {
  std::array<int, 4> m;
  std::array<Z, 4> h;
  for (int j = 0; j < 4; ++j) {
    int width = box[j].hi - box[j].lo;
    assert(width >= 2 && width % 2 == 0);
    m[j] = (box[j].hi + box[j].lo) / 2;
    // 半宽为 pi * width / period，355/113 > pi。
    h[j] = up(Q * 355 * width, 113 * period);
  }
  Z c = cosine(m[0]), s = sine(m[0]);
  Z cd = std::min(Q, std::abs(c) + 1 + h[0]),
    sd = std::min(Q, std::abs(s) + 1 + h[0]);
  Z diag = std::min(2 * Q, up(sd * h[0] + cd * (h[1] + h[3]), Q));
  Z off = std::min(2 * Q, up(cd * h[0] + sd * (h[2] + h[3]), Q));
  std::array<Z, 4> priority{h[0], up(cd * h[1], Q), up(sd * h[2], Q), h[3]};
  split = static_cast<int>(std::max_element(priority.begin(), priority.end()) -
                           priority.begin());
  std::array<Pair, 4> gate;
  std::array<int, 4> phases{m[3] + m[1], m[3] - m[1], m[3] + m[2], m[3] - m[2]};
  for (int k = 0; k < 4; ++k) {
    Z t = k < 2 ? c : s;
    gate[k] = {t * cosine(phases[k]) / Q, t * sine(phases[k]) / Q};
  }
  std::array<Pair, 2> sum{};
  for (int i = 0; i < 19; ++i) {
    Pair amp = gain[0];
    for (int j = 0; j < 2; ++j) {
      const Row &row = rows[2 * i + j];
      Z vr = row.coeff[0][0] * Q, vi = row.coeff[0][1] * Q;
      for (int k = 0; k < 4; ++k) {
        const auto &a = row.coeff[k + 1];
        const auto &g = gate[k];
        vr += a[0] * g[0] - a[1] * g[1];
        vi += a[0] * g[1] + a[1] * g[0];
      }
      vr /= Q;
      vi /= Q;
      Z mag = root_floor(vr * vr + vi * vi);
      Z radius = row.error + up(row.diag * diag + row.off * off, Q);
      Z l = std::max<Z>(0, mag - radius), r = std::min(Q, mag + 1 + radius);
      assert(l <= r);
      amp[0] += gain[j + 1][0] * (l * l / Q) / Q;
      amp[1] += up(gain[j + 1][1] * up(r * r, Q), Q);
    }
    amp[1] = std::min(cap, amp[1]);
    Pair x{amp[0] * amp[0] / Q, up(amp[1] * amp[1], Q)};
    for (int k = 0; k < 2; ++k) {
      sum[k][0] += weights[i][0] * (weights[i][0] < 0 ? x[1] : x[0]);
      sum[k][1] += weights[i][1] * (weights[i][1] < 0 ? x[0] : x[1]);
      if (k == 0)
        x = {x[0] * x[0] / Q, up(x[1] * x[1], Q)};
    }
  }
  for (int k = 0; k < 2; ++k)
    if (sum[k][0] > goal[k][1] * Q || sum[k][1] < goal[k][0] * Q)
      return true;
  return false;
}

int main() {
  Z unit, cycle;
  std::cin >> unit >> cycle >> cap;
  assert(unit == Q && cycle == period && cap > 0 && cap < Q);
  for (auto &a : goal)
    std::cin >> a[0] >> a[1];
  for (auto &a : gain) {
    std::cin >> a[0] >> a[1];
    assert(0 <= a[0] && a[1] <= Q);
  }
  Z total_weight = 0;
  for (auto &a : weights) {
    std::cin >> a[0] >> a[1];
    total_weight += std::max(std::abs(a[0]), std::abs(a[1]));
  }
  assert(total_weight < 1200 * Q);
  for (auto &r : rows) {
    for (auto &t : r.coeff) {
      std::cin >> t[0] >> t[1];
      assert(std::abs(t[0]) <= 2 * Q && std::abs(t[1]) <= 2 * Q);
    }
    std::cin >> r.diag >> r.off >> r.error;
    assert(r.diag <= 4 * Q && r.off <= 4 * Q && r.error < Q);
  }
  for (auto &s : quarter) {
    std::cin >> s;
    assert(s >= 0 && s <= Q);
  }
  assert(std::cin);
  std::vector<Box> todo{{Span{0, period / 4}, Span{-period / 2, period / 2},
                         Span{-period / 2, period / 2}, Span{0, period / 2}}};
  Z leaves = 0;
  while (!todo.empty()) {
    Box box = todo.back();
    todo.pop_back();
    int j = 0;
    if (excluded(box, j)) {
      ++leaves;
      continue;
    }
    int m = (box[j].lo + box[j].hi) / 2;
    assert(box[j].hi - box[j].lo >= 4);
    Box right = box;
    right[j].lo = m;
    box[j].hi = m;
    todo.push_back(right);
    todo.push_back(box);
  }
  std::cout << "一般 U(2) 的认证闭盒数 " << leaves
            << "，至少一项绝对残差严格大于 0.01\n";
}
