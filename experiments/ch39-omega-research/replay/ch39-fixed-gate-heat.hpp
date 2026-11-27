// 两种四角覆盖共用的定点谱投影、加权热矩与向外整数运算。
#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

using Z = std::int64_t;
constexpr Z Q = 1 << 24;
constexpr int period = 1 << 16;
using Pair = std::array<Z, 2>;
struct Row {
  std::array<Pair, 5> coeff;
  std::array<Z, 4> size{};
  Z error;
};
struct Span {
  int lo, hi;
};
using Box = std::array<Span, 4>;
struct Shell {
  std::array<Row, 38> rows;
  std::array<Pair, 19> weights;
  std::array<Pair, 2> goal;
};
std::array<Pair, 3> gain;
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


template <std::size_t N>
bool heat_separates(const Shell &shell, const std::array<Pair, 4> &gate,
                    const std::array<Z, N> &changes) {
  std::array<Pair, 2> sum{};
  for (int i = 0; i < 19; ++i) {
    Pair amp = gain[0];
    for (int j = 0; j < 2; ++j) {
      const Row &row = shell.rows[2 * i + j];
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
      Z motion = 0;
      for (std::size_t k = 0; k < N; ++k)
        motion += row.size[k] * changes[k];
      Z radius = row.error + up(motion, Q);
      Z l = std::max<Z>(0, mag - radius), r = std::min(Q, mag + 1 + radius);
      assert(l <= r);
      amp[0] += gain[j + 1][0] * (l * l / Q) / Q;
      amp[1] += up(gain[j + 1][1] * up(r * r, Q), Q);
    }
    amp[1] = std::min(cap, amp[1]);
    Pair x{amp[0] * amp[0] / Q, up(amp[1] * amp[1], Q)};
    for (int k = 0; k < 2; ++k) {
      sum[k][0] += shell.weights[i][0] * (shell.weights[i][0] < 0 ? x[1] : x[0]);
      sum[k][1] += shell.weights[i][1] * (shell.weights[i][1] < 0 ? x[0] : x[1]);
      if (k == 0)
        x = {x[0] * x[0] / Q, up(x[1] * x[1], Q)};
    }
  }
  for (int k = 0; k < 2; ++k)
    if (sum[k][0] > shell.goal[k][1] * Q || sum[k][1] < shell.goal[k][0] * Q)
      return true;
  return false;
}

template <std::size_t N>
void read_projections(Shell &shell) {
  Z total_weight = 0;
  for (auto &a : shell.weights) {
    std::cin >> a[0] >> a[1];
    total_weight += std::max(std::abs(a[0]), std::abs(a[1]));
  }
  // 带权累加前的绝对上界为 2048 Q²，落在有符号 64 位整数范围内。
  assert(total_weight < 2048 * Q);
  for (auto &r : shell.rows) {
    for (auto &t : r.coeff) {
      std::cin >> t[0] >> t[1];
      assert(std::abs(t[0]) <= 2 * Q && std::abs(t[1]) <= 2 * Q);
    }
    for (std::size_t k = 0; k < N; ++k) {
      std::cin >> r.size[k];
      assert(0 <= r.size[k] && r.size[k] <= 4 * Q);
    }
    std::cin >> r.error;
    assert(0 <= r.error && r.error < Q);
  }
}
