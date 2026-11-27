// 同一二阶酉门的四角覆盖。
#include "ch39-fixed-gate-heat.hpp"
#include <vector>

Shell shell;

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
  return heat_separates(shell, gate, std::array<Z,2>{diag, off});
}

int main() {
  Z unit, cycle;
  std::cin >> unit >> cycle >> cap;
  assert(unit == Q && cycle == period && cap > 0 && cap < Q);
  for (auto &a : shell.goal)
    std::cin >> a[0] >> a[1];
  for (auto &a : gain) {
    std::cin >> a[0] >> a[1];
    assert(0 <= a[0] && a[1] <= Q);
  }
  read_projections<2>(shell);
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
