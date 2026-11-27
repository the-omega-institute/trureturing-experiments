// 两种弦长的四个独立谱相位在两壳上共同取值。
#include "ch39-fixed-gate-heat.hpp"
#include <vector>

int main() {
  Z unit, cycle;
  std::cin >> unit >> cycle >> cap;
  assert(unit == Q && cycle == period && cap > 0 && cap < Q);
  for (auto &bounds : gain) {
    std::cin >> bounds[0] >> bounds[1];
    assert(0 <= bounds[0] && bounds[0] <= bounds[1] && bounds[1] <= Q);
  }
  std::array<Shell,2> shells;
  for (auto &shell : shells) {
    for (auto &bounds : shell.goal)
      std::cin >> bounds[0] >> bounds[1];
    read_projections<4>(shell);
  }
  for (auto &value : quarter) {
    std::cin >> value;
    assert(0 <= value && value <= Q);
  }
  assert(std::cin);
  const Span full{-period/2, period/2};
  std::vector<Box> pending{{full, full, full, full}};
  Z leaves = 0;
  while (!pending.empty()) {
    Box box = pending.back();
    pending.pop_back();
    std::array<Pair,4> phases;
    std::array<Z,4> changes;
    std::array<int,4> widths;
    for (int j = 0; j < 4; ++j) {
      widths[j] = box[j].hi-box[j].lo;
      assert(widths[j]>=2 && widths[j]%2==0);
      int mid = (box[j].lo+box[j].hi)/2;
      phases[j] = {cosine(mid),sine(mid)};
      // 两个单位复数的距离不超过相位差，亦不超过 2。
      changes[j] = std::min(2*Q, up(Q*355*widths[j], 113*period));
    }
    if (heat_separates(shells[0],phases,changes) || heat_separates(shells[1],phases,changes)) {
      ++leaves;
      continue;
    }
    int j = static_cast<int>(std::max_element(widths.begin(),widths.end())-widths.begin());
    assert(widths[j]>=4);
    int mid = (box[j].lo+box[j].hi)/2;
    Box right=box;
    right[j].lo=mid;
    box[j].hi=mid;
    pending.push_back(right);
    pending.push_back(box);
  }
  std::cout << "两种匹配的四维相位环面，认证闭盒数 " << leaves << '\n';
}
