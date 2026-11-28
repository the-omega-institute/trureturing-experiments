// 在每条线的位置确定后，扣除它从逐个不等式的最大值落到这个位置的差额。
#include <cstddef>
#include <cstdint>
#include <vector>

using Integer = std::int64_t;

static bool choose_cells(const Integer* reductions, const int* lines,
                         int depth, int segments, const std::vector<Integer>& budget) {
    if (depth == 4) return true;
    auto count = budget.size();
    std::vector<Integer> next(count);
    auto start = reductions + std::size_t(lines[depth])*segments*count;
    for (int cell=0; cell<segments; ++cell) {
        auto loss = start + std::size_t(cell)*count;
        bool possible = true;
        for (std::size_t r=0; r<count; ++r) {
            next[r] = budget[r] + loss[r];
            if (next[r] < 0) {
                possible = false;
                break;
            }
        }
        if (possible && choose_cells(reductions,lines,depth+1,segments,next)) return true;
    }
    return false;
}

extern "C" void scan(const Integer* values, const Integer* maxima, const Integer* reductions,
                     const Integer* totals, const int* candidates, int groups, int fields,
                     int segments, int* keep) {
    std::vector<Integer> budget(fields);
    for (int group=0; group<groups; ++group) {
        const auto* lines = candidates + 5*group;
        int shared = lines[4];
        for (int r=0; r<fields; ++r) {
            Integer available = totals[shared*fields+r];
            for (int k=0; k<4; ++k) {
                auto i = lines[k];
                available += maxima[i*fields+r] - values[(i*segments+shared)*fields+r];
            }
            budget[r] = available;
        }
        keep[group] = choose_cells(reductions,lines,0,segments,budget);
    }
}
