#include <array>
#include <cstdint>
#include <iostream>

using Integer = std::int64_t;
using Wide = __int128_t;
using Row = std::array<Integer, 3>;
using Ray = std::array<Wide, 3>;

constexpr int vertex_count = 1136;
constexpr int comparison_count = 8;

Ray cross(const Row &p, const Row &q) {
    return {Wide(p[1])*q[2]-Wide(p[2])*q[1], Wide(p[2])*q[0]-Wide(p[0])*q[2],
            Wide(p[0])*q[1]-Wide(p[1])*q[0]};
}

bool feasible(const std::array<Row, comparison_count+3> &constraints) {
    // w_i >= 0 与八个比较切出的归一重心多边形若非空，必有顶点。
    // 每个顶点是两个独立约束平面相交的一条非负射线。
    for (int a=0; a<comparison_count+3; ++a) {
        for (int b=0; b<a; ++b) {
            Ray w = cross(constraints[a], constraints[b]);
            if (w[0] <= 0 && w[1] <= 0 && w[2] <= 0) {
                for (Wide &weight : w) weight = -weight;
            }
            if (w[0] < 0 || w[1] < 0 || w[2] < 0 || (w[0] | w[1] | w[2]) == 0) continue;
            bool allowed = true;
            for (int r=0; r<comparison_count; ++r) {
                const auto &line = constraints[r];
                if (w[0]*line[0]+w[1]*line[1]+w[2]*line[2] > 0) {
                    allowed = false;
                    break;
                }
            }
            if (allowed) return true;
        }
    }
    return false;
}

int main() {
    std::array<std::array<Integer, comparison_count>, vertex_count> data;
    std::array<unsigned, vertex_count> positive{};
    for (int i=0; i<vertex_count; ++i) {
        for (int r=0; r<comparison_count; ++r) {
            Integer x;
            if (!(std::cin >> x) || x < -1000000000001 || x > 1000000000001) return 2;
            data[i][r] = x;
            if (x > 0) positive[i] |= 1u << r;
        }
    }
    // 所有检验量的绝对值至多为 6 (10^12+1)^3，严格小于 2^127。
    std::array<Row, comparison_count+3> constraints{};
    for (int j=0; j<3; ++j) constraints[comparison_count+j][j] = -1;
    std::uint64_t directly=0, by_intersections=0;
    for (int i=0; i<vertex_count; ++i) {
        for (int j=0; j<i; ++j) {
            unsigned mask = positive[i] & positive[j];
            for (int k=0; k<j; ++k) {
                if (mask & positive[k]) {
                    ++directly;
                    continue;
                }
                for (int r=0; r<comparison_count; ++r) {
                    constraints[r] = {data[i][r], data[j][r], data[k][r]};
                }
                if (feasible(constraints)) {
                    std::cerr << "保留的三顶点 " << i << ' ' << j << ' ' << k << '\n';
                    return 1;
                }
                ++by_intersections;
            }
        }
    }
    std::cout << "同一比较排除 " << directly << " 个三角形，整数交点排除 "
              << by_intersections << " 个，合计 " << directly+by_intersections << '\n';
    return directly+by_intersections == std::uint64_t(vertex_count)*(vertex_count-1)*(vertex_count-2)/6 ? 0 : 3;
}
