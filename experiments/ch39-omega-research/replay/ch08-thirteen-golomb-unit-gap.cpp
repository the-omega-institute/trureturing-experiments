#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>

class Search {
    static constexpr int marks = 13;
    static constexpr int max_modulus = 212;
    using Mask = __uint128_t;

    struct Candidate {
        int position;
        Mask distances;
    };

    int modulus;
    int max_distance;
    Mask all_distances;
    std::array<int, marks> points{0, 1};
    std::array<std::array<Mask, max_modulus>, max_modulus> distance_bit{};
    std::uint64_t nodes = 0;

    static int first_available(Mask bits) {
        const auto low = static_cast<std::uint64_t>(bits);
        if (low != 0) {
            return std::countr_zero(low);
        }
        return 64 + std::countr_zero(static_cast<std::uint64_t>(bits >> 64));
    }

    bool enough_gap_length(int count, Mask used) const {
        int gaps_left = marks - count + 1;
        int distance_needed = 0;
        Mask available = all_distances & ~used;
        while (gaps_left > 0 && available != 0) {
            distance_needed += first_available(available);
            available &= available - 1;
            --gaps_left;
        }
        return gaps_left == 0 && distance_needed <= modulus - points[count - 1];
    }

    bool extend(int count, Mask used, const Candidate* candidates,
                int candidate_count, int last_limit) {
        ++nodes;
        if (count == marks) {
            return true;
        }
        if (candidate_count < marks - count || !enough_gap_length(count, used)) {
            return false;
        }
        std::array<Candidate, max_modulus> next_candidates;
        for (int i = 0; i <= candidate_count - (marks - count); ++i) {
            const auto [next, distances] = candidates[i];
            points[count] = next;
            const Mask next_used = used | distances;
            // Reflection exchanging 0 and 1 lets the closing gap be larger.
            const int limit = count == 2 ? modulus - next : last_limit;
            int next_count = 0;
            for (int j = i + 1; j < candidate_count; ++j) {
                const auto [position, earlier_distances] = candidates[j];
                if (position > limit) {
                    break;
                }
                const Mask mutual_distance = distance_bit[next][position];
                if (mutual_distance == 0 || (mutual_distance & earlier_distances) != 0) {
                    continue;
                }
                const Mask new_distances = earlier_distances | mutual_distance;
                if ((new_distances & next_used) == 0) {
                    next_candidates[next_count++] = {position, new_distances};
                }
            }
            if (extend(count + 1, next_used, next_candidates.data(), next_count, limit)) {
                return true;
            }
        }
        return false;
    }

public:
    explicit Search(int dimension) : modulus(dimension), max_distance((dimension - 1) / 2),
                                     all_distances((Mask{1} << (max_distance + 1)) - 2) {
        for (int a = 0; a < modulus; ++a) {
            for (int b = a + 1; b < modulus; ++b) {
                const int distance = std::min(b - a, modulus - (b - a));
                if (distance <= max_distance) {
                    distance_bit[a][b] = Mask{1} << distance;
                }
            }
        }
    }

    void run() {
        constexpr Mask fixed_difference = Mask{1} << 1;
        std::array<Candidate, max_modulus> candidates;
        int size = 0;
        for (int x = 2; x < modulus; ++x) {
            const Mask to_zero = distance_bit[0][x];
            const Mask to_one = distance_bit[1][x];
            if (to_zero && to_one && to_zero != to_one &&
                ((to_zero | to_one) & fixed_difference) == 0) {
                candidates[size++] = {x, to_zero | to_one};
            }
        }
        const bool found = extend(2, fixed_difference, candidates.data(), size, modulus - 1);
        std::cout << modulus << (found ? " SAT " : " UNSAT ") << nodes;
        if (found) {
            for (int point : points) {
                std::cout << ' ' << point;
            }
        }
        std::cout << '\n';
    }
};

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }
    const int dimension = std::atoi(argv[1]);
    if (dimension < 157 || dimension > 212) {
        return 2;
    }
    Search(dimension).run();
}
