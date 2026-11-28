#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>

class Search {
    static constexpr int marks = 13;
    using Mask = __uint128_t;
    static int first_bit(Mask bits) {
        const auto low = static_cast<std::uint64_t>(bits);
        return low ? std::countr_zero(low) : 64 + std::countr_zero(static_cast<std::uint64_t>(bits >> 64));
    }
    int modulus;
    int maximum_distance;
    int first_gap = 0;
    std::array<int, marks> points{};
    std::array<std::array<Mask, 184>, 184> distance_bit{};
    std::uint64_t nodes = 0;

    bool enough_remaining_gap_length(int count, int last, Mask used) const {
        int gaps = marks - count + 1;
        int lower_bound = 0;
        const auto low = (Mask{1} << first_gap) - 1;
        const auto high = (Mask{1} << (maximum_distance + 1)) - 1;
        auto available = high & ~(used | low);
        while (available && gaps > 0) {
            lower_bound += first_bit(available);
            available &= available - 1;
            --gaps;
        }
        return gaps == 0 && lower_bound <= modulus - last;
    }

    bool extend(int count, Mask used) {
        ++nodes;
        const int last = points[count - 1];
        if (count == marks) {
            return true;
        }
        if (!enough_remaining_gap_length(count, last, used)) {
            return false;
        }
        const int limit = modulus - (marks - count) * first_gap;
        for (int next = last + first_gap; next <= limit; ++next) {
            auto next_used = used;
            bool valid = true;
            for (int j = 0; j < count; ++j) {
                const auto mask = distance_bit[points[j]][next];
                if (!mask || (next_used & mask)) {
                    valid = false;
                    break;
                }
                next_used |= mask;
            }
            if (!valid) {
                continue;
            }
            points[count] = next;
            if (extend(count + 1, next_used)) {
                return true;
            }
        }
        return false;
    }

public:
    explicit Search(int dimension) : modulus(dimension), maximum_distance((dimension - 1) / 2) {
        for (int x = 0; x < modulus; ++x) {
            for (int y = x + 1; y < modulus; ++y) {
                const int gap = y - x;
                const int distance = std::min(gap, modulus - gap);
                if (distance <= maximum_distance) {
                    distance_bit[x][y] = Mask{1} << distance;
                }
            }
        }
    }

    void run() {
        for (first_gap = 1; marks * first_gap + marks * (marks - 1) / 2 <= modulus; ++first_gap) {
            points[0] = 0;
            points[1] = first_gap;
            if (extend(2, Mask{1} << first_gap)) {
                std::cout << modulus << " SAT " << nodes;
                for (int point : points) {
                    std::cout << " " << point;
                }
                std::cout << std::endl;
                return;
            }
        }
        std::cout << modulus << " UNSAT " << nodes << std::endl;
    }
};

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }
    const auto dimension = std::atoi(argv[1]);
    if (dimension < 157 || dimension > 183) {
        return 2;
    }
    Search(dimension).run();
}
