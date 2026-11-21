// Exact upward-rounded prefix comparison with modulus 15 omitted.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using Integer = __int128_t;

constexpr std::int64_t scale = 100'000'000'000'000;
constexpr int endpoint = 80'021;
constexpr int width = 34'500;

// Right endpoint of each interval, with its clipping parameter in hundredths.
constexpr std::array schedule{
    std::array{5, 0}, std::array{7, 17}, std::array{11, 11},
    std::array{13, 18}, std::array{17, 20}, std::array{19, 18},
    std::array{23, 24}, std::array{29, 26}, std::array{31, 24},
    std::array{37, 26}, std::array{41, 28}, std::array{43, 27},
    std::array{47, 31}, std::array{53, 29}, std::array{101, 30},
    std::array{211, 32}, std::array{419, 34}, std::array{839, 35},
    std::array{1679, 37}, std::array{3359, 38}, std::array{16001, 38},
    std::array{endpoint, 43},
};

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

Integer ceiling(Integer numerator, Integer denominator) {
    return numerator / denominator + (numerator % denominator > 0);
}

struct Comparison {
    // A*x^2 + B*x + C is an upper bound above width.
    // A*x^2 + B*x is also an upper bound at every positive state.
    Integer quadratic;
    Integer linear;
    Integer constant;
    std::vector<Integer> values;
};

struct GeometricWeights {
    std::vector<Integer> powers;
    std::vector<Integer> point;
};

GeometricWeights weights(int prime) {
    GeometricWeights result{{1, prime}, {0, 0}};
    while (result.powers.back() < scale) {
        result.powers.push_back(result.powers.back() * prime);
        result.point.push_back(ceiling(Integer(scale) * (prime - 1), result.powers.back()));
    }
    return result;
}

Integer positive_expectation(const Comparison& continuation, int prime, int state,
                             const GeometricWeights& geometric) {
    const int tabulated = width / state;
    const int cutoff = std::min(tabulated, int(geometric.powers.size()) - 1);
    const Integer previous = prime - 1;
    const Integer shifted_mean = cutoff * previous + prime;
    // If G is the positive geometric variable, Var(G) = p/(p-1)^2.
    const Integer second = shifted_mean * shifted_mean + prime;
    const Integer constant = cutoff == tabulated ? continuation.constant : 0;
    const Integer tail = continuation.quadratic * state * state * second
        + continuation.linear * state * shifted_mean * previous
        + constant * previous * previous;
    require(tail >= 0, "comparison polynomial has negative tail expectation");
    Integer points = 0;
    for (int multiplier = 2; multiplier <= cutoff; ++multiplier) {
        points += geometric.point[multiplier] * continuation.values[state * multiplier];
    }
    // Both ceilings are upper bounds in the units of the state table.
    return ceiling(points, scale)
        + ceiling(tail, geometric.powers[cutoff] * previous * previous);
}

Integer expectation(Integer fixed, Integer positive, int prime,
                    Integer cap_numerator, Integer cap_denominator) {
    return fixed + ceiling(
        cap_numerator * (prime * positive - fixed), prime * cap_denominator);
}

Comparison extend(const Comparison& continuation, int prime, int clipping, bool skip) {
    const int divisor = prime <= 7 ? prime : 1;
    const int density_denominator = prime <= 7 ? prime * (prime - 1) - 1 : prime - 2;
    const int intercept = divisor + (skip ? prime - 1 : 0);
    const int threshold = 100 * intercept + clipping * density_denominator;
    const Integer cap_numerator = 100 * divisor * (prime - 1);
    const Integer cap_denominator = (100 - clipping) * density_denominator;
    require(!skip || prime == 5, "unexpected omitted slot");
    require(cap_numerator <= prime * cap_denominator, "invalid comparison cap");
    require(threshold <= 100 * divisor * width, "low-state table misses clipping threshold");

    const Integer previous = prime - 1;
    Comparison next{
        continuation.quadratic + ceiling(
            continuation.quadratic * cap_numerator * (3 * prime - 1),
            cap_denominator * previous * previous),
        continuation.linear + ceiling(
            (continuation.linear + scale) * 100 * divisor, cap_denominator),
        continuation.constant - Integer(scale) * threshold / cap_denominator,
        std::vector<Integer>(width + 1),
    };
    const auto geometric = weights(prime);
    for (int state = 1; state <= width; ++state) {
        const Integer fixed = continuation.values[state];
        const Integer positive = positive_expectation(continuation, prime, state, geometric);
        Integer future = expectation(fixed, positive, prime, cap_numerator, cap_denominator);
        const int above = 100 * divisor * state - threshold;
        if (above < 0) {
            const int local_denominator = density_denominator
                - std::max(0, state * divisor - intercept);
            const Integer smaller_cap = expectation(
                fixed, positive, prime, divisor * (prime - 1), local_denominator);
            // Subtract a lower rounding of (bar c - c(a))/bar c.
            const Integer deduction = Integer(scale) * (100 * local_denominator - cap_denominator)
                / (100 * local_denominator);
            future = std::max(smaller_cap, future - deduction);
        }
        next.values[state] = future + ceiling(Integer(scale) * std::max(0, above), cap_denominator);
    }
    return next;
}

std::vector<int> primes_through_endpoint() {
    std::vector<bool> prime(endpoint + 1, true);
    std::vector<int> result;
    for (int p = 2; p <= endpoint; ++p) {
        if (prime[p]) {
            if (p >= 11) {
                result.push_back(p);
            }
            if (p <= endpoint / p) {
                for (int multiple = p * p; multiple <= endpoint; multiple += p) {
                    prime[multiple] = false;
                }
            }
        }
    }
    return result;
}

int main() {
    try {
        std::int64_t tail;
        require(bool(std::cin >> tail) && 0 < tail && tail < scale / endpoint,
                "invalid terminal upper bound");
        Comparison continuation{tail, 0, 0, std::vector<Integer>(width + 1)};
        for (int state = 1; state <= width; ++state) {
            continuation.values[state] = Integer(tail) * state * state;
        }
        const auto primes = primes_through_endpoint();
        for (auto p = primes.rbegin(); p != primes.rend(); ++p) {
            const auto band = std::lower_bound(
                schedule.begin(), schedule.end(), *p,
                [](const auto& row, int prime) { return row[0] < prime; });
            continuation = extend(continuation, *p, (*band)[1], false);
        }
        const Comparison before_seven = extend(continuation, 7, 17, false);
        const Comparison zero = extend(before_seven, 5, 0, false);
        const Comparison positive = extend(before_seven, 5, 0, true);
        const Integer total = ceiling(
            3 * zero.values[1] + 6 * positive_expectation(positive, 3, 1, weights(3)), 5);
        require(0 < total && total * 500 < Integer(scale) * 499,
                "reserve does not exceed 1/500");
        std::cout << primes.size() << ' ' << std::int64_t(total) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
