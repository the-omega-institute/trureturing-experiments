// Exact prefix comparisons for the least-modulus-nine branch.
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
constexpr int clipping_with_twenty_one = 17;
constexpr int clipping_without_twenty_one = 3;

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

int clipping_for_prime(int prime) {
    const auto band = std::lower_bound(schedule.begin(), schedule.end(), prime,
        [](const auto& row, int p) { return row[0] < p; });
    return (*band)[1];
}

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

Comparison extend(const Comparison& continuation, int prime, int clipping, bool omit_ternary_slot,
                  bool omit_square = false, Integer charge = scale) {
    require(!omit_square || prime == 5 || prime == 7, "unexpected omitted square");
    const int divisor = omit_square ? prime * prime : (prime <= 11 ? prime : 1);
    const int density_denominator = prime <= 7 ? divisor * (prime - 1) - 1 : divisor * (prime - 2);
    const int intercept = divisor + (omit_ternary_slot ? divisor / prime * (prime - 1) : 0);
    const int threshold = 100 * intercept + clipping * density_denominator;
    const Integer cap_numerator = 100 * divisor * (prime - 1);
    const Integer cap_denominator = (100 - clipping) * density_denominator;
    require(!omit_ternary_slot || prime == 5 || prime == 7 || prime == 11,
            "unexpected omitted slot");
    require(cap_numerator <= prime * cap_denominator, "invalid comparison cap");
    require(threshold <= 100 * divisor * width, "low-state table misses clipping threshold");

    const Integer previous = prime - 1;
    Comparison next{
        continuation.quadratic + ceiling(
            continuation.quadratic * cap_numerator * (3 * prime - 1),
            cap_denominator * previous * previous),
        continuation.linear + ceiling(
            (continuation.linear + charge) * 100 * divisor, cap_denominator),
        continuation.constant - charge * threshold / cap_denominator,
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
            const Integer deduction = charge * (100 * local_denominator - cap_denominator)
                / (100 * local_denominator);
            future = std::max(smaller_cap, future - deduction);
        }
        next.values[state] = future + ceiling(charge * std::max(0, above), cap_denominator);
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

Integer ternary_average(const Comparison& zero, const Comparison& positive,
                       int cap_numerator = 6, int cap_denominator = 5) {
    return expectation(zero.values[1], positive_expectation(positive, 3, 1, weights(3)),
                       3, cap_numerator, cap_denominator);
}

Integer truncated_expectation(const std::vector<Integer>& future, int prime, int state,
                              const GeometricWeights& geometric, Integer numerator,
                              Integer denominator) {
    const Integer zero_denominator = prime * denominator;
    Integer total = (zero_denominator - numerator) * future[state] / zero_denominator;
    const int cutoff = std::min((int(future.size()) - 1) / state, int(geometric.powers.size()) - 1);
    for (int n = 2; n <= cutoff; ++n) {
        total += numerator * (prime - 1) * future[state * n] / (denominator * geometric.powers[n]);
    }
    return total;
}

void extend_floor(std::vector<Integer>& future, int prime, int clipping,
                  bool omit_square, Integer charge, bool both_branches) {
    const int divisor = omit_square ? prime * prime : (prime <= 7 ? prime : 1);
    const int density_denominator = prime <= 7 ? divisor * (prime - 1) - 1 : prime - 2;
    const int capped_denominator = (100 - clipping) * density_denominator;
    const int base_numerator = divisor * (prime - 1);
    const auto geometric = weights(prime);
    std::vector<Integer> next(future.size());
    for (int state = 1; state < int(future.size()); ++state) {
        const int excess = divisor * (state - 1);
        const int above = 100 * excess - clipping * density_denominator;
        const int denominator = above < 0 ? density_denominator - excess : capped_denominator;
        const int numerator = above < 0 ? base_numerator : 100 * base_numerator;
        Integer total = truncated_expectation(future, prime, state, geometric, numerator, denominator);
        if (both_branches && above < 0) {
            const auto capped = truncated_expectation(
                future, prime, state, geometric, 100 * base_numerator, capped_denominator);
            const auto deduction = ceiling(charge * (100 * denominator - capped_denominator),
                                           100 * denominator);
            total = std::max(total, capped - deduction);
        }
        next[state] = total + charge * std::max(0, above) / capped_denominator;
    }
    future.swap(next);
}

std::vector<Integer> first_branch_floor(const std::vector<int>& primes,
                                       int first_prime, int last_prime, int state_limit,
                                       int clipping_at_seven) {
    // Use the first branch only, with a zero terminal function and bounded states.
    std::vector<int> suffix_primes;
    if (first_prime == 5) {
        suffix_primes.push_back(5);
    }
    suffix_primes.push_back(7);
    suffix_primes.insert(suffix_primes.end(), primes.begin(),
                         std::upper_bound(primes.begin(), primes.end(), last_prime));
    std::vector<Integer> future(state_limit + 1);
    for (auto it = suffix_primes.rbegin(); it != suffix_primes.rend(); ++it) {
        const int prime = *it;
        const int clipping = prime == 7 ? clipping_at_seven : clipping_for_prime(prime);
        extend_floor(future, prime, clipping, false, scale, false);
    }
    return future;
}

Integer ternary_floor(const std::vector<Integer>& future) {
    const auto geometric = weights(3);
    const int cutoff = std::min(int(future.size()) - 1, int(geometric.powers.size()) - 1);
    Integer result = 3 * future[1] / 5;
    for (int n = 2; n <= cutoff; ++n) {
        result += 12 * future[n] / (5 * geometric.powers[n]);
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
        for (auto p = primes.rbegin(); p != primes.rend() && *p > 11; ++p) {
            continuation = extend(continuation, *p, clipping_for_prime(*p), false);
        }
        const auto eleven_without_thirty_three = extend(continuation, 11, clipping_for_prime(11), true);
        continuation = extend(continuation, 11, clipping_for_prime(11), false);
        const Comparison before_seven = extend(continuation, 7, clipping_with_twenty_one, false);
        const Comparison zero = extend(before_seven, 5, 0, false);
        const Comparison positive = extend(before_seven, 5, 0, true);
        const Integer total = ternary_average(zero, positive);
        require(0 < total && total * 500 < Integer(scale) * 499,
                "reserve does not exceed 1/500");
        const auto seven_zero = extend(continuation, 7, clipping_without_twenty_one, false);
        const auto seven_positive = extend(continuation, 7, clipping_without_twenty_one, true);
        const auto aligned_zero = extend(seven_zero, 5, 0, false);
        const auto aligned_positive = extend(seven_positive, 5, 0, true);
        const Integer missing_both = ternary_average(aligned_zero, aligned_positive);
        const auto seven_without_thirty_three = extend(
            eleven_without_thirty_three, 7, clipping_without_twenty_one, true);
        const auto before_five_without_thirty_three = extend(
            seven_without_thirty_three, 5, 0, true);
        const Integer missing_three = ternary_average(aligned_zero, before_five_without_thirty_three);
        const Integer lower = first_branch_floor(primes, 7, 211, 400, clipping_without_twenty_one)[1];
        const Integer zero_tail = ternary_floor(
            first_branch_floor(primes, 5, 1999, 800, clipping_with_twenty_one));
        constexpr Integer charge_at_five_without_twenty_one = 97 * Integer(scale) / 100;
        const auto five_all_zero = extend(
            seven_zero, 5, 0, false, false, charge_at_five_without_twenty_one);
        const auto five_all_positive = extend(
            seven_positive, 5, 0, false, false, charge_at_five_without_twenty_one);
        const Integer limited_ternary_union = ternary_average(five_all_zero, five_all_positive, 36, 31);
        const auto five_no_square_zero = extend(
            seven_zero, 5, 0, false, true, charge_at_five_without_twenty_one);
        const auto five_no_square_positive = extend(
            seven_positive, 5, 0, false, true, charge_at_five_without_twenty_one);
        const Integer no_twenty_one_or_square = ternary_average(
            five_no_square_zero, five_no_square_positive);
        std::vector<Integer> floor_with_terminal(1201);
        for (int state = 1; state < int(floor_with_terminal.size()); ++state) {
            floor_with_terminal[state] = Integer(tail) * state * state;
        }
        for (auto p = primes.rbegin(); p != primes.rend(); ++p) {
            extend_floor(floor_with_terminal, *p, clipping_for_prime(*p), false, scale, true);
        }
        constexpr Integer credit_at_seven = 81 * Integer(scale) / 10'000;
        require(floor_with_terminal[1] >= credit_at_seven, "prime-7 credit exceeds suffix floor");
        auto floor_with_zero_clipping = floor_with_terminal;
        extend_floor(floor_with_zero_clipping, 7, 0, false, scale - credit_at_seven, true);
        extend_floor(floor_with_terminal, 7, clipping_with_twenty_one, true,
                     scale - credit_at_seven, true);
        constexpr Integer credit_at_five = 3486 * Integer(scale) / 100'000;
        require(floor_with_terminal[1] >= credit_at_five, "prime-5 credit exceeds suffix floor");
        require(floor_with_zero_clipping[1] >= credit_at_five,
                "prime-5 credit exceeds suffix floor at zero clipping");
        const auto seven_no_square = extend(
            continuation, 7, clipping_with_twenty_one, false, true, scale - credit_at_seven);
        const auto five_no_square = extend(
            seven_no_square, 5, 0, false, true, scale - credit_at_five);
        const Integer no_small_powers = ternary_average(five_no_square, five_no_square, 54, 47);
        const auto unrestricted_at_seven = extend(
            continuation, 7, 0, false, false, scale - credit_at_seven);
        const auto unrestricted_at_five = extend(
            unrestricted_at_seven, 5, 0, false, false, scale - credit_at_five);
        const Integer early_mixed_bound = ternary_average(
            unrestricted_at_five, unrestricted_at_five);
        require(0 < missing_three && missing_three < missing_both && missing_both < scale
                && 0 <= lower && lower < scale,
                "invalid phase comparison bounds");
        std::cout << primes.size() << ' ' << std::int64_t(total) << ' '
                  << std::int64_t(missing_both) << ' ' << std::int64_t(missing_three) << ' '
                  << std::int64_t(lower) << ' ' << std::int64_t(zero_tail) << ' '
                  << std::int64_t(limited_ternary_union) << ' ' << std::int64_t(no_twenty_one_or_square)
                  << ' ' << std::int64_t(no_small_powers) << ' '
                  << std::int64_t(early_mixed_bound) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
