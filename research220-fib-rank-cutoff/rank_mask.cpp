#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

static std::int64_t mod_pow(std::int64_t base, std::int64_t exponent,
                            std::int64_t modulus) {
  std::int64_t result = 1 % modulus;
  while (exponent != 0) {
    if ((exponent & 1) != 0) {
      result = static_cast<u64>(result) * base % modulus;
    }
    base = static_cast<u64>(base) * base % modulus;
    exponent >>= 1;
  }
  return result;
}

static std::pair<std::int64_t, std::int64_t> fib_pair(std::int64_t n,
                                                       std::int64_t modulus) {
  if (n == 0) return {0, 1};
  const auto [a, b] = fib_pair(n >> 1, modulus);
  const auto c = static_cast<u64>(a) * ((2 * b % modulus - a + modulus) % modulus) % modulus;
  const auto d = (static_cast<u64>(a) * a + static_cast<u64>(b) * b) % modulus;
  if ((n & 1) != 0) return {d, (c + d) % modulus};
  return {c, d};
}

int main(int argc, char **argv) {
  const int cutoff = argc > 1 ? std::atoi(argv[1]) : 7'608'793;
  std::vector<bool> prime(cutoff + 1, true);
  prime[0] = prime[1] = false;
  std::vector<int> primes;
  for (int n = 2; n <= cutoff; ++n) {
    if (!prime[n]) continue;
    primes.push_back(n);
    if (static_cast<std::int64_t>(n) * n <= cutoff) {
      for (std::int64_t m = static_cast<std::int64_t>(n) * n; m <= cutoff; m += n)
        prime[static_cast<std::size_t>(m)] = false;
    }
  }

  std::map<int, std::pair<int, int>> rank_witness;
  for (const int p : primes) {
    if (p == 2 || p == 5) continue;
    const int legendre = mod_pow(5, (p - 1) / 2, p) == 1 ? 1 : -1;
    const int candidate = p - legendre;
    int remaining = candidate;
    std::vector<int> factors;
    for (const int q : primes) {
      if (static_cast<std::int64_t>(q) * q > remaining) break;
      if (remaining % q != 0) continue;
      factors.push_back(q);
      while (remaining % q == 0) remaining /= q;
    }
    if (remaining > 1) factors.push_back(remaining);

    int rank = candidate;
    for (const int q : factors) {
      while (rank % q == 0 && fib_pair(rank / q, p).first == 0) rank /= q;
    }
    if (rank >= 7 && rank <= cutoff && prime[rank]) {
      auto &entry = rank_witness[rank];
      if (entry.second == 0 || p < entry.first) entry.first = p;
      ++entry.second;
    }
  }

  std::cout << "prime_rank_count\t" << rank_witness.size() << "\n";
  for (const auto &[rank, witness] : rank_witness)
    std::cout << rank << '\t' << witness.first << '\t' << witness.second << "\n";
}
