#include "/Users/macstudio/omega/trureturing-01a12631-certificate/.ejm-proof-data/maps.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// A read-only producer diagnostic. This is not Lean acceptance evidence.
using Z = int64_t;
using Children = std::array<uint32_t, 4>;
const std::string source =
    "/Users/macstudio/omega/trureturing-01a12631-certificate/.ejm-proof-data/rounded32/";
std::vector<Children> nodes;
std::vector<bool> cutfree;
std::vector<Z> tensor[2];
int hashes[32];
int weights[20][32]{};
uint64_t splits = 0, freeSplits = 0, positiveSplits = 0, positiveFreeSplits = 0;
Z splitMax = std::numeric_limits<Z>::min();
Z freeMax = std::numeric_limits<Z>::min();
Z freeMin = std::numeric_limits<Z>::max();
uint64_t perDepth[21]{}, positiveDepth[21]{};

template <class T> std::vector<T> read(const std::string &name, size_t count) {
  std::ifstream in(source + name, std::ios::binary);
  std::vector<T> result(count);
  in.read(reinterpret_cast<char *>(result.data()), count * sizeof(T));
  if (!in || in.peek() != EOF) throw std::runtime_error("Bad input: " + name);
  return result;
}

void scan(uint32_t index, int depth) {
  if (index < 92161) return;
  ++splits;
  Z lower = 0;
  for (int k = 0; k < 32; ++k) lower += tensor[k / 16][hashes[k]];
  splitMax = std::max(splitMax, lower);
  positiveSplits += lower >= 1;
  if (cutfree[index]) {
    ++freeSplits;
    ++perDepth[depth];
    freeMax = std::max(freeMax, lower);
    freeMin = std::min(freeMin, lower);
    positiveFreeSplits += lower >= 1;
    positiveDepth[depth] += lower >= 1;
  }
  if (depth >= 20) throw std::runtime_error("Split beyond actual coordinate domain");
  for (int a = 0; a < 4; ++a) {
    for (int k = 0; k < 32; ++k) hashes[k] += (a - 4) * weights[depth][k];
    scan(nodes[index - 92161][a], depth + 1);
    for (int k = 0; k < 32; ++k) hashes[k] -= (a - 4) * weights[depth][k];
  }
}

int main() {
  const auto started = std::chrono::steady_clock::now();
  std::ifstream in(source + "cover.dag.bin", std::ios::binary);
  uint32_t root, count;
  in.read(reinterpret_cast<char *>(&root), 4);
  in.read(reinterpret_cast<char *>(&count), 4);
  nodes.resize(count);
  in.read(reinterpret_cast<char *>(nodes.data()), count * sizeof(Children));
  if (!in || in.peek() != EOF || root != 1646020 || count != 1553860)
    throw std::runtime_error("Bad literal DAG");
  cutfree.resize(count + 92161);
  cutfree[0] = true;
  for (uint32_t j = 0; j < count; ++j) {
    bool free = true;
    for (auto child : nodes[j]) {
      if (child >= j + 92161) throw std::runtime_error("Non-topological DAG");
      free &= cutfree[child];
    }
    cutfree[j + 92161] = free;
  }
  int pow5[11];
  pow5[0] = 1;
  for (int i = 1; i <= 10; ++i) pow5[i] = 5 * pow5[i - 1];
  tensor[0] = read<Z>("centered-0.int64le", pow5[10]);
  tensor[1] = read<Z>("centered-1.int64le", pow5[9]);
  for (int block = 0; block < 2; ++block) {
    int n = block == 0 ? 10 : 9;
    for (int place = 0; place < 16; ++place) {
      hashes[16 * block + place] = pow5[n] - 1;
      for (int i = 0; i < n; ++i) {
        int g = block == 0 ? maps0[place][i] : maps1[place][i];
        weights[g][16 * block + place] = pow5[i];
      }
    }
  }
  scan(root, 0);
  auto seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  std::cout << "producer_only split_occurrences=" << splits
            << " split_max_lower=" << splitMax
            << " positive_split_occurrences=" << positiveSplits
            << " cutfree_split_occurrences=" << freeSplits
            << " cutfree_min_lower=" << freeMin
            << " cutfree_max_lower=" << freeMax
            << " positive_cutfree_occurrences=" << positiveFreeSplits
            << " wall_seconds=" << seconds << '\n';
  for (int d = 0; d <= 20; ++d)
    if (perDepth[d]) std::cout << "depth=" << d << " cutfree=" << perDepth[d]
                             << " positive=" << positiveDepth[d] << '\n';
}
