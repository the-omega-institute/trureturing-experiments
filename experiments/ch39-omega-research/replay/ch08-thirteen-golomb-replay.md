# Thirteen-mark modular Golomb ruler replays

A support of thirteen points in $\mathbb Z_d$ is admissible when its 156 directed differences between distinct points are nonzero and distinct. Equivalently, the 78 unoriented circular distances are distinct and exclude the half-turn in even dimension. This archive preserves the two C++20 searches used in the Appendix D classification. [Book research and witnesses](https://github.com/the-omega-institute/trureturing-books/issues/1676#issuecomment-7234902039)

## Minimum-gap search for $160\le d\le167$

`Search::run` rotates a smallest adjacent gap to start at zero and enumerates its possible lengths. The thirteen adjacent gaps have distinct circular distances, each at least the smallest gap. `Search::extend` adds points in order and rejects a repeated circular distance or a half-turn. `Search::enough_remaining_gap_length` sums the smallest still unused circular distances as a lower bound on the remaining arc length. These are necessary conditions on any completed support. [Original source and proof of pruning](https://github.com/the-omega-institute/trureturing-books/issues/1676#issuecomment-7234357783) · [Published exclusions](https://github.com/the-omega-institute/trureturing-books/issues/1676#issuecomment-7234582276)

## Unit-gap search for $169\le d\le182$ and $184\le d\le192$

In this dimension range the number of nonzero nonunits is below 156. An admissible support therefore has an ordered pair with a unit difference. Translation followed by multiplication by the inverse of that difference places the pair at 0 and 1. Reflection exchanging 0 and 1 puts the larger of the two adjacent gaps last. Equality would repeat a circular distance. `Search::run` and `Search::extend` enforce that normalization and carry the exact distance set of each candidate. `Search::enough_gap_length` bounds the remaining length using the smallest available distinct distances. Each actual arc gap is at least its corresponding circular distance. [Original enumerator and first runs](https://github.com/the-omega-institute/trureturing-books/issues/1676#issuecomment-7234835431) · [Later runs](https://github.com/the-omega-institute/trureturing-books/issues/1676#issuecomment-7234884921)

## Running the archived programs

Run from this directory with a C++20 compiler supporting `__uint128_t`.

```sh
c++ -std=c++20 -O3 ch08-thirteen-golomb-minimum-gap.cpp -o minimum-gap
c++ -std=c++20 -O3 ch08-thirteen-golomb-unit-gap.cpp -o unit-gap
./minimum-gap 160
./unit-gap 183
```

An independent rerun reproduced `160 UNSAT 2002756182`. The positive control at 183 returned `183 SAT 91786858 0 1 3 16 23 28 42 76 82 86 119 137 154`, with all 156 ordered differences distinct. The other archived node counts and witnesses are in the linked research records.
