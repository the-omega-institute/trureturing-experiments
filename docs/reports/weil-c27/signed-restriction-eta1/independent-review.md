# r2511: independent fixed-congruence source and rational review

Reviewed candidate `r2510/fixed_congruence.py`, SHA256 `d436cc16ffb63536a04d04ca2e586fdcca8af4b11a93e4f00de91f7c7180e24e`, and the complete host runner `host-bridge/certify-fixed-congruence-r2512.py`, SHA256 `eee6741b7dc32ef57d49577c090cbba90b0c684260a9e78a520f84cf0538c158`. No blocking mathematical or implementation error was found at these exact bytes. This is a source/helper review, not the actual112 sign certificate. Neither real112 data nor its planner was run by this reviewer; the root candidate was not modified.

## Exact mathematical contract

Let original T range over its same symmetric enclosed family. Define permutation P by P e_i=e_order[i], so Tp=P^T T P has entries T[order[i],order[j]]. The planner's row equilibration is Drow T Drow. Its approximate shears accumulate in U: each remaining column subtracts a linear combination of cached pivot columns, with multipliers from the current approximate Schur state. Reordering into removal order makes P^T U P upper triangular with diagonal exactly1. The final column powers give Dcol. The stored basis is exactly the quantized upper triangle of

\[
X=P^TD_{row}UD_{col}P.
\]

The numerical approximation need not satisfy an exact LDL identity. The only required exact properties are: all entries below the stored triangle are zero, and every stored diagonal is a strictly positive power of2. `enclose_congruence` checks those diagonal integers independently of the planner assertions. Thus X is invertible with determinant equal to the product of its nonzero diagonal entries. Dyadic rounding of off-diagonal entries does not change that conclusion. A nonrepresentable/zero diagonal is rejected rather than silently rounded into a singular basis.

The planner applies row scaling using the original row label, column scaling using the original column label, then permutes both. This is consistent with **first** making Tp and **then** applying the stored X. Interchanging these orders without transporting X would be incorrect; the runner uses the correct order.

## Product and inclusion review

For upper triangular X, the first product is

\[
right_{ij}=\sum_{k=0}^{j} (Tp)_{ik}X_{kj}=(TpX)_{ij}.
\]

The second is

\[
out_{ij}=\sum_{k=0}^{i}X_{ki}right_{kj}=(X^TTpX)_{ij}.
\]

Both ranges and transposes match the source. The implementation lifts original interval endpoints exactly by integer shifts, multiplies by fixed exact dyadic singleton intervals and adds outward. It computes only i<=j and mirrors the resulting interval, which encloses the same exact symmetric realization at both positions. Interval dependencies can make the enclosure wider but cannot invalidate inclusion. It does not clip, narrow or replace an input interval by its midpoint. Arbitrary planner errors or midpoints chosen from the same data remain sound because the subsequently fixed X is exact and invertible, and the entire sign proof uses outward intervals.

4096-bit multiplication limits additional rounding. It does not recover uncertainty already in the1024-bit source intervals. No successful final pivot or spectrum is predicted.

## Runner review

The reviewed runner pins the core, reciprocal source, corrected small reducer, congruence source and exact reduction bytes. It checks the stored basis grid4096, original grid1024 and dimension112; verifies order is a complete permutation and that flattening the proposed blocks equals that order; and calls `_check_blocks(..., complete=True)` before the real elimination.

It constructs `permuted=[[original[i][j] for j in order] for i in order]`. Its matrix positions are the transformed basis columns. Passing `labels=order` and the same original-label blocks names these columns by their planner provenance; the new result schema explicitly calls them transformed coordinates. Those names do not assert that a basis column is an untouched original coordinate vector. Subsequent `_eliminate` maps those names to the current positions and requires no particular pivot signs.

After elimination the runner requires active empty and the counted dimension112. It only then transports the full inertia back by Sylvester congruence and the inherited paid2/boundary19 elimination. A generic interval boundary produces a failure record with no full inertia. If the inherited source relation forces at least56 negative directions in T112, its sanity check correctly rejects a contradictory complete count. Source bindings, original physical allowance and existing premises are retained; no eta-budget, RH, compilation or official acceptance is inferred.

The actual result depends on the runner retaining normal Python assertions for its authentication/fresh-output checks; it is written for ordinary `python3`, not `python3 -O`. Core interval, symmetry, basis diagonal, coverage and pivot checks use explicit `require`. The reviewed execution path does not request optimized Python. No new user approval is required by this review.

## Independent synthetic verification

`check_fixed_congruence.py` ran once with actual process exit0. It imports only the four source modules; it never opens a scientific packet and does not call `propose_basis`. Its references use full Fraction matrix products with no triangular multiplication shortcuts, plus an exact coefficient hull for each entry of the symmetric interval family. Its inertia reference computes exact characteristic coefficients and uses real-root Descartes counts, independently of the interval Schur recurrence.

The six fixtures have dimensions1–5, nontrivial permutations, signed off-diagonal shears and separately varying row/column power-of-two scales. One uses the actual1024-to4096 grid change. One requires a prescribed two-coordinate indefinite pivot, exercising transformed-coordinate labels that differ from their positions. The checks comprise1250 exact product containment entries,64 exact box-hull entries,138 exhaustive symmetric-box vertices,6 independent inertia comparisons,6 complete interval counts,18 diagonal-scale checks,12 exact floor checks including negative values, and7 malformed/singular/grid/symmetry rejections. All passed. Native tool wall time was approximately0.117 seconds; the test's own measured elapsed time is in `synthetic-result.json`.

These tests verify the fixed-basis helper and its integration with the generic kernel on small independent operands. They do not certify the real planner choice, conditioning,112 inertia, restricted domination, tau budget or RH. Root may run the reviewed single actual112 computation under its current resource admission without repeating the large Gram.
