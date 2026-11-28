#include <algorithm>
#include <cstddef>
#include <cstdint>

#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace {
constexpr int columns = 15;
constexpr int rows = 16;
constexpr int control_count = 16;
constexpr int parts = 2;
constexpr std::int64_t shift = 1 << 15;
constexpr std::int64_t unit = std::int64_t{1} << 38;

// Each operand is a bounded signed 32-bit integer. Even without cancellation,
// 32 products bounded by 2^55 and one unit term fit in signed 64-bit arithmetic.
void complex_dot(const std::int32_t* left, const std::int32_t* right,
                 std::int64_t& real, std::int64_t& imag) {
#if defined(__aarch64__)
    int64x2_t real_sum = vdupq_n_s64(0);
    int64x2_t imag_sum = vdupq_n_s64(0);
    for (int index = 0; index < rows; index += 4) {
        const int32x4x2_t l = vld2q_s32(left + parts * index);
        const int32x4x2_t r = vld2q_s32(right + parts * index);
        real_sum = vmlal_s32(real_sum, vget_low_s32(l.val[0]), vget_low_s32(r.val[0]));
        real_sum = vmlsl_s32(real_sum, vget_low_s32(l.val[1]), vget_low_s32(r.val[1]));
        real_sum = vmlal_high_s32(real_sum, l.val[0], r.val[0]);
        real_sum = vmlsl_high_s32(real_sum, l.val[1], r.val[1]);
        imag_sum = vmlal_s32(imag_sum, vget_low_s32(l.val[0]), vget_low_s32(r.val[1]));
        imag_sum = vmlal_s32(imag_sum, vget_low_s32(l.val[1]), vget_low_s32(r.val[0]));
        imag_sum = vmlal_high_s32(imag_sum, l.val[0], r.val[1]);
        imag_sum = vmlal_high_s32(imag_sum, l.val[1], r.val[0]);
    }
    real = vaddvq_s64(real_sum);
    imag = vaddvq_s64(imag_sum);
#else
    real = 0;
    imag = 0;
    for (int index = 0; index < parts * rows; index += parts) {
        const std::int64_t lr = left[index];
        const std::int64_t li = left[index + 1];
        const std::int64_t rr = right[index];
        const std::int64_t ri = right[index + 1];
        real += lr * rr - li * ri;
        imag += lr * ri + li * rr;
    }
#endif
}

std::int64_t defect_units(std::int64_t residual) {
    const std::int64_t magnitude = residual < 0 ? -residual : residual;
    return (std::min(magnitude, unit) + shift - 1) / shift;
}
}

extern "C" void max_squared_defects(const std::int32_t* left,
                                    const std::int32_t* controls,
                                    std::int64_t* output,
                                    std::size_t count) {
    for (std::size_t case_index = 0; case_index < count; ++case_index) {
        const auto* case_left = left + case_index * columns * rows * parts;
        const auto* case_controls = controls + case_index * control_count * rows * columns * parts;
        std::int64_t maximum = 0;
        for (int control_index = 0; control_index < control_count; ++control_index) {
            const auto* matrix = case_controls + control_index * rows * columns * parts;
            alignas(64) std::int32_t transposed[columns][rows * parts];
            for (int row = 0; row < rows; ++row) {
                for (int column = 0; column < columns; ++column) {
                    const auto* entry = matrix + parts * (row * columns + column);
                    transposed[column][row * parts] = entry[0];
                    transposed[column][row * parts + 1] = entry[1];
                }
            }
            std::int64_t total = 0;
            for (int row = 0; row < columns; ++row) {
                for (int column = 0; column < columns; ++column) {
                    std::int64_t real, imag;
                    complex_dot(case_left + row * rows * parts, transposed[column], real, imag);
                    const auto re = defect_units(real - (row == column ? unit : 0));
                    const auto im = defect_units(imag);
                    total += re * re + im * im;
                }
            }
            maximum = std::max(maximum, total);
        }
        output[case_index] = maximum;
    }
}
