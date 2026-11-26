// Finite torus sum for the one-loop normal, using pointwise spectral bounds.
#include <algorithm>
#include <cmath>
#include <cstdint>

extern "C" void project_nodes(int count, const std::uint16_t *angles,
                              const double *weyl, const double *table,
                              double *sums) {
    double positive = 0, negative = 0, weak_upper = 0;
    constexpr int charges[15] = {1,1,-4,2,1,1,-4,2,1,1,-4,2,-3,-3,6};
    constexpr int weak[15] = {1,-1,0,0,1,-1,0,0,1,-1,0,0,1,-1,0};
    constexpr int color[4] = {1,1,-1,-1};
    for (int i = 0; i < count; ++i) {
        double low = weyl[i], high = weyl[i], y = 0, w = 0, dy = 0, dw = 0;
        double c[3] = {}, dc[3] = {};
        for (int j = 0; j < 15; ++j) {
            const double *entry = table + 4*angles[15*i+j];
            low *= entry[0];
            high *= entry[1];
            y += charges[j]*entry[2];
            dy += std::abs(charges[j])*entry[3];
            w += weak[j]*entry[2];
            dw += std::abs(weak[j])*entry[3];
            if (j < 12) {
                c[j/4] += color[j%4]*entry[2];
                dc[j/4] += entry[3];
            }
        }
        // The margin exceeds roundoff in all fifteen slopes and their sums.
        auto floor_square = [](double value, double radius) {
            const double x = std::max(0.0, std::abs(value)-radius-5e-7);
            return x*x;
        };
        double gain = 0.536*floor_square(y,dy)/60;
        for (int a = 0; a < 3; ++a) {
            const int b = (a+1)%3;
            gain += floor_square(c[a]-c[b], dc[a]+dc[b])/48;
        }
        const double max_w = std::abs(w)+dw+5e-7;
        const double w2 = max_w*max_w/12;
        const double loss = 1.516*w2;
        const double residual = gain-loss-1e-11*(gain+loss);
        positive += low*std::max(0.0,residual);
        negative += high*std::max(0.0,-residual);
        weak_upper += high*w2;
    }
    sums[0] = positive;
    sums[1] = negative;
    sums[2] = weak_upper;
}
