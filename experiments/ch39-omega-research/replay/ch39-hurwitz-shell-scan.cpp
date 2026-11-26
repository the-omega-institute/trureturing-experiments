// Enumerate in-box primitive norm-five rotors in signed 64-bit arithmetic.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using Int = std::int64_t;
using Wide = __int128_t;
constexpr std::array<Int,3> center = {467453,638110,1058300};
constexpr Int denominator = 1000000000;

Int root(Int x) {
    Int r = static_cast<Int>(std::sqrt(static_cast<double>(x)));
    while (r*r > x) --r;
    while ((r+1)*(r+1) <= x) ++r;
    return r;
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    int depth = std::atoi(argv[1]);
    if (depth < 0 || depth > 24) return 1;
    Int norm = 1;
    for (int j=0;j<depth;++j) norm *= 5;
    Wide lower = Wide(denominator)*denominator, upper = lower;
    for (Int c : center) {
        lower += Wide(c*999)*(c*999);
        upper += Wide(c*1001)*(c*1001);
    }
    const Wide scaled = Wide(norm)*denominator*denominator;
    const Int first = root(static_cast<Int>(scaled/upper));
    const Int last = root(static_cast<Int>(scaled/lower));
    constexpr int period = 3*5*7*11;
    std::array<bool,period> square{};
    for (int i=0;i<period;++i) square[i*i%period] = true;
    std::vector<std::vector<int>> admissible(period);
    for (int n=0;n<period;++n) {
        for (int t=0;t<period;++t) {
            int r = (n+period-(t*t%period))%period;
            if (square[r]) admissible[n].push_back(t);
        }
    }
    for (Int a=first+(first%2==0);a<=last;a+=2) {
        if (Wide(a)*a*upper < scaled) continue;
        std::array<Int,3> lo,hi;
        for (int j=0;j<3;++j) {
            lo[j] = (a*center[j]*999+denominator-1)/denominator;
            hi[j] = a*center[j]*1001/denominator;
            lo[j] += lo[j]%2;
        }
        for (Int b=lo[0];b<=hi[0];b+=2) {
            const Int tail = norm-a*a-b*b;
            if (tail < lo[2]*lo[2]) continue;
            Int c_lo = std::max(lo[1],root(std::max<Int>(0,tail-hi[2]*hi[2])));
            Int c_hi = std::min(hi[1],root(tail-lo[2]*lo[2]));
            // Set c=2t, d=2s. For every retained residue test the
            // complete integer square; the modular test is only a sieve.
            const Int t_lo = (c_lo+1)/2, t_hi = c_hi/2;
            const auto& residues = admissible[(tail/4)%period];
            for (Int start = (t_lo/period)*period;start<=t_hi;start+=period) {
                auto it = std::lower_bound(residues.begin(),residues.end(),t_lo-start);
                for (;it!=residues.end() && start+*it<=t_hi;++it) {
                    Int c = 2*(start+*it);
                    Int d2 = tail-c*c;
                    if (d2%16!=0 && d2%16!=4) continue;
                    Int d = root(d2);
                    if (d < lo[2] || d > hi[2] || d*d!=d2) continue;
                    if (a%5==0 && b%5==0 && c%5==0 && d%5==0) continue;
                    std::cout << a << ' ' << b << ' ' << c << ' ' << d << '\n';
                }
            }
        }
    }
}
