/* 尺度 M=10^7。所有区间端点都向外取整，平方比较使用 128 位整数。 */
#include <stdint.h>
#include <stdlib.h>

typedef int64_t Integer;
typedef __int128 Wide;
static const Integer M = 10000000;

static Integer minimum(Integer x, Integer y) { return x < y ? x : y; }
static Integer maximum(Integer x, Integer y) { return x > y ? x : y; }

static Integer floor_div(Integer x, Integer y) {
    return x / y - (x % y < 0);
}

static Integer ceil_div(Integer x, Integer y) {
    return -floor_div(-x, y);
}

static Integer product_low(Integer a, Integer b, Integer c, Integer d) {
    return minimum(minimum(a*c, a*d), minimum(b*c, b*d));
}

static Integer product_high(Integer a, Integer b, Integer c, Integer d) {
    return maximum(maximum(a*c, a*d), maximum(b*c, b*d));
}

int supports(Integer count, Integer segments, const Integer *form,
             const Integer *spans, const Integer *q) {
    Integer pl=q[0], ph=q[1], ll=q[2], lh=q[3], td=q[4], rd=q[5], loss=q[6], sine=q[7];
    if (td <= 0 || count <= 0 || count > 80 || segments <= 0 || segments > 512) return 0;
    for (int k=0; k<8; k++) if (q[k] < -4*M || q[k] > 4*M) return 0;
    Integer base[80], parallel_low[80], parallel_high[80], norm_squared[80];
    for (int i=0; i<count; i++) {
        const Integer *f = form+5*i;
        base[i] = f[0] - ceil_div(f[1]*M, td);
        parallel_low[i] = product_low(pl, ph, f[2], f[3]);
        parallel_high[i] = product_high(pl, ph, f[2], f[3]);
        norm_squared[i] = (Integer)(((Wide)f[4]*rd*rd + M-1)/M);
    }
    for (int j=0; j<segments; j++) {
        int separated=0;
        for (int i=0; i<count; i++) {
            const Integer *v=spans+2*(count*j+i);
            Integer low=floor_div(parallel_low[i]+product_low(ll,lh,v[0],v[1]),M);
            Integer high=ceil_div(parallel_high[i]+product_high(ll,lh,v[0],v[1]),M);
            Integer abs_high=maximum(llabs(low),llabs(high));
            Integer abs_low=maximum(0,maximum(low,-high));
            Integer margin=base[i]+low-ceil_div(abs_high*loss,M);
            if (margin <= 0) continue;
            Integer tangent_squared=maximum(0,norm_squared[i]-abs_low*abs_low);
            if ((Wide)margin*margin*M*M > (Wide)tangent_squared*sine*sine) {
                separated=1;
                break;
            }
        }
        if (!separated) return 0;
    }
    return 1;
}
