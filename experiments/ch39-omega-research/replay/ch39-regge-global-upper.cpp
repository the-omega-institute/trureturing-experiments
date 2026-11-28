// Certify a local upper bound on the two normalized Freudenthal simplex charts.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <utility>

using Vector = std::array<double,10>;
#include "regge-generated.h"

constexpr double gate = 1.0/16.0;
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53);

double up(double x) { return std::nextafter(x, INFINITY); }
double down(double x) { return std::nextafter(x, -INFINITY); }

struct Interval {
    double lo, hi;
    Interval(double point = 0): lo(point), hi(point) {}
    Interval(double a, double b): lo(a), hi(b) { assert(lo <= hi); }
    double magnitude() const { return std::max(std::abs(lo),std::abs(hi)); }
};

Interval operator+(Interval a, Interval b) { return {down(a.lo+b.lo),up(a.hi+b.hi)}; }
Interval operator-(Interval a, Interval b) { return {down(a.lo-b.hi),up(a.hi-b.lo)}; }
Interval operator*(Interval a, Interval b) {
    const std::array<double,4> values = {a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    auto extrema = std::minmax_element(values.begin(),values.end());
    return {down(*extrema.first),up(*extrema.second)};
}
Interval operator/(Interval a, Interval b) {
    assert(b.lo > 0);
    return a * Interval(down(1/b.hi),up(1/b.lo));
}
Interval root(Interval a) {
    assert(a.lo >= 0);
    return {std::max(0.,down(std::sqrt(a.lo))),up(std::sqrt(a.hi))};
}

// Monotone table of atan(y/(1-y)), y = x/(1+x). Each tabulated
// endpoint lies between atan_lower[j] and atan_lower[j] + angle_unit.
Interval atan_point(double x) {
    if (x < 0) {
        Interval opposite=atan_point(-x);
        return {-opposite.hi,-opposite.lo};
    }
    Interval y = Interval(x) / (Interval(1)+x);
    int i = std::max(0,static_cast<int>(std::floor(bins*y.lo)));
    int j = std::min(bins,static_cast<int>(std::ceil(bins*y.hi)));
    return {atan_lower[i],atan_lower[j]+angle_unit};
}
Interval angle(Interval cotangent) {
    return {atan_point(cotangent.lo).lo, atan_point(cotangent.hi).hi};
}

struct Face { std::array<int,3> edges; Interval shift; };
std::array<Face,10> faces;

int edge(int i,int j) { return i*(9-i)/2 + j-i-1; }

void initialize_faces() {
    const Interval pi_over_six = Interval(atan_lower[bins],atan_lower[bins]+angle_unit)/3;
    for(int i=0;i<5;i++) for(int j=i+1;j<5;j++) {
        std::array<int,3> v{};
        int count=0;
        for(int k=0;k<5;k++) if(k!=i && k!=j) v[count++]=k;
        std::array<int,3> gaps = {v[1]-v[0],v[2]-v[1],5+v[0]-v[2]};
        bool valence_six = *std::max_element(gaps.begin(),gaps.end())==3;
        faces[edge(i,j)]={{edge(v[0],v[1]),edge(v[0],v[2]),edge(v[1],v[2])},
                          valence_six ? pi_over_six : Interval(0)};
    }
}

struct Level {
    Vector radii{}, gradient_tail{};
    double determinant_tail=0;
    int axis=0;
};
std::array<Level,53> levels;

void initialize_levels(int fixed_edge) {
    Vector radii;
    radii.fill(7./16);
    radii[fixed_edge]=0;
    for(auto& level:levels) {
        level={};
        level.radii=radii;
        level.axis=std::max_element(radii.begin(),radii.end())-radii.begin();
        for(int i=0;i<10;i++) {
            Interval sum;
            for(int j=0;j<10;j++) for(int k=0;k<10;k++)
                sum=sum+Interval(third_bound[i][j][k])*radii[j]*radii[k]/2;
            level.gradient_tail[i]=sum.hi;
            level.determinant_tail=up(level.determinant_tail+up(up(sum.hi*radii[i])/3));
        }
        radii[level.axis]/=2;
    }
}

enum class Outcome { split, empty, bounded };

Outcome classify(const Vector& q,int depth) {
    const Level& level=levels[depth];
    const Vector& r=level.radii;
    std::array<double,66> jet;
    jets(q,jet);
    Vector gradient_error=level.gradient_tail;
    double det_error=level.determinant_tail;
    int entry=11;
    for(int i=0;i<10;i++) {
        det_error=up(det_error+up(std::abs(jet[i+1])*r[i]));
        for(int j=i;j<10;j++) {
            double h=std::abs(jet[entry++]);
            gradient_error[i]=up(gradient_error[i]+up(h*r[j]));
            if(i!=j) gradient_error[j]=up(gradient_error[j]+up(h*r[i]));
            double term=up(up(h*r[i])*r[j]);
            det_error=up(det_error+(i==j ? up(term/2) : term));
        }
    }
    double max_det=up(jet[0]+det_error);
    if(max_det<gate) return Outcome::empty;
    double feasible_min=std::max(gate,down(jet[0]-det_error));

    // Sylvester test at the center (these degree-three dyadics have
    // denominator at most (2*1024)^3 and absolute intermediates below 16).
    double u=q[0],v=q[1],w=q[2];
    double a=(u+v-q[4])/2,b=(u+w-q[5])/2,c=(v+w-q[7])/2;
    bool valid_center=jet[0]>1./128 && u*v>a*a &&
        u*v*w+2*a*b*c-u*c*c-v*b*b-w*a*a>0;
    double segment_min=valid_center ? std::min(feasible_min,jet[0]) : feasible_min;
    Interval segment_det(segment_min,max_det);
    Interval center_score;
    std::array<Interval,10> gradient{};

    for(int f=0;f<10;f++) {
        const auto& indices=faces[f].edges;
        std::array<double,3> x,g;
        for(int k=0;k<3;k++) x[k]=q[indices[k]];
        double heron=2*(x[0]*x[1]+x[1]*x[2]+x[2]*x[0])
                    -(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        double linear=0, radius_sum=0;
        for(int k=0;k<3;k++) {
            g[k]=x[(k+1)%3]+x[(k+2)%3]-x[k];
            radius_sum+=r[indices[k]];
            linear+=2*std::abs(g[k])*r[indices[k]];
        }
        // Exact dyadic degree-two operations here, including the remainder.
        double error=linear+radius_sum*radius_sum;
        double upper=std::min(3.,heron+error); // 16 times area squared
        double lower=std::max(4*segment_min,heron-error);
        if(upper<lower) return Outcome::empty;
        if(!valid_center) continue;
        Interval area=root(Interval(lower,upper))/4;
        Interval cotangent=(Interval(jet[1+f])+Interval(-gradient_error[f],gradient_error[f]))
                           /(Interval(2)*area*root(segment_det));
        Interval deficit=angle(cotangent)-faces[f].shift;
        Interval center_area=root(Interval(heron))/4;
        Interval center_cot=Interval(jet[1+f])/(Interval(2)*center_area*root(jet[0]));
        center_score=center_score+center_area*(angle(center_cot)-faces[f].shift);
        for(int k=0;k<3;k++) {
            Interval partial=Interval(g[k]-radius_sum,g[k]+radius_sum)/(Interval(16)*area);
            int i=indices[k];
            gradient[i]=gradient[i]+deficit*partial;
        }
    }
    if(valid_center) {
        double upper=center_score.hi;
        for(int i=0;i<10;i++) upper=up(upper+up(gradient[i].magnitude()*r[i]));
        Interval threshold=(Interval(391)/200)*root(root(feasible_min)/24);
        if(upper<threshold.lo) return Outcome::bounded;
    }
    return Outcome::split;
}

struct Counts { long long split=0,empty=0,bounded=0; int deepest=0; };

void traverse(Vector& center,int depth,Counts& counts) {
    assert(depth<static_cast<int>(levels.size()));
    counts.deepest=std::max(counts.deepest,depth);
    Outcome result=classify(center,depth);
    if(result==Outcome::empty) ++counts.empty;
    else if(result==Outcome::bounded) ++counts.bounded;
    else {
        ++counts.split;
        int i=levels[depth].axis;
        double middle=center[i],step=levels[depth].radii[i]/2;
        center[i]=middle-step; traverse(center,depth+1,counts);
        center[i]=middle+step; traverse(center,depth+1,counts);
        center[i]=middle;
    }
}

int main() {
    initialize_faces();
    for(int fixed_edge : {0,1}) {
        initialize_levels(fixed_edge);
        constexpr int initial_depth=14;
        constexpr int pieces=1<<initial_depth;
        long long split=0,empty=0,bounded=0;
        int deepest=0;
        #pragma omp parallel for schedule(dynamic,4) reduction(+:split,empty,bounded) reduction(max:deepest)
        for(int id=0;id<pieces;id++) {
            Vector center; center.fill(9./16); center[fixed_edge]=1;
            for(int depth=0;depth<initial_depth;depth++) {
                int i=levels[depth].axis;
                center[i]+=((id>>depth)&1 ? .5 : -.5)*levels[depth].radii[i];
            }
            Counts count;
            traverse(center,initial_depth,count);
            split+=count.split; empty+=count.empty; bounded+=count.bounded;
            deepest=std::max(deepest,count.deepest);
        }
        assert(empty+bounded==split+pieces);
        std::printf("Chart %d: %lld empty, %lld bounded, maximum depth %d\n",
                    fixed_edge,empty,bounded,deepest);
        std::fflush(stdout);
    }
}
