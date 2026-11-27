// 已向外取整的一维表。逐支撑深度优先，只保留当前树枝。
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Integer = std::int64_t;
constexpr int layers = 19;
using Range = std::array<Integer, 2>;

template<class T> void read_array(T *target, std::size_t n) {
    std::cin.read(reinterpret_cast<char *>(target), sizeof(T)*n);
    if (!std::cin) throw std::runtime_error("区间表不完整");
}

struct Cover {
    int k, leaves, projections;
    std::vector<Range> intervals, initial, base;
    std::vector<Integer> drifts;
    std::vector<int> active, nodes;
    Integer supports=0, examined=0, kinetic=0;

    Cover(int size,int bins,int axes) : k(size),leaves(bins),projections(axes),
        intervals(layers*2*bins*axes), initial(axes),base(axes),drifts(layers*axes),
        active(size),nodes(size,1) {
        read_array(initial.data(),initial.size());
        read_array(drifts.data(),drifts.size());
        read_array(intervals.data(),intervals.size());
    }
    const Range *at(int i,int node) const {
        return &intervals[(i*2*leaves+node)*projections];
    }
    bool check() {
        ++examined;
        std::array<const Range*,layers> entries{};
        for(int i=0;i<k;++i) entries[i]=at(active[i],nodes[i]);
        for(int j=0;j<projections;++j){
            Integer low=base[j][0], high=base[j][1];
            for(int i=0;i<k;++i){
                low+=entries[i][j][0];
                high+=entries[i][j][1];
            }
            if(high<0 || (j<6 && low>0)){
                if(j>=6) ++kinetic;
                return true;
            }
        }
        int split=-1;
        Integer best=-1;
        for(int i=0;i<k;++i){
            if(nodes[i]>=leaves) continue;
            const auto width=[&](int j){return entries[i][j][1]-entries[i][j][0];};
            const auto score=width(0)+20*width(2)+150*width(5);
            if(score>best){best=score;split=i;}
        }
        if(split<0){
            std::cout<<"尚未排除的层和对数格";
            for(int i=0;i<k;++i) std::cout<<" "<<active[i]+1<<"/"<<nodes[i]-leaves;
            std::cout<<std::endl;
            for(int j=0;j<projections;++j){
                Integer low=base[j][0],high=base[j][1];
                for(int i=0;i<k;++i){
                    low+=entries[i][j][0];
                    high+=entries[i][j][1];
                }
                std::cout<<" "<<j<<"/"<<low/4294967296.0<<","<<high/4294967296.0;
            }
            std::cout<<std::endl;
            return false;
        }
        int old=nodes[split];
        nodes[split]=2*old;
        if(!check()) return false;
        ++nodes[split];
        if(!check()) return false;
        nodes[split]=old;
        return true;
    }
    bool all(){
        std::iota(active.begin(),active.end(),0);
        while(true){
            ++supports;
            base=initial;
            for(int i=0;i<layers;++i){
                if(std::binary_search(active.begin(),active.end(),i)) continue;
                for(int j=0;j<projections;++j){
                    base[j][0]-=drifts[i*projections+j];
                    base[j][1]+=drifts[i*projections+j];
                }
            }
            if(!check()) return false;
            int i=k-1;
            while(i>=0 && active[i]==layers-k+i) --i;
            if(i<0) break;
            ++active[i];
            for(int j=i+1;j<k;++j) active[j]=active[j-1]+1;
        }
        std::cout<<"支撑、全部闭盒、途中正性排除数 "<<supports<<" "<<examined<<" "<<kinetic<<std::endl;
        return true;
    }
};

int main(){
    std::array<Integer,3> dims{};
    read_array(dims.data(),dims.size());
    if(dims[0]<1 || dims[0]>layers || dims[1]<2 || dims[1]>65536 ||
       dims[2]<7 || dims[2]>128 || (dims[1]&(dims[1]-1))) return 2;
    Cover work(dims[0],dims[1],dims[2]);
    return work.all()?0:1;
}
