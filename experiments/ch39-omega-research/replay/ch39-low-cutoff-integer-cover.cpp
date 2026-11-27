// 已向外取整的一维表。逐支撑深度优先，只保留当前树枝。
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
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
    std::vector<std::array<int,2>> domain, bins;
    Integer supports=0, examined=0, kinetic=0;

    Cover(int size,int bin_count,int axes) : k(size),leaves(bin_count),projections(axes),
        intervals(layers*2*bin_count*axes), initial(axes),base(axes),drifts(layers*axes),
        active(size),nodes(size,1) {
        bins.resize(2*leaves);
        for(int i=0;i<leaves;++i) bins[i+leaves]={i,i};
        for(int i=leaves-1;i>0;--i) bins[i]={bins[2*i][0],bins[2*i+1][1]};
        read_array(initial.data(),initial.size());
        read_array(drifts.data(),drifts.size());
        read_array(intervals.data(),intervals.size());
    }
    const Range *at(int i,int node) const {
        return &intervals[(i*2*leaves+node)*projections];
    }
    bool check() {
        ++examined;
        for(int i=0;i<static_cast<int>(domain.size());++i){
            const auto range=bins[nodes[i]];
            if(range[0]>domain[i][1] || range[1]<domain[i][0]) return true;
        }
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
    void set_base(){
        base=initial;
        for(int i=0;i<layers;++i){
            if(std::find(active.begin(),active.end(),i)!=active.end()) continue;
            for(int j=0;j<projections;++j){
                base[j][0]-=drifts[i*projections+j];
                base[j][1]+=drifts[i*projections+j];
            }
        }
    }
    bool selected(){
        Integer cases;
        read_array(&cases,1);
        if(cases<1 || cases>100000) throw std::runtime_error("支撑数越界");
        const int max_k=k;
        for(Integer n=0;n<cases;++n){
            Integer size;
            read_array(&size,1);
            if(size<1 || size>max_k) throw std::runtime_error("支撑维数越界");
            k=static_cast<int>(size);
            active.resize(k); nodes.assign(k,1); domain.resize(k);
            std::array<bool,layers> used{};
            for(int i=0;i<k;++i){
                std::array<Integer,3> row{};
                read_array(row.data(),3);
                if(row[0]<1 || row[0]>layers || used[row[0]-1] || row[1]<0 ||
                   row[1]>row[2] || row[2]>=leaves) throw std::runtime_error("支撑或坐标区间无效");
                active[i]=static_cast<int>(row[0]-1);
                used[active[i]]=true;
                domain[i]={static_cast<int>(row[1]),static_cast<int>(row[2])};
            }
            ++supports;
            set_base();
            if(!check()) return false;
        }
        std::cout<<"所选支撑、全部闭盒、途中排除数 "<<supports<<" "<<examined<<" "<<kinetic<<std::endl;
        return true;
    }
    bool all(){
        std::iota(active.begin(),active.end(),0);
        while(true){
            ++supports;
            set_base();
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

int main(int argc,char** argv){
    std::array<Integer,3> dims{};
    read_array(dims.data(),dims.size());
    if(dims[0]<1 || dims[0]>layers || dims[1]<2 || dims[1]>65536 ||
       dims[2]<7 || dims[2]>128 || (dims[1]&(dims[1]-1))) return 2;
    Cover work(dims[0],dims[1],dims[2]);
    if(argc==1) return work.all()?0:1;
    if(argc==2 && std::string(argv[1])=="--selected") return work.selected()?0:1;
    return 2;
}
