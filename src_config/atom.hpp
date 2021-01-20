#ifndef __ATOM_H__
#define __ATOM_H__

#include <array>
const std::array<std::array<int,3>,4> A = {{{{1,1,1}},{{-1,-1,1}},{{1,-1,-1}},{{-1,1,-1}}}};
const std::array<std::array<int,3>,4> B = {{{{-1,-1,-1}},{{1,1,-1}},{{-1,1,1}},{{1,-1,1}}}};
struct Atom{
    Atom(bool state, std::array<int,3> coord_in):is_full(state){coord = coord; for(auto& el:NN)el=nullptr;}
    ~Atom(){}
    void addNN(int idx, Atom* NNptr){NN[idx]=NNptr;}
    bool is_full{false};
    std::array<int,3> coord;
    std::array<Atom*,4> NN;
};
inline bool operator<(const Atom& lhs, const Atom& rhs){
    return lhs.coord<rhs.coord;
}
inline bool operator>(const Atom& lhs, const Atom& rhs){
    return lhs.coord>rhs.coord;
}
inline bool operator==(const Atom& lhs, const Atom& rhs){
    return lhs.coord==rhs.coord;
}
#endif // __ATOM_H__