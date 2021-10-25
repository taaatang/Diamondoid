#pragma once

#include <array>
#include <vector>
using VecI=std::vector<int>;

/// Two nonequivalent carbon atoms in diamond crystal
enum class Type{A,B};
/// Type A atom neighbor relative coordinates
const std::array<std::array<int,3>,4> ANN = {{{{1,1,1}},{{-1,-1,1}},{{1,-1,-1}},{{-1,1,-1}}}};
/// Type B atom neighbor relative coordinates
const std::array<std::array<int,3>,4> BNN = {{{{-1,-1,-1}},{{1,1,-1}},{{-1,1,1}},{{1,-1,1}}}};
/// Type A atom origin coordinates
const std::array<int,3> A0 = {{0,0,0}};
/// Type B atom origin coordinates
const std::array<int,3> B0 = {{1,1,1}};
/// Basis vector for diamond crystal
const std::array<int,3> a = {{2,2,0}}, b = {{0,2,2}}, c = {{2,0,2}};
struct Atom{
    Atom() = default;
    Atom(Type type_in, std::array<int,3> coord_in, int idx_in=-1):idx(idx_in), type(type_in){coord = coord_in; for(auto& el:NN)el=nullptr;}
    ~Atom() = default;
    void classNN(std::vector<int>& filled, std::vector<int>& empty);
    int countFilledNN() const {int count = 0; for(const auto& nn:NN){if(nn and nn->idx!=-1)count++;} return count;};
    int pos{0}; ///< postion in Lattice::latt vector
    int idx{-1}; ///< belongs to idx-th molecule. idx=-1 means empty
    bool is_visited{false};
    int clusterIdx{0};///< label same cluster
    Type type{Type::A};///< two carbon atoms in a unit cell
    std::array<int,3> coord{0,0,0}; ///< coord in a,b,c. a=[4,0,0], b=[0,4,0], c=[0,0,4]
    std::array<Atom*,4> NN{nullptr, nullptr, nullptr, nullptr}; ///< pointer to nearest neighbors
};

void unvisit(std::vector<Atom>* mol,  const VecI &visited_idx){
    for(auto idx:visited_idx){
        mol->at(idx).is_visited=false;
    }
}

void unvisit(std::vector<Atom>* mol){
    for(auto& atom:(*mol)){
        atom.is_visited=false;
    }
}

/// Classify neighboring sites to filled and empty states
void Atom::classNN(std::vector<int>& filled, std::vector<int>& empty){
    filled.clear();
    empty.clear();
    for (int i = 0; i < 4; i++) {
        if (NN[i]) {
            if (NN[i]->idx == -1) empty.push_back(i);
            else filled.push_back(i);
        }
    }
}

inline bool operator<(const Atom& lhs, const Atom& rhs){
    return lhs.coord<rhs.coord;
}
inline bool operator>(const Atom& lhs, const Atom& rhs){
    return lhs.coord>rhs.coord;
}
inline bool operator==(const Atom& lhs, const Atom& rhs){
    return lhs.coord==rhs.coord;
}