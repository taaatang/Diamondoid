#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <vector>
#include <algorithm>
#include <queue>
#include <assert.h>
#include "atom.hpp"
#include "algebra.hpp"

struct Lattice{
    Lattice(int na, int nb, int nc);
    ~Lattice(){}

    void clustering( );
    void linkNN(Atom& atom, const std::array<std::array<int,3>,4>& NN);
    void linkNN();

    int Na{0}, Nb{0}, Nc{0};
    Atom* center;
    std::vector<Atom> latt;
    int empty{0}, occupied{0};
};

Lattice::Lattice(int na, int nb, int nc):Na(na),Nb(nb),Nc(nc){
    latt.clear();
    for(int ia=-Na; ia<=Na; ia++)for(int ib=-Nb; ib<=Nb; ib++)for(int ic=-Nc; ic<=Nc; ic++){
        Atom tmpa(Type::A, A0+ia*a+ib*b+ic*c);
        Atom tmpb(Type::B, B0+ia*a+ib*b+ic*c);
        latt.push_back(tmpa);
        latt.push_back(tmpb);
    }
    std::sort(latt.begin(), latt.end());
    for(int i = 0; i < latt.size(); i++) latt[i].pos = i;
    Atom tmp(Type::A, A0);
    auto low = std::lower_bound(latt.begin(), latt.end(), tmp);
    assert(low!=latt.end());
    center = &(*low);
    linkNN();
}

void Lattice::clustering( ) {
    unvisit(&latt);
    empty = 0;
    occupied = 0;
    for (auto& atom : latt) {
        if (atom.is_visited) continue;
        int count = 0;
        std::queue<Atom*> q;
        q.push(&atom);
        atom.is_visited = true;
        bool is_empty = (atom.idx==-1);
        while (!q.empty()) {
            auto a = q.front();
            q.pop();
            ++count;
            for (auto& an:a->NN) {
                if (!an or an->is_visited) continue;
                if (is_empty == (an->idx==-1)) {
                    an->is_visited = true;
                    q.push(an);
                }
            }
        }
        if (is_empty) {
            ++empty;
            std::cout<<"empty cluster "<<empty<<" size: "<<count<<"\n";
        } else {
            ++occupied;
            std::cout<<"occupied cluster "<<occupied<<" size: "<<count<<"\n";
        }
    }
}

void Lattice::linkNN(Atom& atom, const std::array<std::array<int,3>,4>& NN){
    Atom tmp = atom;
    for(int idx = 0; idx < 4; idx++){
        tmp.coord = atom.coord + NN[idx];
        auto low = std::lower_bound(latt.begin(), latt.end(),tmp);
        if(low==latt.end() or (*low).coord>tmp.coord) atom.NN[idx] = nullptr;
        else atom.NN[idx]= &(*low);
    }
}

void Lattice::linkNN(){
    for(auto& atom:latt){
        if(atom.type==Type::A)linkNN(atom,ANN);
        else if(atom.type==Type::B)linkNN(atom,BNN);
    }
}

#endif // __LATTICE_H__