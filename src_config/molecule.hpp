#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include <vector>
#include <utility>
#include "atom.hpp"
#include "lattice.hpp"
#include "random.hpp"

class Molecule{
public:
    Molecule(int idx_in){idx=idx_in;}
    ~Molecule(){};

    void setid(int idx_in){this->idx = idx_in;}
    std::vector<Atom*> linkedToCluster() const;
    int countBondCur();
    int countBondNext() const;
    bool isSurface() const;
    int getRandRep() const {return diceI(repNum-1);}
    void move();
    void putBack();
    void linkedAndSurf(std::vector<int>& linked, std::vector<int>& surf);
    virtual bool tryFill(Atom* dest, int repidx)=0;
    
    int idx{0};
    int repNum;
    std::vector<int> freeBondNum;
    std::vector<double> weight;
    std::vector<Atom*> cur, next;
    int bondNum{0};
};

class Adamantane:public Molecule{
public:
    Adamantane(int idx_in);
    ~Adamantane(){};
    bool tryFill(Atom* dest, int repidx);
};

int Molecule::countBondCur(){
    int count = 0;
    for(auto& atom:cur){
        for(auto& nn:atom->NN){
            if(nn->idx!=this->idx and nn->idx!=-1) count++;
        }
    }
    bondNum = count;
    return count;
}

int Molecule::countBondNext() const {
    int count = 0;
    for(auto& atom:next){
        for(auto& nn:atom->NN){
            if(nn->idx!=this->idx and nn->idx!=-1) count++;
        }
    }
    return count;
}

bool Molecule::isSurface() const {
    for(auto& atom:cur)for(auto& nn:atom->NN)if(nn->idx==-1)return true;
    return false;
}

void Molecule::move(){
    // std::cout<<"cur:"<<cur.size()<<",next:"<<next.size()<<"\n";
    cur = std::move(next);next.clear();
    // std::cout<<"cur:"<<cur.size()<<",next:"<<next.size()<<"\n";
    assert(cur.size()==10);
    for(auto& atom:cur) atom->idx = idx;
}

void Molecule::putBack(){
    for(auto& atom:cur) atom->idx = this->idx;
}

void Molecule::linkedAndSurf(std::vector<int>& linked, std::vector<int>& surf){
    linked.clear(); surf.clear();
    for(const auto& atom:cur){
        bool is_linked = false;
        for(const auto& nn:atom->NN){
            if(nn->idx==-1){if(nn->countFilledNN()==1) surf.push_back(nn->pos);}
            else if (!is_linked and nn->idx!=this->idx){linked.push_back(atom->pos); is_linked=true;}
        }
    }
}

std::vector<Atom*> Molecule::linkedToCluster() const {
    std::vector<Atom*> result;
    for(auto atom:cur){
        for(auto nn:atom->NN){
            if(nn->idx!=this->idx and nn->idx!=-1){
                result.push_back(atom);
                break;
            }
        }
    }
    return result;
}

Adamantane::Adamantane(int idx_in):Molecule(idx_in){
    repNum = 2;
    freeBondNum = std::vector<int> {1,2};
    weight = std::vector<double> {4.0/16.0, 12.0/16.0};
}

bool Adamantane::tryFill(Atom* dest, int repidx){
    std::vector<int> filledNN, emptyNN;
    dest->classNN(filledNN, emptyNN);
    if(repidx==0){
        if(filledNN.size()>1) return false;
        next.clear();
        next.push_back(dest);
        int i0;
        if(!filledNN.empty()) i0 = filledNN[0];
        else {int i = diceI(3); i0 = emptyNN[i]; emptyNN.erase(emptyNN.begin()+i);}
        for(int i = 0; i < 3; i++){
            int i1 = emptyNN[i];
            int i2 = emptyNN[(i+1)%3];
            Atom* nn1 = dest->NN[i1];
            Atom* nn2 = nn1->NN[i0];
            Atom* nn3 = nn2->NN[i2];
            if(nn2->idx==-1 and nn3->idx==-1){next.push_back(nn1);next.push_back(nn2);next.push_back(nn3);}
            else return false;
        }
        assert(next.size()==10);
        return true;
    }else{
        int fillsize = filledNN.size();
        if(fillsize>2) return false;
        next.clear();
        next.push_back(dest);
        for(int i = 0; i<(2-fillsize); i++){
            int j = diceI(emptyNN.size()-1);
            filledNN.push_back(emptyNN[j]);
            emptyNN.erase(emptyNN.begin()+j);
        }

        int i0 = filledNN[0], i1 = filledNN[1], i2 = emptyNN[0], i3 = emptyNN[1];
        auto nn = dest->NN[i2]; next.push_back(nn);
        auto nn1 = nn->NN[i0]; if(nn1->idx==-1) next.push_back(nn1); else return false;
        auto nn2 = nn1->NN[i3]; if(nn2->idx==-1) next.push_back(nn2); else return false;
        auto nn3 = nn2->NN[i1]; if(nn3->idx==-1) next.push_back(nn3); else return false; 
        nn1 = nn->NN[i1]; if(nn1->idx==-1) next.push_back(nn1); else return false;
        nn2 = nn1->NN[i3]; if(nn1->idx==-1) next.push_back(nn2); else return false;
        
        nn = dest->NN[i3]; next.push_back(nn);
        nn1 = nn->NN[i0]; if(nn1->idx==-1) next.push_back(nn1); else return false;
        nn1 = nn->NN[i1]; if(nn1->idx==-1) next.push_back(nn1); else return false;
        assert(next.size()==10);
        return true;
    }
}


#endif // __MOLECULE_H__