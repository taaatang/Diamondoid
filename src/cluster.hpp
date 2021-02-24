#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "../utils/utils.hpp"

#include "atom.hpp"
#include "molecule.hpp"
#include "lattice.hpp"
#include "random.hpp"

const int INIT_MAX_TRY = 5000;
const int MAX_TRY = 2;

class Cluster{
public:
    Cluster(int na, int nb, int nc);
    ~Cluster(){};

    Atom* getRandPos();

    bool add(Molecule* mol, int tryNum = 1);
    void push_back(Molecule* mol){structure.push_back(mol);}
    void computeSurf();
    void computePos();
    void evalPos();
    void addPos(Molecule* mol);
    void rmvPos(Molecule* mol);
    void recenter();
    int countBond();

    void setTemp(double T) { temp = T; }
    void init(std::vector<Molecule*>& mols);
    void singleStep();
    void compute();

    void clustering(int size = 100);

    void saveCoords(std::string filename);
    void saveSurface(std::string filename);
    void saveVacancy(std::string filename);

private:
    double temp;
    Lattice Latt;
    int curMolIdx;
    int totBond{0};
    std::vector<Molecule*> structure;    
    std::vector<int> surface;
    std::vector<int> positions;
    std::vector<int> compatible, uncompatible;
    std::vector<int> vacancy;
};

void Cluster::clustering(int size) {
    unvisit(&Latt.latt);
    int empty = 0;
    int occupied = 0;
    vacancy.clear();
    for (auto& atom : Latt.latt) {
        if (atom.is_visited) continue;
        int count = 0;
        VecI pos;
        std::queue<Atom*> q;
        q.push(&atom);
        atom.is_visited = true;
        bool is_empty = (atom.idx==-1);
        while (!q.empty()) {
            auto a = q.front();
            q.pop();
            ++count;
            if (is_empty and (int)pos.size()<size) {
                pos.push_back(a->pos);
            }
            for (auto& an:a->NN) {
                if (!an or an->is_visited) continue;
                if (is_empty == (an->idx==-1)) {
                    an->is_visited = true;
                    q.push(an);
                }
            }
        }
        if (is_empty) {
            if ((int)pos.size()<size) {
                vacancy.insert(vacancy.end(), pos.begin(), pos.end());
            }
            ++empty;
            // std::cout<<"empty cluster "<<empty<<" size: "<<count<<"\n";
        } else {
            ++occupied;
            // std::cout<<"occupied cluster "<<occupied<<" size: "<<count<<"\n";
        }
    }
}

Cluster::Cluster(int na, int nb, int nc):Latt(na,nb,nc),curMolIdx(-1){
    positions.push_back(Latt.center->pos);
}

void Cluster::init(std::vector<Molecule*>& mols){
    std::cout<<"\nBegin init...\n";
    int count = 0;
    for(auto& mol:mols) {
        mol->idx = count; count++;
        if(!add(mol, INIT_MAX_TRY)){
            std::cout<<"Failed to add a molecule!\n";
            exit(1);
        }
        else structure.push_back(mol);
    }
    computeSurf();
    computePos();
    std::cout<<"Finished init. Total Bond: "<<countBond()<<".\n";
}

void Cluster::singleStep(){
    auto i = diceI(surface.size()-1);
    auto molIdx = surface[i];
    auto mol = structure.at(molIdx);
    rmvPos(mol);
    if(add(mol, MAX_TRY)){
        computeSurf();
    }else{
        mol->putBack();
        addPos(mol);
    }
}

Atom* Cluster::getRandPos(){
    return &Latt.latt.at(positions.at(diceI(positions.size()-1)));
}

bool Cluster::add(Molecule* mol, int tryNum){
    // mol->setid(structure.size());
    int count = 0;
    bool flag = false;
    while(count<tryNum){
        count++;
        // std::cout<<"Try add mol "<<mol->idx<<", count:"<<count<<"\n";
        auto dest = getRandPos();
        // if(mol->tryFill(dest,mol->getRandRep())){
        if(mol->tryFill(dest)){
            auto bondNum = mol->countBondNext();
            if(bondNum>mol->bondNum){flag=true;break;}
            else{
                if(diceD()<std::exp((bondNum-mol->bondNum)/temp)){flag=true;break;}
            }
        }
    }
    if(flag) {
        mol->move();addPos(mol); 
        // std::cout<<"Succeed Jump!\n";
    }
    // else std::cout<<"Failed Jump!\n";
    return flag;
}

void Cluster::computeSurf(){
    surface.clear();
    for(int i=0; i < (int)structure.size(); i++){
        if(structure[i]->isSurface()) surface.push_back(i);
    }
}

void Cluster::computePos() {
    positions.clear();
    for (auto& mol : structure) {
        auto surf = mol->surf();
        positions.insert(positions.end(), surf.begin(), surf.end());
    }
    std::sort(positions.begin(), positions.end());
    std::unique(positions.begin(), positions.end());
}

void Cluster::evalPos() {
    compatible.clear();
    uncompatible.clear();
    Molecule* mol = structure.at(0);
    for (auto pos : positions) {
        if (mol->tryAll(&Latt.latt.at(pos))) {
            compatible.push_back(pos);
        } else {
            uncompatible.push_back(pos);
        }
    }
}

void Cluster::addPos(Molecule* mol){
    std::vector<int> linked, surf;
    mol->linkedAndSurf(linked,surf);
    for(auto val:linked){
        auto low = std::lower_bound(positions.begin(),positions.end(),val);
        if(low!=positions.end()) positions.erase(low);
    }
    positions.insert(positions.end(),surf.begin(),surf.end());
    std::sort(positions.begin(),positions.end());
    // auto last = std::unique(positions.begin(),positions.end());
    // positions.erase(last,positions.end());
}

void Cluster::rmvPos(Molecule* mol){
    std::vector<int> linked, surf;
    mol->linkedAndSurf(linked,surf); 
    for(auto val:surf){
        auto low = std::lower_bound(positions.begin(),positions.end(),val);
        if(low!=positions.end()) positions.erase(low);
    }
    positions.insert(positions.end(),linked.begin(),linked.end());
    std::sort(positions.begin(),positions.end());
    for(auto& atom:mol->cur) atom->idx = -1;
}

int Cluster::countBond(){
    totBond = 0;
    for(auto& mol:structure)totBond+=mol->countBondCur();
    totBond /= 2;
    return totBond;
}

void Cluster::saveCoords(std::string filename){
    if(structure.empty()) return;
    std::ofstream outfile;
    save<int>(structure[0]->cur[0]->coord.data(),3,&outfile,filename);
    for(int i = 1; i < (int)structure[0]->cur.size(); i++)save<int>(structure[0]->cur[i]->coord.data(),3,&outfile,filename,true);
    for(int m = 1; m < (int)structure.size(); m++){
        for(int i = 0; i<(int)structure[m]->cur.size();i++)save<int>(structure[m]->cur[i]->coord.data(),3,&outfile,filename,true);
    }
}

void Cluster::saveSurface(std::string filename) {
    std::ofstream outfile;
    save<int>(nullptr, 0, &outfile, filename + "_comp.dat"); 
    for (int i = 0; i < (int)compatible.size(); ++i) {
        save<int>(Latt.latt.at(compatible[i]).coord.data(), 3, &outfile, filename+"_comp.dat", true); 
    }
    save<int>(nullptr, 0, &outfile, filename + "_uncomp.dat");
    for (int i = 0; i < (int)uncompatible.size(); ++i) {
        save<int>(Latt.latt.at(uncompatible[i]).coord.data(), 3, &outfile, filename+"_uncomp.dat", true); 
    }

    // if (!compatible.empty()) {
    //     save<int>(Latt.latt.at(compatible[0]).coord.data(), 3, &outfile, filename + "_comp.dat");
    //     for (int i = 1; i < compatible.size(); ++i) {
    //        save<int>(Latt.latt.at(compatible[i]).coord.data(), 3, &outfile, filename+"_comp.dat", true); 
    //     }
    // }
    // if (!uncompatible.empty()) {
    //     save<int>(Latt.latt.at(uncompatible[0]).coord.data(), 3, &outfile, filename + "_uncomp.dat");
    //     for (int i = 1; i < uncompatible.size(); ++i) {
    //        save<int>(Latt.latt.at(uncompatible[i]).coord.data(), 3, &outfile, filename+"_uncomp.dat", true); 
    //     }
    // }
}

void Cluster::saveVacancy(std::string filename) {
    std::ofstream outfile;
    save<int>(nullptr, 0, &outfile, filename);
    for (int i = 0; i < (int)vacancy.size(); ++i) {
        save<int>(Latt.latt.at(vacancy[i]).coord.data(), 3, &outfile, filename, true); 
    }
}
#endif // __CLUSTER_H__