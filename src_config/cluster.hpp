#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "../src/utils.hpp"

#include "atom.hpp"
#include "molecule.hpp"
#include "lattice.hpp"
#include "random.hpp"

const int MAX_TRY = 10;
const double TEMP = 0.2;

class Cluster{
public:
    Cluster(int na, int nb, int nc);
    ~Cluster(){};

    Atom* getRandPos();

    bool add(Molecule* mol);
    void push_back(Molecule* mol){structure.push_back(mol);}
    void computeSurf();
    void addPos(Molecule* mol);
    void rmvPos(Molecule* mol);
    void recenter();
    int countBond();

    void init(std::vector<Adamantane>& mols);
    void singleStep();
    void compute();

    void saveCoords(std::string filename);

private:
    Lattice Latt;
    int curMolIdx;
    int totBond{0};
    std::vector<Molecule*> structure;    
    std::vector<int> surface;
    std::vector<int> positions;
};

Cluster::Cluster(int na, int nb, int nc):Latt(na,nb,nc),curMolIdx(-1){
    positions.push_back(Latt.center->pos);
}

void Cluster::init(std::vector<Adamantane>& mols){
    std::cout<<"\nBegin init...\n";
    int count = 0;
    for(auto& mol:mols) {
        mol.idx = count; count++;
        if(!add(&mol)){std::cout<<"Failed to add a molecule!\n"; exit(1);}
        else structure.push_back(&mol);
    }
    computeSurf();
    std::cout<<"Finished init. Total Bond: "<<countBond()<<".\n";
}

void Cluster::singleStep(){
    auto i = diceI(surface.size()-1);
    auto molIdx = surface[i];
    auto mol = structure.at(molIdx);
    rmvPos(mol);
    if(add(mol)){
        computeSurf();
    }else{
        mol->putBack();
        addPos(mol);
    }
}

Atom* Cluster::getRandPos(){
    return &Latt.latt.at(positions.at(diceI(positions.size()-1)));
}

bool Cluster::add(Molecule* mol){
    // mol->setid(structure.size());
    int count = 0;
    bool flag = false;
    while(count<MAX_TRY){
        count++;
        // std::cout<<"Try add mol "<<mol->idx<<", count:"<<count<<"\n";
        auto dest = getRandPos();
        if(mol->tryFill(dest,mol->getRandRep())){
            auto bondNum = mol->countBondNext();
            if(bondNum>mol->bondNum){flag=true;break;}
            else{
                if(diceD()<std::exp((bondNum-mol->bondNum)/TEMP)){flag=true;break;}
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
    for(int i=0; i<structure.size(); i++){
        if(structure[i]->isSurface()) surface.push_back(i);
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
    for(int i = 1; i<structure[0]->cur.size(); i++)save<int>(structure[0]->cur[i]->coord.data(),3,&outfile,filename,true);
    for(int m = 1; m<structure.size(); m++){
        for(int i = 0; i<structure[m]->cur.size();i++)save<int>(structure[m]->cur[i]->coord.data(),3,&outfile,filename,true);
    }
}
#endif // __CLUSTER_H__