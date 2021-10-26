#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "utils.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "lattice.hpp"
#include "random.hpp"

const int INIT_MAX_TRY = 5000;
const int MAX_TRY = 2;

class Cluster{
public:
    Cluster(int na, int nb, int nc);

    Atom* getRandPos();

    bool add(Molecule* mol, int tryNum = 1, double d = -1.0);
    void push_back(Molecule* mol) { structure.push_back(mol); }
    void computeSurf();
    void computePos();
    void evalPos();
    void addPos(Molecule* mol);
    void rmvPos(Molecule* mol);
    int countBond();
    void getBox(Molecule* excluded, arr<3> &minimum, arr<3> &maximum) const;
    double getVol() const { return volume; }

    void setTemp(double T) { temp = T; }
    void setPressure(double P) { pressure = P; }
    void init(std::vector<Molecule*>& mols);
    void singleStep(double d = -1.0);

    void clustering(int size = 100);

    void saveCoords(const std::string& filename);
    void saveSurface(const std::string& filename);
    void saveVacancy(const std::string& filename);

private:
    double temp{1.0};
    double pressure{0.0};
    double volume;
    Lattice Latt;
    int curMolIdx;
    int totBond{0};
    std::vector<Molecule*> structure;    
    std::vector<int> surface;
    std::vector<int> positions;
    std::vector<int> compatible, uncompatible;
    std::vector<int> vacancy;
};

void boundingBox(const std::vector<Atom*>& molecule, arr<3>& minimum, arr<3>& maximum, int empty = -100) {
    if (molecule.empty()) {
        return;
    }
    if (minimum[0] == empty) {
        minimum = molecule[0]->coord;
    }
    if (maximum[0] == empty) {
        maximum = molecule[0]->coord;
    }
    for (const auto& atom : molecule) {
        min(minimum, atom->coord);
        max(maximum, atom->coord);
    }
}

/// Find vacancy
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
    setPressure(0.0);
    for(auto &mol : mols) {
        mol->idx = count; 
        count++;
        if(!add(mol, INIT_MAX_TRY)) {
            std::cout << "Failed to add a molecule!\n";
            exit(1);
        }
        else structure.push_back(mol);
    }
    computeSurf();
    computePos();
    arr<3> bmin;
    arr<3> bmax;
    getBox(nullptr, bmin, bmax);
    auto box = bmax - bmin;
    volume = getVolume(box);
    std::cout << "Finished init. Total Bond: " << countBond() << ". Cluster box " << box << std::endl;
}

void Cluster::singleStep(double d) {
    // random pick from surface molecules
//    auto i = diceI(int(surface.size())-1);
//    auto molIdx = surface[i];
    auto molIdx = diceI(int(structure.size()) - 1); ///< random pick from all atoms
    auto mol = structure.at(molIdx);
    rmvPos(mol);
    if(add(mol, MAX_TRY, d)) {
        computeSurf();
    } else {
        mol->putBack();
        addPos(mol);
    }
}

/// Get a random empty atom position on the surface of the cluster
Atom* Cluster::getRandPos() {
    return &Latt.latt.at(positions.at(diceI(positions.size()-1)));
}

/// Try add molecule on to the surface.
bool Cluster::add(Molecule* mol, int tryNum, double d) {
    // mol->setid(structure.size());
    int count = 0;
    bool flag = false;
    while(count < tryNum) {
        count++;
        // std::cout<<"Try add mol "<<mol->idx<<", count:"<<count<<"\n";
        auto dest = getRandPos();
        // if(mol->tryFill(dest,mol->getRandRep())){
        if(mol->tryFill(dest)) {
            // limit the distance a molecule can move in a single step
            if (d > 0 && mol->jumpDistance() > d) {
                continue;
            }
            arr<3> minB, maxB;
            arr<3> minM, maxM;
            minM.fill(-10000);
            maxM.fill(-10000);
            getBox(mol, minB, maxB);
//            std::cout << "minB:" << minB << ", maxB:" << maxB;
            boundingBox(mol->next, minM, maxM, -10000);
//            std::cout << "\nnext molecule positions:\n";
//            for (const auto& atom : mol->next) {
//                std::cout << "atom:" << atom->coord << std::endl;
//            }
//            std::cout << ", minM:" << minM << ", maxM:" << maxM << std::endl << std::endl;
            min(minB, minM);
            max(maxB, maxM);
            auto box = maxB - minB;
            auto volumeNext = getVolume(box);
            auto bondNumNext = mol->countBondNext();
            double Ei = pressure * volume - mol->bondNum;
            double Ef = pressure * volumeNext - bondNumNext;
//            std::cout << "Ei: " << Ei << "Vi, :" << volume << ", Bi:" << mol->bondNum << "-> Ef: " << Ef << ", Vf:" << volumeNext << ", Bf:" << bondNumNext << std::endl;
            if(diceD() < std::exp((Ei - Ef) / temp)) {
                flag = true;
                volume = volumeNext;
                break;
            }
        }
    }
    if(flag) {
        mol->move();
        addPos(mol); 
        // std::cout<<"Succeed Jump!\n";
    }
    // else std::cout<<"Failed Jump!\n";
    return flag;
}

/// Find all molecules on the surface
void Cluster::computeSurf(){
    surface.clear();
    for(int i=0; i < (int)structure.size(); i++){
        if(structure[i]->isSurface()) surface.push_back(i);
    }
}

/// Find all surface positions (available to attach new molecule)
void Cluster::computePos(){
    positions.clear();
    for (auto& mol : structure) {
        auto surf = mol->surf();
        positions.insert(positions.end(), surf.begin(), surf.end());
    }
    std::sort(positions.begin(), positions.end());
    std::unique(positions.begin(), positions.end());
}

/// Find compatible/uncompatible surface positions
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

/// Update empty surface positions after adding a molecule
void Cluster::addPos(Molecule* mol){
    std::vector<int> linked, surf;
    mol->linkedAndSurf(linked,surf);
    for(auto val:linked){
        auto low = std::lower_bound(positions.begin(),positions.end(),val);
        if(low!=positions.end() and *low == val) positions.erase(low);
    }
    positions.insert(positions.end(),surf.begin(),surf.end());
    std::sort(positions.begin(),positions.end());
    // auto last = std::unique(positions.begin(),positions.end());
    // positions.erase(last,positions.end());
}

/// Remove empty surface positions after remove a molecule
void Cluster::rmvPos(Molecule* mol){
    std::vector<int> linked, surf;
    mol->linkedAndSurf(linked,surf); 
    for(auto val:surf){
        auto low = std::lower_bound(positions.begin(),positions.end(),val);
        if(low!=positions.end() and *low == val) positions.erase(low);
    }
    positions.insert(positions.end(),linked.begin(),linked.end());
    std::sort(positions.begin(),positions.end());
    for(auto& atom:mol->cur) atom->idx = -1;
}

/// Count total C-C bonds
int Cluster::countBond(){
    totBond = 0;
    for (const auto &mol : structure) totBond += mol->countBondCur();
    totBond /= 2;
    return totBond;
}

/// Save all atom coordinates
void Cluster::saveCoords(const std::string& filename){
    if(structure.empty()) return;
    std::ofstream outfile;
    save<int>(structure[0]->cur[0]->coord.data(),3,&outfile,filename);
    for(int i = 1; i < (int)structure[0]->cur.size(); i++) {
        save<int>(structure[0]->cur[i]->coord.data(),3,&outfile,filename,true);
    }
    for(int m = 1; m < (int)structure.size(); m++) {
        for(auto & i : structure[m]->cur) {
            save<int>(i->coord.data(),3,&outfile,filename,true);
        }
    }
}

/// Save compatible/incompatible surface positions
void Cluster::saveSurface(const std::string& filename) {
    std::ofstream outfile;
    save<int>(nullptr, 0, &outfile, filename + "_comp.dat"); 
    for (int i : compatible) {
        save<int>(Latt.latt.at(i).coord.data(), 3, &outfile, filename+"_comp.dat", true);
    }
    save<int>(nullptr, 0, &outfile, filename + "_uncomp.dat");
    for (int i : uncompatible) {
        save<int>(Latt.latt.at(i).coord.data(), 3, &outfile, filename+"_uncomp.dat", true);
    }
}

/// Save vacancy coordinates
void Cluster::saveVacancy(const std::string& filename) {
    std::ofstream outfile;
    save<int>(nullptr, 0, &outfile, filename);
    for (int i : vacancy) {
        save<int>(Latt.latt.at(i).coord.data(), 3, &outfile, filename, true);
    }
}

/// Find bounding box for all atoms
void Cluster::getBox(Molecule *excluded, arr<3> &minimum, arr<3> &maximum) const {
    int empty = -1000000;
    minimum.fill(empty);
    maximum.fill(empty);
    for (const auto &mol: structure) {
        if (mol != excluded) {
            boundingBox(mol->cur, minimum, maximum, empty);
        }
    }
}