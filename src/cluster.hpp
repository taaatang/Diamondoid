#pragma once

#include <vector>

#include "utils.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "lattice.hpp"
#include "random.hpp"

class Cluster{
public:
    Cluster(int na, int nb, int nc);

    Atom* getRandPos();

    bool add(int step, Molecule* mol, int tryNum = 1, double d = -1.0);
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
    void singleStep(int step, double d = -1.0);

    void clustering(int size = 100);

    void saveCoords(const std::string& filename);
    void saveSurface(const std::string& filename);
    void saveVacancy(const std::string& filename);

    bool detect_config_saved(int step) {
        return step > 10000;
    }

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
    bool m_is_config_saved{false};
};

void boundingBox(const std::vector<Atom*>& molecule, arr<3>& minimum, arr<3>& maximum, int empty = -100);