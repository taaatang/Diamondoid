#pragma once

#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include "atom.hpp"
#include "lattice.hpp"
#include "random.hpp"

class Molecule {
public:
    Molecule(int idx_in) {
        idx = idx_in;
        idmap = VecI(4, 0);
    }

    void setid(int idx_in) { this->idx = idx_in; }

    std::vector<Atom *> linkedToCluster() const; ///< find atoms in molecule linked to the bulk
    int countBondCur(); ///< count C-C bond in current configuration
    int countBondNext() const; ///< count C-C bond in next configuration
    bool isSurface() const; ///< judge if molecule has empty nearest C position
    int getRandRep() const;

    void move(); ///< move to next position
    void putBack(); ///< go back to current position
    void linkedAndSurf(VecI &linked, VecI &surf); ///< linked: atom in cur linked to bulk; surf: empty C only linked to cur;
    VecI surf();

    void setNN(int from, int bondid, int to);

    void setNN(int from, VecI vto);

    void setmolid() {
        for (int i = 0; i < (int) mol->size(); i++) {
            mol->at(i).pos = i;
        }
    }

    void genidmap(int root, Atom *dest);

    bool tryFill(Atom *dest);

    bool tryAll(Atom *dest);

    std::array<double, 3> computeCM(std::vector<Atom *> &atoms);

    double jumpDistance();
    // virtual bool tryFill(Atom* dest, int repidx)=0;

    void saveCoords(const std::string& filename, bool is_app=false);

    int idx{0};
    std::vector<Atom> *mol{nullptr};
    VecI idmap; // map mol bond idx to lattice bondidx. idmap[mol_bidx] = latt_bidx
    std::vector<Atom *> cur, next;
    std::array<double, 3> cm; // center of mass
    int bondNum{0};
};

class Adamantane : public Molecule {
public:
    Adamantane(int idx_in, std::vector<Atom> *mol_);
    // bool tryFill(Atom* dest, int repidx);
    inline static bool mol_is_set{false};
};

class Diamantane : public Molecule {
public:
    Diamantane(int idx_in, std::vector<Atom> *mol);
    inline static bool mol_is_set{false};
};

class Triamantane : public Molecule {
public:
    Triamantane(int idx_in, std::vector<Atom> *mol_);
    inline static bool mol_is_set{false};
};

class Tetramantane : public Molecule {
public:
    Tetramantane(int idx_in, std::vector<Atom> *mol_);
    inline static bool mol_is_set{false};
};

class Pentamantane : public Molecule {
public:
    Pentamantane(int idx_in, std::vector<Atom> *mol_);
    inline static bool mol_is_set{false};
};

class Pentamantane1212 : public Molecule {
public:
    Pentamantane1212(int idx_in, std::vector<Atom> *mol_);
    inline static bool mol_is_set{false};
};

/// @brief given molecule type and number of carbon atoms, return closest number of molecules
int getMolNum(const std::string& mol, int atomN);
