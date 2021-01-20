#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <vector>
#include "atom.hpp"

struct Lattice{
    Lattice(){}
    ~Lattice(){}
    int Na{0}, Nb{0}, Nc{0};
    Atom* center;
    std::vector<Atom> latt;
}

#endif // __LATTICE_H__