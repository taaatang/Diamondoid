#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include <vector>
#include <utility>
#include "atom.hpp"

class Molecule
{
private:
    /* data */
public:
    Molecule(/* args */);
    ~Molecule();
    void move(){cur = std::move(next); next.clear();}
    std::vector<Atom*> cur, next;
    int bondNum{0};
};

Molecule::Molecule(/* args */)
{
}

Molecule::~Molecule()
{
}
#endif // __MOLECULE_H__