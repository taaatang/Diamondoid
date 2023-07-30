#pragma once

#include <vector>
#include <algorithm>
#include <queue>
#include <cassert>
#include "atom.hpp"
#include "algebra.hpp"

struct Lattice {
    Lattice(int na, int nb, int nc);

    void linkNN(Atom &atom, const std::array<std::array<int, 3>, 4> &NN);

    void linkNN();

    int Na{0}, Nb{0}, Nc{0};
    Atom *center;
    std::vector<Atom> latt;
};