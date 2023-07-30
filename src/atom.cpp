#include "atom.hpp"

Atom::Atom(Type type_in, std::array<int,3> coord_in, int idx_in): idx(idx_in), type(type_in) {
	coord = coord_in; 
	for (auto& el:NN) {
		el=nullptr;
	}
}

void Atom::classNN(std::vector<int>& filled, std::vector<int>& empty) {
	filled.clear();
	empty.clear();
	for (int i = 0; i < 4; i++) {
		if (NN[i]) {
			if (NN[i]->idx == -1) empty.push_back(i);
			else filled.push_back(i);
		}
	}
}

int Atom::countFilledNN() const {
	int count = 0; 
	for (const auto& nn:NN) {
		if (nn and nn->idx != -1) count++;
	} 
	return count;
};

void unvisit(std::vector<Atom>* mol,  const VecI &visited_idx){
    for(auto idx:visited_idx){
        mol->at(idx).is_visited=false;
    }
}

void unvisit(std::vector<Atom>* mol){
    for(auto& atom:(*mol)){
        atom.is_visited=false;
    }
}

bool operator<(const Atom& lhs, const Atom& rhs){
    return lhs.coord<rhs.coord;
}

bool operator>(const Atom& lhs, const Atom& rhs){
    return lhs.coord>rhs.coord;
}

bool operator==(const Atom& lhs, const Atom& rhs){
    return lhs.coord==rhs.coord;
}