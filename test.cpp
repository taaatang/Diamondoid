#include <iostream>
#include "src/algebra.hpp"
#include "src/randoms.hpp"
#include "src/MD.hpp"
int main(){
    // std::cout<<"Mol size:"<<Adamantane.size()<<"\n";
    // for(auto a:Adamantane)std::cout<<a;
    // ConstructMol();
    // std::cout<<"IMat:"<<ComputeInert();
    Initialization();
    while(NotFinished()){
        SingleStep();
    }
    return 0;
}