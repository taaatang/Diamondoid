#include <iostream>
#include "src/algebra.hpp"
#include "src/randoms.hpp"
#include "src/MD.hpp"
extern int stepCount;
extern const int tstepNum = 10000;
int main(){
    // std::cout<<"Mol size:"<<Adamantane.size()<<"\n";
    // for(auto a:Adamantane)std::cout<<a;
    // std::cout<<Adamantane[0]+Adamantane[1];
    // std::cout<<norm(Adamantane[0])<<"\n";
    // std::cout<<"Random:";
    // for(int i = 0; i < 10; i++) std::cout<<RandR()<<" ";
    // std::cout<<"\nRandom vec:"<<VRand()<<"\n";
    // ConstructMol();
    // std::cout<<"IMat:"<<ComputeInert();
    Initialization();
    while(stepCount < tstepNum){
        SingleStep();
        stepCount++;
    }
    // Vec3d a = {1,2,3};
    // Vec3d b = a;
    // std::cout<<a-b;
    return 0;
}