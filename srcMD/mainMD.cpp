#include <iostream>
#include "algebra.hpp"
#include "randoms.hpp"
#include "MD.hpp"
int main(){
    // std::cout<<"Mol size:"<<Adamantane.size()<<"\n";
    // for(auto a:Adamantane)std::cout<<a;
    // ConstructMol();
    // std::cout<<"IMat:"<<ComputeInert();
    Initialization();
    while(NotFinished()){
        SingleStep();
    }
    InitMeasure();
    while(NotFinishedMeasure()){
        SingleStep();
        if(isMeasureRdf())EvalRdf();
    }
    saveRdf();
    return 0;
}