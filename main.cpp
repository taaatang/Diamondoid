#include <iostream>
#include <vector>
#include <string>
#include "utils/utils.hpp"
#include "src/random.hpp"
#include "src/molecule.hpp"
#include "src/cluster.hpp"

const int HEAT_STEP = 100000;
const int COOL_STEP = 800000;
const int MAX_MC_STEP = 100000;
const int TOT_STEP = HEAT_STEP + COOL_STEP + MAX_MC_STEP;

constexpr int MEASURE_STEP = 100;
const int PRINT_STEP = (HEAT_STEP + COOL_STEP + MAX_MC_STEP)/10;

const double INIT_TEMP = 3.0;
const double FINAL_TEMP = 0.1;
const int TAdjust_STEP = 1000;

const int N_mol = 39;

int main() {
    Timer timer;
    std::string dir = "data/Pentamantane_"+std::to_string(N_mol);
    mkdir_fs(dir);

    // // linear temperature
    // auto temp = [=](int step) {
    //     return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) / MAX_MC_STEP * step;
    // };

    // exp temperature
    auto temp = [=](int step) {
        if (step < HEAT_STEP) return INIT_TEMP;
        else if (step > HEAT_STEP+COOL_STEP) return FINAL_TEMP;
        else return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) / COOL_STEP * (step-HEAT_STEP); 
        // else return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) * (1.0 - std::exp((double)(step-HEAT_STEP)/COOL_STEP))/(1.0 - std::exp(1.0));
    };

    std::vector<Atom> Mol;
    timer.tik();
    std::vector<Pentamantane> Mols(N_mol,Pentamantane(0,&Mol));
    std::vector<Molecule*> Molsptr; for(auto& mol:Mols)Molsptr.push_back(&mol);

    // INIT
    Cluster clus(40,40,40);
    clus.setTemp(temp(0));
    clus.init(Molsptr);
    timer.tok();
    std::cout<<"Initialization time:"<<timer.elapse()/1000.0<<"s.\n\n";

    // // HEAT
    // timer.tik();
    // for (int count = 0; count < HEAT_STEP; ++count) {
    //     clus.singleStep();
    // }
    // timer.tok();
    // std::cout<<"Heating time:"<<timer.elapse()/1000.0<<"s.\n\n";

    // MC
    timer.tik();
    std::vector<int> bondNum;
    for (int stepCount = 1; stepCount <= TOT_STEP; ++stepCount){

        clus.singleStep();

        if(stepCount%MEASURE_STEP==0){
            int mstep = stepCount / MEASURE_STEP;
            std::string outfile = dir + "/step_"+std::to_string(mstep)+".dat";
            auto bnum = clus.countBond();
            bondNum.push_back(bnum);
            clus.saveCoords(outfile);
        }

        if(stepCount%PRINT_STEP==0){
            auto bnum = clus.countBond();
            std::cout<<"Step: "<<stepCount<<" / "<<TOT_STEP<<", T:"<<temp(stepCount)<<", Total Bonds: "<<bnum<<".\n";
        }

        if(stepCount%TAdjust_STEP==0) {
            clus.setTemp(temp(stepCount));
        }

    }

    // save measured bond num
    std::ofstream outf;
    save<int>(bondNum.data(),bondNum.size(),&outf,dir+"/bondNum.dat");
    timer.tok();
    std::cout<<"\nMonte Carlo time:"<<timer.elapse()/1000.0<<"s.\n";

    // evaluate compatible/uncompatible surface positions
    timer.tik();
    std::cout<<"\nBegin evaluate surface positions...\n";
    clus.computePos();
    clus.evalPos();
    clus.saveSurface(dir + "/surface");
    std::cout<<"Surface evaluation time:"<<timer.elapse()/1000.0<<"s.\n";
    timer.tok();

    // clustering occupied/unoccupied sites. discovering vacancy
    timer.tik();
    clus.clustering();
    clus.saveVacancy(dir + "/vacancy.dat");
    timer.tok();
    std::cout<<"Clustering time:"<<timer.elapse()/1000.0<<"s.\n";

  
    
    // clus.saveCoords(outfile);
    return 0;
}