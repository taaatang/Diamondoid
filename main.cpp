#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "src/cluster.hpp"
#include "src/paras.hpp"

int main() {
    Parameters para(".", {"input.txt"});

    const std::string dataDir = para.maps.at("data dir");
    const std::string molName = para.maps.at("molecule name");
    const std::string runid = para.maps.at("runid");
    //int molNum = para.mapi.at("molecule num");
    const int atomNum = para.mapi.at("atom num");
    const bool isRandseed = para.mapi.at("random seed");
    const int HEAT_STEP = para.mapi.at("heat step");
    const int ANNEAL_STEP = para.mapi.at("anneal step");
    const int COOL_STEP = para.mapi.at("cool step");
    const int MEASURE_STEP = para.mapi.at("measure step");
    const int TAdjust_STEP = para.mapi.at("TAdjust step");
    const double INIT_TEMP = para.mapd.at("init T");
    const double FINAL_TEMP = para.mapd.at("final T");
    const int TOT_STEP = HEAT_STEP + ANNEAL_STEP + COOL_STEP;
    const int PRINT_STEP = (HEAT_STEP + ANNEAL_STEP + COOL_STEP)/10;
    const double PRESSURE = para.mapd.at("pressure");
    const double JUMP_LIMIT = para.mapd.at("jump limit");
    // temperature
    auto temp = [=](int step) {
        if (step < HEAT_STEP) return INIT_TEMP;
        else if (step >= HEAT_STEP + ANNEAL_STEP) return FINAL_TEMP;
        else return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) / ANNEAL_STEP * (step - HEAT_STEP);
        // else return INIT_TEMP + (FINAL_TEMP - INIT_TEMP) * (1.0 - std::exp((double)(step-HEAT_STEP)/COOL_STEP))/(1.0 - std::exp(1.0));
    };

    auto measureBond = [=](int step) {
        return step % MEASURE_STEP == 0;
    };

    auto measureStruc = [=](int step) {
        int equil_step = HEAT_STEP + ANNEAL_STEP;
        int m1 = equil_step / 2;
        int m2 = COOL_STEP / 40;
        if (step <= equil_step) {
            return step % m1 == 0;
        } else {
            return (step - equil_step) % m2 == 0;
        }
    };

    auto molNum = getMolNum(molName, atomNum);
    std::cout << "Begin simulation for:" << molNum << " x " << molName + "mantane" + " ~ " << atomNum << " x "
              << "C atom" << std::endl;
    int digit = 4;
    auto TPD_LABEL =
            "T_" + tostr(FINAL_TEMP, digit) + "/P_" + tostr(PRESSURE, digit) + "/D_" + tostr(JUMP_LIMIT, digit);
    std::cout << TPD_LABEL << std::endl;
    Timer timer;
    std::string dir =
            dataDir + "/" + molName + "mantane_" + std::to_string(molNum) + "/" + TPD_LABEL + "/run" + runid;
    mkdir_fs(dir);
    para.print(dir + "/para.txt");
    randomSeed(isRandseed);


    std::vector<Atom> Mol;
    timer.tik();
    std::vector<Molecule *> Molsptr;
    for (int i = 0; i < molNum; ++i) {
        if (molName == "Ada") {
            Molsptr.push_back(new Adamantane(0, &Mol));
        } else if (molName == "Dia") {
            Molsptr.push_back(new Diamantane(0, &Mol));
        } else if (molName == "Tria") {
            Molsptr.push_back(new Triamantane(0, &Mol));
        } else if (molName == "Tetra") {
            Molsptr.push_back(new Tetramantane(0, &Mol));
        } else if (molName == "Penta") {
            Molsptr.push_back(new Pentamantane(0, &Mol));
        } else if (molName == "Penta1212") {
            Molsptr.push_back(new Pentamantane1212(0, &Mol));
        } else {
            std::cout << "molecule " << molName << "mantane not defined!\n";
            exit(1);
        }
    }
    std::cout << molNum << " " << molName << "mantanes are created!\n";

    // INIT
    Cluster clus(40, 40, 40);
    std::cout << "init temp = " << temp(0) << std::endl;
    clus.setTemp(temp(0));
    clus.init(Molsptr);
    timer.tok();
    std::cout << "Initialization time:" << timer.elapse() / 1000.0 << "s.\n\n";
    int mstep = 0;
    {// save initial configuration
        std::string outfile = dir + "/step_" + std::to_string(mstep) + ".dat";
        clus.saveCoords(outfile);
    }
    // MC
    timer.tik();
    std::vector<int> bondNum;
    clus.setPressure(PRESSURE);
    double dMax = JUMP_LIMIT * std::sqrt(10.0 / Mol.size());
    std::cout << "Jump limit: " << dMax << std::endl;
    for (int stepCount = 1; stepCount <= TOT_STEP; ++stepCount) {

        clus.singleStep(dMax);

        if (measureBond(stepCount)) {
            auto bnum = clus.countBond();
            bondNum.push_back(bnum);
        }
        if (measureStruc(stepCount)) {
            ++mstep;
            std::string outfile = dir + "/step_" + std::to_string(mstep) + ".dat";
            clus.saveCoords(outfile);

            clus.computePos();
            clus.evalPos();
            clus.saveSurface(dir + "/step_" + std::to_string(mstep) + "_surface");

            clus.clustering();
            clus.saveVacancy(dir + "/step_" + std::to_string(mstep) + "_vacancy.dat");
        }

        if (stepCount % PRINT_STEP == 0) {
            auto bnum = clus.countBond();
            std::cout << "Step: " << stepCount << " / " << TOT_STEP << ", T:" << temp(stepCount)
                      << ", Total Bonds: " << bnum << ". Total Volume: " << clus.getVol() << "\n";
        }

        if (stepCount % TAdjust_STEP == 0) {
            clus.setTemp(temp(stepCount));
        }

    }

    // save measured bond num
    std::ofstream outf;
    save<int>(bondNum.data(), bondNum.size(), &outf, dir + "/bondNum.dat");
    timer.tok();
    std::cout << "\nMonte Carlo time:" << timer.elapse() / 1000.0 << "s.\n";

    // evaluate compatible/uncompatible surface positions
    // timer.tik();
    // std::cout<<"\nBegin evaluate surface positions...\n";
    // clus.computePos();
    // clus.evalPos();
    // clus.saveSurface(dir + "/surface");
    // std::cout<<"Surface evaluation time:"<<timer.elapse()/1000.0<<"s.\n";
    // timer.tok();

    // clustering occupied/unoccupied sites. discovering vacancy
    // timer.tik();
    // clus.clustering();
    // clus.saveVacancy(dir + "/vacancy.dat");
    // timer.tok();
    // std::cout<<"Clustering time:"<<timer.elapse()/1000.0<<"s.\n";
    // clus.saveCoords(outfile);

    for (auto &p : Molsptr) {
        delete p;
    }
    return 0;
}
