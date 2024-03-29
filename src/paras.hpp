#pragma once

#include <iostream>
#include <memory>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <type_traits>
#include <assert.h>

template<typename T>
using map = std::map<std::string,T>;

/**
 * @brief Parameter class. Handling the input
 * 
 */
class Parameters{

public:
    Parameters(){}
    Parameters(std::string configFile){config(configFile);}
    Parameters(std::string InputDir, std::vector<std::string> files);
    ~Parameters(){};
    void config(std::string configFile);
    void clear();
    bool read(const std::string& filename);
    void print(std::string filename) const;
    void print(std::ostream& os) const;

    std::string rootDataPath,project;
    map<int> mapi;
    map<double> mapd;
    map<std::string> maps;
    map<std::vector<std::string>> mapvecs;
    map<std::vector<double>> mapvecd;
    map<std::vector<std::vector<double>>> maparrd;
};

inline void removeComment(std::string& s, char delim){
    std::istringstream ins (s);
    std::string tmp;
    std::getline(ins,tmp,delim);
    s = std::move(tmp);
}

template<typename T>
std::vector<T> readVec(std::string& s, char sep=','){
    std::istringstream ins(s);
    std::vector<T> result;
    std::string s_val;
    while(std::getline(ins,s_val,sep)){
        T val;std::istringstream ins_val(s_val);ins_val>>val;
        result.push_back(val);
    }
    return result;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& myvec){
    if(!std::is_same<T,std::vector<double>>::value)os<<"\n";
    if(std::is_floating_point<T>::value)os<<std::setprecision(2)<<std::fixed;
    for(int i = 0; i< (int)myvec.size()-1; ++i){
        os<<myvec.at(i);
        if(!std::is_same<T,std::vector<double>>::value)os<<", ";
    }
    if(!myvec.empty())os<<myvec.back();
    if(std::is_floating_point<T>::value)os.unsetf(std::ios::floatfield); 
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const map<T>& mymap){
    for(const auto& pair:mymap){
        os<<pair.first<<":"<<pair.second<<"\n";
    }
    return os;
}