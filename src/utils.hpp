#pragma once

#include <iostream>
#include <chrono>
#include <fstream>
#include <cassert>

class Timer{
    std::chrono::system_clock::time_point tik_{}, tok_{};
    std::chrono::duration<double> duration_{};
    bool is_tok{false};
public:
    Timer(){tik();}

    void tik(){tik_ = std::chrono::system_clock::now(); is_tok = false;}
    void tok(){tok_ = std::chrono::system_clock::now(); is_tok = true;}
    // return duration in unit of ms
    double elapse();
};

std::string tostr(double val, int digit);

std::string tostr(int val);

void assert_msg(bool condition, const std::string& msg);

void mkdir_fs(const std::string& dir);

template <class T>
inline void save(T *d_pt, size_t size, std::ofstream *f_pt, const std::string& filename, bool is_app=false){
    if(is_app)f_pt->open(filename, std::ios::binary|std::ios::app);
    else f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->write(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        // if(is_app)std::cout<<"Data appended to "<<filename<<std::endl;
        // else std::cout<<"Data wrote to "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}