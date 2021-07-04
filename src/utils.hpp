#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <chrono>
#include <filesystem>
#include <iomanip>
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
    double elapse(){
        if (!is_tok) tok();
        duration_ = tok_ - tik_;
        is_tok = false;
        return duration_.count()*1000.0;
    }
};

inline void assert_msg(bool condition, const std::string& msg){
    if(!condition){
        std::cout<<msg<<std::endl;
        exit(1);
    }
}

inline void mkdir_fs(const std::string& dir) {
    std::filesystem::path p(dir);
    std::filesystem::create_directories(p);
    bool succeed = std::filesystem::is_directory(p);
    assert_msg(succeed, dir + " failed to creat!");
}

template <class T>
inline void save(T *d_pt, int size, std::ofstream *f_pt, const std::string& filename, bool is_app=false){
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
#endif // __UTILS_H__