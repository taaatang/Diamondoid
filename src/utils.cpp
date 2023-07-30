#include <sstream>
#include <filesystem>
#include <iomanip>

#include "utils.hpp"

double Timer::elapse() {
	if (!is_tok) tok();
	duration_ = tok_ - tik_;
	is_tok = false;
	return duration_.count()*1000.0;
}

std::string tostr(double val, int digit){
    std::ostringstream strTmp;
    strTmp<<std::fixed<<std::setprecision(digit)<<val;
    return strTmp.str();
}

std::string tostr(int val){
    return std::to_string(val);
}

void assert_msg(bool condition, const std::string& msg){
    if(!condition){
        std::cout<<msg<<std::endl;
        exit(1);
    }
}

void mkdir_fs(const std::string& dir) {
    std::filesystem::path p(dir);
    std::filesystem::create_directories(p);
    bool succeed = std::filesystem::is_directory(p);
    assert_msg(succeed, dir + " failed to creat!");
}