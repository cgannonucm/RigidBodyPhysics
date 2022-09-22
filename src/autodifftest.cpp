/*
#include <iostream>
#include <autodiff/reverse/var.hpp>
#include <Eigen/Eigen/Dense>

using namespace autodiff;

double partial(var const &u, var const &x)
{
    auto l = derivatives(u,wrt(x));
    std :: cout << l[0] << std::endl;
    return l[0];
}

var f(var &x, var &y, var &z){
    return x + y + z;
}

int main() {

    long double y_d = 0;

    var x = y_d;
    var y = 0;
    var z = 0;

    var u = f(x,y,z);

    std :: vector<var *> args = {&x,&y,&z};

    for (var * arg : args){
        auto d = partial(u,*arg);
        std :: cout << d << std::endl;
    }

    auto l = partial(u,x); 

    //std :: cout << l[0] << std::endl;
    std :: cout << "Test" << std::endl;
    return 0;
};*/