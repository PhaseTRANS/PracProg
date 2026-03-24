#include<cmath>
#include<iostream>

int main(int argc, char** argv) {
    double base = std::stoi(argv[1]);
    double exp = std::stoi(argv[2]);

    double res = pow(base, exp);
    std::cout << "base^pow=" << res << std::endl; 
    return 0;
}