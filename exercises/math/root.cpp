#include<cmath>
#include<iostream>

int main(int argc, char** argv) {
    int n = std::stoi(argv[1]);
    double root2 = std::sqrt(n);
    std::cout << "sqrt(n) =" << root2 << std::endl;
    return 0;
}