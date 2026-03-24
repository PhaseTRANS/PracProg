#include<iostream>
#include"vec.h"

int main() {
    vec v = vec(1, 2, 3);
    vec w = vec(0.5, -0.5, 1.5);

    std::cout << "----Print vector----" << std::endl;
    v.print("v = ");
    w.print("w = ");
    std::cout << "Overloaded <<:" << v << std::endl;
    std::cout << "Overloaded <<:" << w << std::endl;

    std::cout << "----Add----" << std::endl;
    vec a = v + w;
    a.print("v + w = ");

    std::cout << "----subtract----" << std::endl;
    vec b = v - w;
    b.print("v - w = ");

    std::cout << "----multiply----" << std::endl;
    vec c = 2 * v;
    c.print("2 * v = ");

    std::cout << "----norm----" << std::endl;
    std::cout << "norm = " << v.norm() << std::endl;

    std::cout << "----cross-product----" << std::endl;
    vec x = v.cross(w);
    x.print("v x w = ");

    std::cout << "----dot-product----" << std::endl;
    double y = v.dot(w);
    std::cout << "v * w = " << y << std::endl;
    return 0;
}