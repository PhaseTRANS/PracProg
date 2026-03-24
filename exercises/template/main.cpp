#include<iostream>
#include"vec.h"
#include<complex>

template<typename T>
void test_vec(const std::string& label) {  // function which tests the vec class for the given type
    std::cout << "-----------" << label << "--------------\n";

    T x1 = T(1);
    T y1 = T(2);
    T z1 = T(3);

    T x2 = T(5);
    T y2 = T(-5);
    T z2 = T(15);

    vec<T> v{x1,y1,z1};
    vec<T> w{x2,y2,z2};

    std::cout << "----Print vector----\n";
    v.print("v = ");
    w.print("w = ");
    std::cout << "<< : " << v << "\n";
    std::cout << "<< : " << w << "\n";

    std::cout << "----Add----\n";
    auto a = v + w;
    a.print("v + w = ");

    std::cout << "----Subtract----\n";
    auto b = v - w;
    b.print("v - w = ");

    std::cout << "----Multiply----\n";
    auto c = T(2) * v;
    c.print("2 * v = ");

    std::cout << "----Norm----\n";
    std::cout << "norm = " << v.norm() << "\n";

    std::cout << "----Cross-product----\n";
    auto x = v.cross(w);
    x.print("v x w = ");

    std::cout << "----Dot-product----\n";
    auto y = v.dot(w);
    std::cout << "v * w = " << y << "\n\n";
}


int main() {
    // Testing for all scalar types that support arithmetic operations
    test_vec<double>("Double vector");
    test_vec<long double>("Long Double vector");
    test_vec<float>("Float vector");
    test_vec<int>("Int vector");
    test_vec<std::complex<double>>("Complex vector");

    return 0;
}