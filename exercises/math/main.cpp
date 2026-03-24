#include <iostream>
#include "sfuns.h"
#include<complex>

using complex=std::complex<double>;
constexpr double  π = 3.14159265358979324;
constexpr double  E = 2.71828182845904523;
constexpr complex I = complex(0,1);

int main() {
    std::cout << "--------------PART 1-----------------" << std::endl;
    std::cout << "log(I)=" << std::log(I)   << " Exact: (0, 1.5707963267948966)" << "\n";
    std::cout << "   I^I=" << std::pow(I,I) << " Exact: (0.20787957635076193, 0)" << "\n";
    std::cout << "   π^I=" << std::pow(π,I) << " Exact: (0.41329211610159433, 0.9105984992126147)" << "\n";
    std::cout << "   E^I=" << std::pow(E,I) << " Exact: (0.5403023058681398, 0.8414709848078965)" << "\n";
    std::cout << "   sqrt(2)=" << std::sqrt(2) << " Exact: 1.4142135623730951" << "\n";
    std::cout << "   2^(1/5)=" << std::pow(2, 1/5) << " Exact: 1.148698354997035" << "\n";
    std::cout << "   e^pi=" << std::pow(E, π) << " Exact: 23.140692632779267" << "\n";
    std::cout << "   pi^e=" << std::pow(π, E) << " Exact: 22.45915771836104" << "\n";


    std::cout << "--------------PART 2-----------------" << std::endl;

    for (int i = 1; i < 10 + 1; i++) {
        std::cout << "loggamma(" << i << ")=" << sfuns::lngamma(i) << '\n';
    }

    std::cout << "-----------------Exact results--------------------" << std::endl;
    std::cout << "loggamma(1) = " << "0.000000" << std::endl;
    std::cout << "loggamma(2) = " << "0.000000" << std::endl;
    std::cout << "loggamma(3) = " << "0.693147" << std::endl;
    std::cout << "loggamma(4) = " << "1.791759" << std::endl;
    std::cout << "loggamma(5) = " << "3.178053" << std::endl;
    std::cout << "loggamma(6) = " << "4.787491" << std::endl;
    std::cout << "loggamma(7) = " << "6.579251" << std::endl;
    std::cout << "loggamma(8) = " << "8.525161" << std::endl;
    std::cout << "loggamma(9) = " << "10.60460" << std::endl;
    std::cout << "loggamma(10) = " << "12.80182" << std::endl;


}

