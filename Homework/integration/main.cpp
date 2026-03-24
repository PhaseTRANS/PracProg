#include<iostream>
#include<cmath>
#include<functional>
#include<numeric>
#include<numbers>
#include<fstream>
#include<vector>
#include"matrix.hpp"


std::pair<double, double> integrate(std::function<double(double)> f, double a, double b,
    double acc=0.001, double eps=0.001, double f2=std::nan("1"), double f3=std::nan("1")) { // NaN indicates first call 
   
    double h=b-a;
    if(std::isnan(f2)){ 
        std::cout << "working" << std::endl;
        f2=f(a+2*h/6); f3=f(a+4*h/6);
    } // first call, no points to reuse
    double f1=f(a+h/6), f4=f(a+5*h/6);
    double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
    double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule

    double err=std::abs(Q-q);
    if(err <= acc+eps*std::abs(Q)) return {Q,err};
    else {
        auto [Q1, err1] = integrate(f,a,(a+b)/2,acc/std::sqrt(2),eps,f1,f2);
        auto [Q2, err2] = integrate(f,(a+b)/2,b,acc/std::sqrt(2),eps,f3,f4);
        return {Q1+Q2,std::sqrt(err1*err1+err2*err2)};
        }
};

// integrator transforming into clenshaw-curtis
// it uses the integrator above on the transformed function
auto integrate_cc(std::function<double(double)> f, double a, double b,
    double acc = 0.001, double eps = 0.001) {

    if (a == -INFINITY && b == INFINITY) {
        auto integrand = [&](double t) {
            return f(t / (1 - t*t)) * (1 + t*t) / ((1-t*t) * (1-t*t));
        };
        return integrate(integrand, -1.0, 1.0, acc, eps);
    }
    if (a == -INFINITY && std::isfinite(b)) {
        auto integrand = [&](double t) {
            return f(b - (1-t)/t) / (t*t);
        };
        return integrate(integrand, 0.0, 1.0, acc, eps);
    }
    if (std::isfinite(a) && b == INFINITY) {
        auto integrand = [&](double t) {
            return f(a + (1-t)/t) / (t*t);
        };
        return integrate(integrand, 0.0, 1.0, acc, eps);
    }

    // Both finite: apply Clenshaw-Curtis transformation
    double mid  = (a + b) / 2;
    double half = (b - a) / 2;
    auto integrand = [&](double theta) {
        return f(mid + half * std::cos(theta)) * std::sin(theta) * half;
    };
    return integrate(integrand, 0.0, M_PI, acc, eps);
};

double my_erf(double z, double eps=0.001, double acc=0.001) {
    if (z < 0) {
        return -my_erf(-z, eps, acc);
    }
    if (0 <= z && z <= 1) {
        auto F = [](double x) {
            return std::exp((-x) * x);
        };
        auto [val_F, err_F] =  integrate(F, 0.0, z, acc, eps);
        return 2.0 / std::sqrt(M_PI) * val_F;
    }
    if (z > 1) {
        auto T = [z](double t) {
            return std::exp(-(z+(1-t)/t) * (z+(1-t)/t))/t/t ;
        };
        auto [val_T, err_T] = integrate(T, 0.0, 1.0, acc, eps);
        return 1 - 2 / std::sqrt(M_PI) * val_T;
    }
    else return 0.0;
};

std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> result;

    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (stop - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
}

int main() {
    std::cout << "-----------PART A-------------" << std::endl;
    std::cout << "Testing integration..."  << std::endl;

    // Define functions to test on 
    auto A = [](double x) {
        return std::sqrt(x);
    };
    auto B = [](double x) {
        return 1/std::sqrt(x);
    };
    auto C = [](double x) {
        return std::sqrt(1 - x * x);
    };
    auto D = [](double x) {
        return std::log(x) / std::sqrt(x);
    };
    double lower = 0.0;
    double upper = 1.0;

    // Evaluate definite integral
    auto [int_A, err_A] = integrate(A, lower, upper);
    auto [int_B, err_B] = integrate(B, lower, upper);
    auto [int_C, err_C] = integrate(C, lower, upper);
    auto [int_D, err_D] = integrate(D, lower, upper);
    std::cout << "integrate(A, 0, 1) = " << int_A << " analytical result: " << "0.66667" << std::endl;
    std::cout << "integrate(B, 0, 1) = " << int_B << " analytical result: " << "2" << std::endl;
    std::cout << "integrate(C, 0, 1) = " << int_C << " analytical result: " << "0.7854" << std::endl;
    std::cout << "integrate(D, 0, 1) = " << int_D << " analytical result: " << "-4" << std::endl;

    std::cout << "Testing error function..."  << std::endl;
    double erf_1 = my_erf(1.0);
    double erf_2 = my_erf(2.0);
    double erf_3 = my_erf(-1.0);
    std::cout << "erf(1) = " << erf_1 << "Tabulated value: erf(1)=0.84270079" << std::endl;
    std::cout << "erf(2) = " << erf_2 << "Tabulated value: erf(2)=0.99532" << std::endl;
    std::cout << "erf(-1) = " << erf_3 << "Tabulated value: erf(-1)=-0.84270079" << std::endl;

    std::vector zs = linspace(-10, 10, 300);
    // writing error function data to fil
    std::ofstream erf_data("erf.dat");
    for (int i = 0; i < (int)zs.size(); i++) {
        double z = zs[i];
        double val = erf(z);
        erf_data << z << " " << val << std::endl;
    }
    erf_data.close();

    // making accuracy plot
    double acc = 0.1;
    double erf1_exact = 0.84270079294971486934;
    pp::vector accs(13);
    std::ofstream acc_data("acc.dat");
    for (int i=0; i<(int)accs.size(); i++) {
        accs[i] = acc;
        double value_0 = my_erf(1.0, 0.0, acc);
        acc_data << acc << " " << value_0 << " " << std::abs(value_0 - erf1_exact) << std::endl;
        acc /= 10;
    }
    acc_data.close();

    std::cout << "-----------PART B-------------" << std::endl;
    auto [intcc_A, errcc_A] = integrate_cc(A, lower, upper);
    auto [intcc_B, errcc_B] = integrate_cc(B, lower, upper);
    auto [intcc_C, errcc_C] = integrate_cc(C, lower, upper);
    auto [intcc_D, errcc_D] = integrate_cc(D, lower, upper);
    // print results and compare to normale integrator
    std::cout << "integrate_cc(A, 0, 1) = " << intcc_A << " analytical result: " << "0.66667" << " |with trans. - without trans.|= " << std::abs(intcc_A-int_A) << std::endl;
    std::cout << "integrate_cc(B, 0, 1) = " << intcc_B << " analytical result: " << "2" << " |with trans. - without trans.|= " << std::abs(intcc_B-int_B) << std::endl;
    std::cout << "integrate_cc(C, 0, 1) = " << intcc_C << " analytical result: " << "0.7854" << " |with trans. - without trans.|= " << std::abs(intcc_C-int_C) << std::endl;
    std::cout << "integrate_cc(D, 0, 1) = " << intcc_D << " analytical result: " << "-4" << " |with trans. - without trans.|= " << std::abs(intcc_D-int_D) << std::endl;

    // counting the integration calls
    int ncalls = 0;
    // function which counts its evaluations
    auto f = [&ncalls](double z) {
        ncalls++;
        return z * z;
    };
    auto [int_f, err_f] = integrate_cc(f, lower, upper);
    std::cout << "integrate_cc(z^2, 0, 1) = " << int_f << " Evaluations on integrand: " << ncalls << " Evaluations using scipy: 21" << std::endl;

    // Evaluate infinite limits implementation
    auto g = [&](double x) {
        ncalls++;
        return std::exp(- x * x);
    };
    double upper_inf = INFINITY;
    double lower_inf = -INFINITY;
    ncalls = 0;
    auto [int_g, err_g] = integrate_cc(g, lower_inf, upper_inf);
    std::cout << "integrate_cc(exp(-x^2), -inf, inf) = " << int_g << " Evaluations on integrand: " << ncalls << " Evaluations using scipy: 277" << std::endl;
    ncalls = 0;
    auto [int1_g, err1_g] = integrate_cc(g, lower, upper_inf);
    std::cout << "integrate_cc(exp(-x^2), 0, inf) = " << int1_g << " Evaluations on integrand: " << ncalls << " Evaluations using scipy: 135" << std::endl;
    ncalls = 0;
    auto [int2_g, err2_g] = integrate_cc(g, lower_inf, lower);
    std::cout << "integrate_cc(exp(-x^2), -inf, 0) = " << int2_g << " Evaluations on integrand: " << ncalls << " Evaluations using scipy: 135" << std::endl;

    std::cout << "-----------PART C-------------" << std::endl;
    // first complicated integral \int_0^pi sin(x)/exp(-x^2) dx
    auto c = [](double x) {
        return std::sin(x)/std::exp(- x * x);
    };
    auto [int_c, err_c] = integrate_cc(c, 0.0, M_PI);
    std::cout << "int_0^pi sin(x)/exp(-x^2) = " << int_c << " integration error: " << err_c << " actual error: " << std::abs(570.4357348443238 - int_c) << std::endl;

    // second complicated integral \int_0^inf x/exp(x) dx
    auto d = [](double x) {
        return x/std::exp(x);
    };
    auto [int_d, err_d] = integrate_cc(d, 0.0, INFINITY);
    std::cout << "int_0^inf x/exp(x) = " << int_d << " integration error: " << err_d << " actual error: " << std::abs(1-int_d) << std::endl;

    return 0;
}