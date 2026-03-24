#include<iostream>
#include<limits>
#include<cmath>
#include<iomanip>

// function to compare double by machien precision
bool approx(double a, double b, double acc=1e-9, double eps=1e-9) {
    double absolute = std::abs(a - b);
    double relative = absolute / std::max(std::abs(a), std::abs(b));
    if (absolute <= acc || relative <= eps) {
        return true;
    }
    else {
        return false;
    }
}

int main() {
    // PART 1
    std::cout << "-----------PART 1------------" << std::endl;

    float f=1.0; 
    while(float (1.0+f) != 1.0){
        f/=2.0;
    } f*=2.0;

    double d=1.0; 
    while(double (1.0+d) != 1.0) {
        d/=2.0;
    } d*=2.0;

    long double l=1.0; 
    while((long double) (1.0+l) != 1.0) {
        l/=2.0;
    } l*=2.0;

    // printing simulated epsilon
    std::printf("float eps=%g\n",f);
    std::printf("double eps=%g\n",d);
    std::printf("long double eps=%Lg\n",l);

    // printing computer epsilon
    std::cout << std::numeric_limits<float>::epsilon() << "\n";
    std::cout << std::numeric_limits<double>::epsilon() << "\n";
    std::cout << std::numeric_limits<long double>::epsilon() << "\n";

    // calculating machine epsilon
    std::cout << "Calculated epsilon double (2^-52)=" << std::pow(2, -52) << std::endl;
    std::cout << "Calculated epsilon float (2^-23)=" << std::pow(2, -23) << std::endl;

    // PART 2
    std::cout << "-----------PART 2------------" << std::endl;

    double epsilon = std::pow(2, -52);
    double tiny = epsilon / 2;
    double a = 1 + tiny + tiny;
    double b = tiny + tiny + 1;

    std::cout << "a==b ? " << (a==b ? "true":"false") << std::endl;
    std::cout << "a>1  ? " << (a>1  ? "true":"false") << std::endl;
    std::cout << "b>1  ? " << (b>1  ? "true":"false") << std::endl;

    std::cout << std::fixed << std::setprecision(17);
    std::cout << "tiny=" << tiny << std::endl;
    std::cout << "1+tiny+tiny=" << a << std::endl;
    std::cout << "tiny+tiny+1=" << b << std::endl;

    // PART 3
    std::cout << "-----------PART 3------------" << std::endl;

    double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
    double d2 = 8*0.1;
    std::cout << "d1==d2? " << (d1==d2 ? "true":"false") << std::endl; 
    std::cout << std::fixed << std::setprecision(17);
    std::cout << "d1=" << d1 << std::endl;
    std::cout << "d2=" << d2 << std::endl;

    std::cout << "d1==d2? " << approx(d1, d2) << std::endl; 
    return 0;
}