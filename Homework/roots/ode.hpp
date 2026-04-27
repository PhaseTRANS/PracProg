#include<functional>
#include<vector>
#include"matrix.hpp"


namespace pp{
    // rkstep12: single Runge-Kutta step with error estimate
std::pair<pp::vector, pp::vector> rkstep12(
    std::function<pp::vector(double, pp::vector)> f,    /* the f from dy/dx=f(x,y) */
    double x,                                           /* current value of the variable */
    pp::vector y,                                       /* current value y(x) */
    double h                                            /* step size */
)
{
    pp::vector k0 = f(x, y);                            /* embedded lower order formula (Euler) */
    pp::vector k1 = f(x + h/2, y + k0*(h/2));           /* higher order formula (midpoint) */
    pp::vector yh = y + k1*h;                           /* y(x+h) estimate */
    pp::vector dy = (k1 - k0)*h;                        /* error estimate */
    return {yh, dy};
}

// driver: adaptive step-size ODE integrator
std::pair<std::vector<double>, std::vector<pp::vector>> driver(
    std::function<pp::vector(double, pp::vector)> F,    /* the f from dy/dx=f(x,y) */
    std::pair<double, double> interval,                 /* (initial-point, final-point) */
    pp::vector yinit,                                   /* y(initial-point) */
    double h   = 0.125,                                 /* initial step-size */
    double acc = 0.01,                                  /* absolute accuracy goal */
    double eps = 0.01                                   /* relative accuracy goal */
)
{
    auto [a, b] = interval;
    double x = a;
    pp::vector y = yinit;

    std::vector<double> xlist;  xlist.push_back(x);
    std::vector<pp::vector> ylist;  ylist.push_back(y);

    do {
        if (x >= b) return {xlist, ylist};    /* job done */
        if (x + h > b) h = b - x;             /* last step should end at b */

        auto [yh, dy] = rkstep12(F, x, y, h);

        double tol = (acc + eps * yh.norm()) * std::sqrt(h / (b - a));
        double err = dy.norm();

        if (err <= tol) {   // accept step
            x += h;  y = yh;
            xlist.push_back(x);
            ylist.push_back(y);
        }

        if (err > 0) h *= std::min(std::pow(tol / err, 0.25) * 0.95, 2.0);
        else         h *= 2;

    } while (true);
}

}