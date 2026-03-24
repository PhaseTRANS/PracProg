#include<cmath>
#include<functional>
#include<vector>
#include<algorithm>
#include<fstream>
#include"matrix.hpp"

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

int main() {
    //----------PART A------------
    std::cout << "--------------PART A---------------" << std::endl;

    // example from https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
    auto F = [](double x, pp::vector y) -> pp::vector {
        double b = 0.25;
        double c = 5.0;
        pp::vector dydt(2);
        dydt[0] =  y[1];
        dydt[1] = -b * y[1] - c * std::sin(y[0]);
        return dydt;
    };

    // initial conditions
    pp::vector yinit = {M_PI - 0.1, 0.0};
    auto [xs, ys] = driver(F, {0.0, 10}, yinit);

    // build result file
    std::ofstream ode_result("ode_result.dat");
    for (int i = 0; i < (int)xs.size(); i++) {
        double x = xs[i];
        double theta = ys[i][0];
        double omega = ys[i][1];
        ode_result << x << " " << theta << " " << omega << std::endl;
    }
    ode_result.close();

    //-------------PART B------------
    std::cout << "-----------PART B-------------" << std::endl;

    // function for the orbit equation
    auto makeOrbit = [](double eps) {
        return [eps](double phi, pp::vector y) -> pp::vector {
            pp::vector dydt(2);
            dydt[0] =  y[1];
            dydt[1] = 1 - y[0] + eps * y[0] * y[0];
            return dydt;
        };
    };

    // inital condition to change in the loop
    std::vector<double> epss    = {0.0,  0.0, 0.01};
    std::vector<double> uprimes = {0.0, -0.5, -0.5};

    // vectors to store orbits
    std::vector<std::vector<double>> results_phi;
    std::vector<std::vector<pp::vector>> results_u;

    for (int i=0; i<(int)epss.size(); i++) {
        pp::vector init = {1.0, uprimes[i]};
        auto F = makeOrbit(epss[i]);
        auto [phis, us] = driver(F, {0.0, 8*M_PI}, init, 0.01, 1e-2, 1e-2);
        results_phi.push_back(phis);
        results_u.push_back(us);
    }
    
    // build orbit file
    std::ofstream orbits("orbits.dat");
    for (int i=0; i<(int)epss.size(); i++) {
        orbits << "# eps=" << epss[i] << " uprime=" << uprimes[i] << "\n";
        for (int j=0; j<(int)results_phi[i].size(); j++) {
            orbits << results_phi[i][j] << " "          // phi
                << results_u[i][j][0] << " "            // u
                << results_u[i][j][1] << std::endl;     // u'
        }
        orbits << "\n\n";  // devide the data into blocks for plotting
    }
    orbits.close();

    //-----------PART C------------
    std::cout << "------------PART C--------------" << std::endl;
    // I dont know if I want to do this

    return 0;
}