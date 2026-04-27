#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<functional>
#include"matrix.hpp"
#include"qr.hpp"
#include"ode.hpp"

pp::matrix& jacobian(  // to allocate just one matrix and update it we use a reference instead (idea from AI)
    std::function<pp::vector(pp::vector)> f,
    pp::vector x,
    pp::vector fx,
    pp::vector dx,
    pp::matrix& J) {
        if(dx.size() == 0) dx = x.map([](double xi) { double eps = std::max(std::abs(xi), 1.0); return eps * std::pow(2, -26); });
        if(fx.size() == 0) fx = f(x);
        for(int j=0; j < x.size(); j++){
            x[j] += dx[j];
            pp::vector df = f(x) - fx;
            for(int i=0; i < x.size(); i++) J[i,j] = df[i] / dx[j];
            x[j] -= dx[j];
            }
        return J;
}

pp::vector newton(
    std::function<pp::vector(pp::vector)> f, /* the function to find the root of */
	pp::vector start,        /* the start point */
	double acc=1e-2,     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	pp::vector dx=pp::vector(),      /* optional δx-vector for calculation of jacobian */
	double lambmin=0.01) {
    pp::vector x(start);
    pp::vector fx=f(x), z, fz;
    pp::matrix J(start.size(), start.size());  // create matrix once
    do{ /* Newton's iterations */
        if(fx.norm() < acc) break; /* job done */
        jacobian(f,x,fx,dx,J);  // update it
        pp::QR QRJ;
        auto [Q, R] = QRJ.decomp(J);
        pp::vector Dx = QRJ.solve(Q, R, -fx); /* Newton's step */
        double lamb=1;
        do{ /* linesearch */
            z=x+lamb*Dx;
            fz=f(z);
            if( fz.norm() < (1-lamb/2)*fx.norm() ) break;
            if( lamb < lambmin ) break;
            lamb/=2;
            }while(true);
        x=z; fx=fz;
        }while(true);
    return x;
}

pp::vector newton_interp(
    std::function<pp::vector(pp::vector)> f,
    pp::vector start,
    double acc     = 1e-2,
    pp::vector dx  = pp::vector(),
    double lambmin = 1.0/128) {
    pp::vector x  = start;
    pp::vector fx = f(x);
    pp::matrix J(start.size(), start.size());  // create matrix once

    do {
        if (fx.norm() < acc) break;

        jacobian(f, x, fx, dx, J);
        pp::QR     QRJ;
        auto [Q, R]   = QRJ.decomp(J);
        pp::vector Dx = QRJ.solve(Q, R, -fx);

        double phi0  =  0.5 * fx.norm() * fx.norm();
        double dphi0 = -fx.norm() * fx.norm();        

        // search function
        std::function<std::pair<double, pp::vector>(double, pp::vector)>
        search = [&](double lamb, pp::vector fz) -> std::pair<double, pp::vector>
        {
            double phi_lamb = 0.5 * fz.norm() * fz.norm();

            if (phi_lamb < phi0 * (1.0 - lamb / 2.0)) return {lamb, fz};
            if (lamb < lambmin) return {lamb, fz};

            double c = (phi_lamb - phi0 - dphi0 * lamb) / (lamb * lamb);
            double lamb_new = -dphi0 / (2.0 * c);

            lamb_new = std::max(lambmin, std::min(0.9 * lamb, lamb_new));

            return search(lamb_new, f(x + lamb_new * Dx)); // recurse
        };

        auto [lamb, fz] = search(1.0, f(x + Dx));
        x  = x + lamb * Dx;
        fx = fz;

    } while (true);

    return x;
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
    // ------------PART A--------------
    std::cout << "-----------PART A-------------" << std::endl;
    // analytical gradient of rosenbrock
    auto rosenbrock = [](pp::vector v) {
        double x = v[0];
        double y = v[1];
        double gx = -2 * (1 - x) + 200 * (y - x*x) * (-2 * x);
        double gy = 200 * (y - x*x);
        pp::vector result {gx, gy};
        return result;
        
    };
    // analytical gradient of himmelblau
    auto himmelblau = [](pp::vector v) {
        double x = v[0];
        double y = v[1];
        double gx = 4 * x * (x * x + y - 11) + 2 * (x + y * y - 7);
        double gy = 2 * (x * x + y - 11) + 4 * y * (x + y * y - 7);
        pp::vector result {gx, gy};
        return result;
    };

    // debugging newton method by using rosenbrock and himmelblau
    pp::vector start_rosen{-3.0, -4.0};
    pp::vector roots_rosen = newton(rosenbrock, start_rosen);
    std::cout << "df(x)=0 for rosenbrock at x= ";
    roots_rosen.print(" ");
    std::cout << "Analytical result = (1, 1)" << std::endl;

    pp::vector start_himmel{-0.5, -1.0};
    pp::vector roots_himmel = newton(himmelblau, start_himmel);
    std::cout << "f(x)=0 for himmelblau at x= ";
    roots_himmel.print(" ");
    std::cout << "Analytical result = (-0.270845, -0.923039)" << std::endl;

    // ---------PART B------------
    std::cout << "---------PART B------------" << std::endl;
    auto FE = [](pp::vector E) {
        double e = E[0];
        auto wave = [e](double r, pp::vector f) {
            pp::vector dfdr(2);
            dfdr[0] = f[1];
            dfdr[1] = -2.0 * (e * f[0] + 1/r * f[0]);
            return dfdr;
        };
        double rmin = 0.01;
        double rmax = 8.0;
        pp::vector f_init {rmin - rmin*rmin, 1.0 - 2.0*rmin};

        auto [rs, fs] = pp::driver(wave, {rmin, rmax}, f_init);

        std::vector r = linspace(0, 8, (int)rs.size());
        std::ofstream outfile("wave.dat");
        for (int j = 0; j < (int)rs.size(); j++) {
            double exact = r[j] * std::exp(-r[j]);
            outfile << rs[j] << " " << fs[j][0] << " " << exact << std::endl;
        }
        outfile.close();

        return pp::vector{fs[fs.size()-1][0]};
    };
    pp::vector start_wave {-0.9};
    pp::vector roots_wave = newton(FE, start_wave);

    roots_wave.print("E=");
    std::cout << "Exact result: E=-1/2" << std::endl;

    // ---------PART C------------
    std::cout << "---------PART C------------" << std::endl;
    pp::vector roots_wave_new = newton_interp(FE, start_wave);
    std::cout << "Results with Quadratic Interpolation line-search" << std::endl;
    roots_wave_new.print("E=");
    std::cout << "Exact result: E=-1/2" << std::endl;

    return 0;
}