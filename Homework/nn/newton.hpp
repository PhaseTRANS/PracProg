#include<cmath>
#include<functional>
#include<fstream>
#include"matrix.hpp"
#include"qr.hpp"

namespace pp {
    pp::vector gradient(
        std::function<double(pp::vector)> f,
        pp::vector x
    ) {
        double fx = f(x);
        pp::vector gf ((int)x.size());

        for (int i=0; i<(int)x.size(); i++) {
            double dxi = (1 + std::abs(x[i])) * std::pow(2, -26);
            x[i] += dxi;
            gf[i] = (f(x) - fx) / dxi;
            x[i] -= dxi;
        }
        return gf;
    }

    pp::matrix hessian(
        std::function<double(pp::vector)> f,
        pp::vector x
    ) {
        int x_size = x.size();
        pp::matrix H(x_size, x_size);
        pp::vector gfx = gradient(f, x);

        for (int j=0; j<x_size; j++) {
            double dxj = (1 + std::abs(x[j])) * std::pow(2, -13);
            x[j] += dxj;
            pp::vector dgf = gradient(f, x) - gfx;
            for (int i=0; i<x_size;i++) {
                H[i, j] = dgf[i] / dxj;
            }
            x[j] -=dxj;
        }
        return H;
    }

    pp::vector newton(
        std::function<double(pp::vector)> f,
        pp::vector x,
        double acc=1e-3,
        int maxiter=1000,
        std::function<pp::vector(std::function<double(pp::vector)>, pp::vector)> gradient=gradient,
        std::function<pp::matrix(std::function<double(pp::vector)>, pp::vector)> hessian=hessian
    ) {
        int counter = 0;
        while (counter < maxiter) {
            pp::vector g = gradient(f, x);
            if (g.norm() < acc) break;
            pp::matrix H = hessian(f, x);
            pp::QR QRH;
            auto [Q, R] = QRH.decomp(H);
            pp::vector dx = QRH.solve(Q, R, -g);
            double lambda = 1;
            while (lambda >= 1.0/1024) {
                if (f(x + lambda * dx) < f(x)) break;
                lambda /= 2;
            }
            x = x + lambda * dx;
            counter++;
        }
        return x;
    }
}