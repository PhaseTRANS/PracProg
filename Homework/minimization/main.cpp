#include<iostream>
#include<vector>
#include<cmath>
#include<functional>
#include<fstream>
#include"matrix.hpp"
#include"qr.hpp"

// Python implementations
// def newton(φ,x,acc=1e-3) :
//     while True :                    # Newton iterations
//         g = gradient(φ,x)
//         if g.norm() < acc : break   # job done
//         H = hessian(φ,x)
//         dx = QRdecomposition(H).solve(-g)
//         λ = 1 ;
//         while λ ≥ 1/1024 :       # backtracking linesearch
//             if φ(x+λ*dx) < φ(x) : break # good step
//             λ /= 2
//         x=x+λ*dx
//     return x

// def gradient(φ,x) :
//     φx = φ(x)
//     gφ = vector(len(x))
//     for i in range(len(x)) :
//         dxi = (1+abs(x[i]))*2**(-26)
//         x[i]+=dxi
//         gφ[i]=(φ(x)-φx)/dxi
//         x[i]-=dxi
//     return gφ

// def hessian(φ,x) :
//     H = matrix(len(x),len(x))
//     gφx = gradient(φ,x)
//     for j in range(len(x)) :
//         dxj=(1+abs(x[j]))*2**(-13)
//         x[j]+=dxj
//         dgφ=gradient(φ,x)-gφx
//         for i in range(len(x)) : H[i,j]=dgφ[i]/dxj
//         x[j]-=dxj
//     return H

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

pp::vector gradient_central(
    std::function<double(pp::vector)> f,
    pp::vector x
) {
    pp::vector gf ((int)x.size());

    for (int i=0; i<(int)x.size(); i++) {
        double dxi = (1 + std::abs(x[i])) * std::pow(2, -26);
        x[i] += dxi;
        double dx_plus = f(x);
        x[i] -= dxi;
        x[i] -= dxi;
        double dx_minus = f(x);
        x[i] += dxi;
        gf[i] = (dx_plus - dx_minus) / (2 * dxi);
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

pp::matrix hessian_central(
    std::function<double(pp::vector)> f,
    pp::vector x
) {
    int n = x.size();
    pp::matrix H(n, n);
    pp::vector gfx = gradient_central(f, x);

    for (int j = 0; j < n; j++) {
        double dxj = (1 + std::abs(x[j])) * std::pow(2, -13);
        x[j] += dxj;
        pp::vector g_plus = gradient_central(f, x);
        x[j] -= dxj;
        x[j] -= dxj;
        pp::vector g_minus = gradient_central(f, x);
        x[j] += dxj;

        for (int i = 0; i < n; i++) {
            H[i, j] = (g_plus[i] - g_minus[i]) / (2 * dxj);
        }
    }
    return H;
}

pp::vector newton(
    std::function<double(pp::vector)> f,
    std::function<pp::vector(std::function<double(pp::vector)>, pp::vector)> gradient,
    std::function<pp::matrix(std::function<double(pp::vector)>, pp::vector)> hessian,
    pp::vector x,
    double acc=1e-3,
    int maxiter=1000
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

int main() {
    // PART A
    std::cout << "-----------PART A-------------" << std::endl;
    int counter_rosen = 0;
    auto rosenbrock = [&counter_rosen](pp::vector v) {
        double x = v[0];
        double y = v[1];
        double f = (1 - x) * (1 - x) + 100 * (y - x * x) * (y - x * x);
        counter_rosen++;
        return f;
    };
    int counter_himmel = 0;
    auto himmelblau = [&counter_himmel](pp::vector v) {
        double x = v[0];
        double y = v[1];
        double f = (x * x + y - 11) * (x * x + y - 11) + (x + y * y - 7) * (x + y * y - 7);
        counter_himmel++;
        return f;
    };

    // find minimum of rosenbrock
    pp::vector rosen_start {0.5, 1.5};
    pp::vector min_rosen = newton(rosenbrock, gradient, hessian, rosen_start);

    // find minimum of himmelblau
    pp::vector himmel_start {6.3, 6.4};
    pp::vector min_himmel = newton(himmelblau, gradient, hessian, himmel_start);

    min_rosen.print("Rosenbrock minimum at: ");
    std::cout << "Found with " << counter_rosen << " evaluations. Analytic minimum rosenbrock: (1.0, 1.0)" << std::endl;
    min_himmel.print("Himmelblau minimum at: ");
    std::cout << "Found with " << counter_himmel << " evaluations. Analytic minimum Himmelblau: (3.0, 2.0)" << std::endl;

    // PART B
    std::cout << "----------PART B--------------" << std::endl;
    auto breit_wigner = [](double e, pp::vector p) {
        double A = p[0];
        double m = p[1];
        double gamma = p[2];
        double F = A / ((e - m) * (e - m) + gamma * gamma / 4);
        return F;
    };

    // read data from higgs.dat
    std::ifstream file("higgs.dat");
    pp::vector E, S, S_err;
    double e, s, s_err;
    while (file >> e >> s >> s_err) {
        E.data.push_back(e);
        S.data.push_back(s);
        S_err.data.push_back(s_err);
    }
    file.close();

    // check if the data is loaded correctly
    for (int i=0; i<E.size(); i++) {
        std::cout << E[i] << " " << S[i] << " " << S_err[i] << std::endl;
    }
    std::cout << "Loaded " << E.size() << " data points" << std::endl;

    auto resid_error = [breit_wigner, E, S, S_err](pp::vector p) {
        double error = 0;
        for (int i=0; i<E.size(); i++) {
            double F = breit_wigner(E[i], p);
            error += ((F - S[i]) / S_err[i]) * ((F - S[i]) / S_err[i]);
        }
        return error;
    };

    // make fit
    pp::vector p_guess {1.0, 120.0, 1.0};
    pp::vector p_opt = newton(resid_error, gradient, hessian, p_guess, 0.001, 1000);
    p_opt.print("(A, m, Gamma) = ");

    // safe to file
    std::vector Es = linspace(100, 170, 200);
    std::ofstream outfile("higgs_fit.dat");
    for (int i=0; i<(int)Es.size(); i++) {
        outfile << Es[i] << " " << S[i] << " " << S_err[i] << " " << breit_wigner(Es[i], p_opt) << std::endl;
    }
    outfile.close();

    // PART C
    std::cout << "----------PART C-------------" << std::endl;
    // do the same as in part a with central difference
    pp::vector min_rosen_central = newton(rosenbrock, gradient_central, hessian_central, rosen_start);
    pp::vector min_himmel_central = newton(himmelblau, gradient_central, hessian_central, himmel_start);

    min_rosen_central.print("Rosenbrock minimum at: ");
    min_himmel_central.print("Himmelblau minimum at: ");

    // evaluating which of the methods varies the less from the exact result
    pp::vector exact_rosen {1.0, 1.0};
    pp::vector exact_himmel {3.0, 2.0};
    pp::vector difference_rosen_foward = min_rosen - exact_rosen;
    pp::vector difference_himmel_foward = min_himmel - exact_himmel;
    difference_rosen_foward.print("Rosenbrock: Exact - forward = ");
    difference_himmel_foward.print("Himmelblau: Exact - forward = ");
    pp::vector difference_rosen_central= min_rosen_central - exact_rosen;
    pp::vector difference_himmel_central = min_himmel_central - exact_himmel;
    difference_rosen_central.print("Rosenbrock: Exact - central = ");
    difference_himmel_central.print("Himmelblau: Exact - central = ");

    return 0;
}