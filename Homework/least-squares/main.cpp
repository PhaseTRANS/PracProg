#include"matrix.hpp"
#include"qr.hpp"
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

// least squares fit function using QR-decomposition
std::pair<pp::vector, pp::matrix> lsfit(
    std::vector<std::function<double(double)>> fs, 
    pp::vector x, 
    pp::vector y, 
    pp::vector dy) {
        // build A
        int n = x.size();
        int m = fs.size();

        pp::matrix A(n, m);        
        for (int i = 0; i < m; i++) {
            std::function f = fs[i];
            for (int j = 0; j < n; j++) {
                A[j, i] = f(x[j]);  // row j, column i
            }
        }

        // use qr-decomposition
        pp::QR qr;
        auto [Q, R] = qr.decomp(A);
        pp::vector c = qr.solve(Q, R, y);

        // calculate covariance matrix
        pp::matrix I(R);
        I.setid();
        pp::matrix RI = qr.inverse(I, R);
        pp::matrix Sigma = RI * RI.transpose();

        return {c, Sigma};
}

int main() {
    // PART A
    std::cout << "----- PART A-----" << std::endl;
    // initialize data
    pp::vector x {1.0, 2.0, 3.0, 4.0, 6.0, 9.0, 10.0, 13.0, 15.0};
    pp::vector y {117.0, 100.0, 88.0, 72.0, 53.0, 29.5, 25.2, 15.2, 11.1};
    pp::vector dy {6, 5, 4, 4, 4, 3, 3, 2, 2};

    // using std::vector since pp::vector is not templated for function
    std::vector<std::function<double(double)>> fs = {
        [](double z) { return 1.0; },
        [](double z) { return -z; }
	};

    // solve least-squares
    pp::vector logy(y);
    pp::vector logdy(y);
    for (int i=0; i<y.size(); i++) {
        logy[i] = std::log(y[i]);
        logdy[i] = dy[i]/y[i];
    }

    auto [c, Sigma] = lsfit(fs, x, logy, logdy);

    // safe data to file
    std::ofstream data("decay.dat");
    for (int i=0; i<x.size(); i++) {
        data << x[i] << " " << " " << logy[i] << " " << logdy[i] << std::endl;
    }
    data.close();

    // safe fit to file
    std::ofstream fit("fit.dat");
    double x_min = 1.0, x_max = 15.0;
    int npoints = 200;
    for (int k = 0; k < npoints; k++) {
        double xi = x_min + k * (x_max - x_min) / (npoints - 1);
        double yi = 0;
        for (int i = 0; i < c.size(); i++)  // cs = your coefficients
            yi += c[i] * fs[i](xi);
        fit << xi << " " << yi << std::endl;
    }
    fit.close();

    // estimate half-life
    double lamb = c[1];
    double halflife = std::log(2) / lamb;
    std::cout << "Theoretical Value: (3.6316 +/- 0.0014)days" << std::endl;
    std::cout  << "T_½ = " << halflife << "days" << std::endl;

    // PART B
    std::cout << "----- PART B-----" << std::endl;
    Sigma.print("Sigma = ");
    double dhalflife = std::sqrt(Sigma[1][1]);
    std::cout  << "T_½ = (" << halflife << " +/- " << dhalflife << ")days" << std::endl;

    // PART C
    std::cout << "----- PART C-----" << std::endl;

    // Read existing file and essentially delete its content
    std::ifstream fit_in("fit.dat");
    std::vector<std::pair<double,double>> fit_data;
    double xi, yi;
    while (fit_in >> xi >> yi)
        fit_data.push_back({xi, yi});
    fit_in.close();

    // Rewrite file with the extra columns
    std::ofstream fit_out("fit.dat");
    for (auto& [xi, yi] : fit_data) {
        double yi_up = 0, yi_down = 0;
        for (int i = 0; i < c.size(); i++) {
            yi_up   += (c[i] + std::sqrt(Sigma[i][i])) * fs[i](xi);
            yi_down += (c[i] - std::sqrt(Sigma[i][i])) * fs[i](xi);
        }
        fit_out << xi << " " << yi << " " << yi_up << " " << yi_down << "\n";
    }
    fit_out.close();
    return 0;
}