#include<iostream>
#include<functional>
#include<cmath>
#include<random>
#include<fstream>
#include"matrix.hpp"

// helpers
double random_double() {
    // Random device for seeding
    std::random_device rd;

    // Mersenne Twister engine
    std::mt19937 gen(rd());

    // Uniform distribution between 0 and 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
};
pp::vector linspace(double start, double end, int n) {
    pp::vector result(n);
    if (n <= 1) {
        result[0] = start;
        return result;
    }
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++) {
        result[i] = start + i * step;
    }
    return result;
};

// MC integrator
std::pair<double, double> plainmc(
    std::function<double(pp::vector)> f, 
    pp::vector a, 
    pp::vector b, 
    int N) {
    int dim = a.size(); double V = 1; 
    for (int i=0; i<dim; i++) {
        V*=b[i]-a[i];
    }
    double sum = 0, sum2 = 0;
    pp::vector x(dim);
    for(int i=0; i<N; i++){
        for(int k=0; k<dim; k++) {
            x[k] = a[k] + random_double() * (b[k] - a[k]);
        }
    double fx = f(x); sum += fx; sum2 += fx * fx;
    }
    double mean=sum/N, sigma=std::sqrt(sum2/N-mean*mean);
    return {mean * V, sigma * V / std::sqrt(N)};
};

// helper to create van der corput numbers

double van_der_corput(int n, int base) {
    double x = 0, denom = 1;
    while (n > 0) {
        denom *= base;
        x += (n % base) / denom;
        n /= base;
    }
    return x;
};

// helper to generate vector of primes
pp::vector prime_numbers(int n) {
    pp::vector primes(n);
    double candidate = 2;
    int found = 0;                          
    while (found < n) {
        bool is_prime = true;
        for (int i = 0; i < found; i++) {   
            if (primes[i] * primes[i] > candidate) break;
            if ((int)candidate % (int)primes[i] == 0) {
                is_prime = false;
                break;                     
            }
        }
        if (is_prime) {
            primes[found] = candidate;
            found++;                       
        }
        candidate++;
    }
    return primes;
}

std::pair<double, double> quasimc(
    std::function<double(pp::vector)> f,
    pp::vector a,
    pp::vector b,
    int N) {
    int dim = a.size();
    double V = 1;
    for (int i = 0; i < dim; i++) V *= b[i] - a[i];
    pp::vector primes = prime_numbers(dim);
    double sum = 0, sum2 = 0;
    pp::vector x(dim);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < dim; k++) {
            x[k] = a[k] + van_der_corput(i + 1, primes[k]) * (b[k] - a[k]);
        }
        double fx = f(x); sum += fx; sum2 += fx * fx;
    }
    double mean = sum / N;
    double sigma = std::sqrt(sum2 / N - mean * mean);
    return {mean * V, sigma * V / std::sqrt(N)};
};

// recursive stratified mc integration copied from lecture note and translated into C++
std::pair<double, double> strata(
      std::function<double(pp::vector)> f,
      pp::vector a,
      pp::vector b,
      int N,
      double acc = 0.01,
      double eps = 0.01) {

      int dim = a.size();
      auto [area, err] = plainmc(f, a, b, N);

      double tol = acc + std::abs(area) * eps;
      if (err < tol) return {area, err};  // also return the error!

      // find dimension with largest variance
      int idiv = 0;
      double maxvar = err;  // use error as proxy for variance

      // split in half along dimension idiv
      double midpoint = (a[idiv] + b[idiv]) / 2;

      // Left half: original a, new b (midpoint)
      pp::vector b_left = b;
      b_left[idiv] = midpoint;
      auto [area_left, err_left] = strata(f, a, b_left, N, acc, eps);

      // Right half: new a (midpoint), original b
      pp::vector a_right = a;
      a_right[idiv] = midpoint;
      auto [area_right, err_right] = strata(f, a_right, b, N, acc, eps);

      return {area_left + area_right, std::sqrt(err_left*err_left + err_right*err_right)};
  };

int main() {
    std::cout << "----------PART A------------" << std::endl;
    // define function for unit circle
    auto circ1 = [](pp::vector v) {
        double dist_sq = 0;
        for (int i = 0; i < v.size(); i++) {
            dist_sq += v[i] * v[i];
        }
        return (dist_sq <= 1.0) ? 1.0 : 0.0;
    };

    // make MC simulation to calculate the area of unit circle for debugging
    int N = 100000;
    pp::vector a = {-2, -2};
    pp::vector b = {2, 2};
    auto [circ_area, err] = plainmc(circ1, a, b, N);
    std::cout << "Area of unit circle: " << circ_area << " ≈ π" << std::endl; 

    // Calculate area of unit circle with different N
    std::ofstream circle_data("circle.dat");
    pp::vector Ns = linspace(100, 10000, 100);
    for (int i=0; i<Ns.size(); i++) {
        auto [area, err] = plainmc(circ1, a, b, Ns[i]);
        circle_data << Ns[i] << " " << err << " " << area << std::endl;
    }
    circle_data.close();

    // Calculate the complicated given integral
    auto gamma = [](pp::vector v) {
        double cosines = 1;
        for (int i=0; i<v.size(); i++) {
            cosines *= std::cos(v[i]);
        }
        double numerator = (1 - cosines) * std::pow(M_PI, 3);
        double value = 1/numerator;
        return value;
    };
    pp::vector uppers = {M_PI, M_PI, M_PI};
    pp::vector lowers = {0, 0, 0};
    auto [int_gamma, err_gamma] = plainmc(gamma, lowers, uppers, N);
    std::cout << "Volume of complicated function: " << int_gamma << " ≈ 1.393" << std::endl;

    std::cout << "-----------PART B------------" << std::endl;
    // running the new function on unit circle for comparison
    std::ofstream quasi_data("quasi.dat");
    for (int i=0; i<Ns.size(); i++) {
        auto [area, err] = quasimc(circ1, a, b, Ns[i]);
        quasi_data << Ns[i] << " " << err << " " << area << std::endl;
    }
    quasi_data.close();

    std::cout << "------------PART C------------" << std::endl;
    std::ofstream strata_data("strata.dat");
    for (int i=0; i<Ns.size(); i++) {
        auto [area, err] = strata(circ1, a, b, Ns[i]);
        strata_data << Ns[i] << " " << err << " " << area << std::endl;
    }
    strata_data.close();


    return 0;
}