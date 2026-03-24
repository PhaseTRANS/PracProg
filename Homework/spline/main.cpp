#include<stdexcept>
#include<vector>
#include<cmath>
#include<fstream>
#include"matrix.hpp"

pp::vector linspace(double start, double stop, int num) {
    pp::vector result(num);

    if (num == 1) {
        result[0] = start;
        return result;
    }

    double step = (stop - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

int binsearch(const pp::vector& x, double z) {
    /* locates the interval for z by bisection */
    if (z < x[0] || z > x[x.size() - 1]) {
        throw std::invalid_argument("invalid");
    }
    int i = 0, j = x.size() - 1;
    while (j - i > 1) {
        int mid = (i + j) / 2;
        if (z > x[mid]) i = mid; else j = mid;
    }
    return i;
}

struct lspline {
    pp::vector x, y;

    double linterp(const pp::vector& x, const pp::vector& y, double z) {
        int i = binsearch(x, z);
        double dx = x[i + 1] - x[i];
        if (!(dx > 0)) {
            throw std::invalid_argument("uups...");
        }
        double dy = y[i + 1] - y[i];
        return y[i] + dy / dx * (z - x[i]);
    }

    double linterpInteg(const pp::vector& x, const pp::vector& y, double z) {
        int i = binsearch(x, z);  // find which interval z falls in
        
        double sum = 0;
        
        // sum up the complete trapezoids before interval i
        for (int k = 0; k < i; k++) {
            double dx = x[k+1] - x[k];
            sum += dx * (y[k] + y[k+1]) / 2.0;  // trapezoid area
        }
        
        // add the partial trapezoid from x[i] to z
        double dx = z - x[i];
        double y_z = linterp(x, y, z);       
        sum += dx * (y[i] + y_z) / 2.0;         
        
        return sum;
    }

    // wrappers
    double evaluate(double z) {
        return linterp(x, y, z);
    }
    double integrate(double z) {
        return linterpInteg(x, y, z);
    }
};

struct qspline {
    pp::vector x, y, b, c;

    void qinterp(const pp::vector& xs, const pp::vector& ys) {
        int n = xs.size();
        x = xs;
        y = ys;
        b.resize(n - 1);
        c.resize(n - 1);

        c[0] = 0;
        for (int i = 0; i < n - 2; i++) {
            double dx  = x[i+1] - x[i];
            double dy  = y[i+1] - y[i];
            double dx1 = x[i+2] - x[i+1];
            double dy1 = y[i+2] - y[i+1];
            c[i+1] = (dy1/dx1 - dy/dx + c[i]*dx) / dx1;
        }

        pp::vector c2(n - 1);
        c2[n-2] = 0;
        for (int i = n - 3; i >= 0; i--) {
            double dx  = x[i+1] - x[i];
            double dy  = y[i+1] - y[i];
            double dx1 = x[i+2] - x[i+1];
            double dy1 = y[i+2] - y[i+1];
            c2[i] = (dy/dx - dy1/dx1 + c2[i+1]*dx1) / dx;
        }
        for (int i = 0; i < n - 1; i++) {
            c[i] = (c[i] + c2[i]) / 2.0; 
        }

        // compute b coefficients from c
        for (int i = 0; i < n - 1; i++) {
            double dx = x[i+1] - x[i];
            double dy = y[i+1] - y[i];
            b[i] = dy/dx - c[i]*dx;
        }
    }

    double evaluate(double z) {
        /* evaluate the spline: y[i] + b[i]*(z-x[i]) + c[i]*(z-x[i])^2 */
        qinterp(x, y);
        int i = binsearch(x, z);
        double dx = z - x[i];
        return y[i] + b[i]*dx + c[i]*dx*dx;
    }

    double derivative(double z) {
        /* derivative: b[i] + 2*c[i]*(z-x[i]) */
        int i = binsearch(x, z);
        double dx = z - x[i];
        return b[i] + 2*c[i]*dx;
    }

    double integrate(double z) {
        /* analytic integral from x[0] to z */
        int i = binsearch(x, z);
        double sum = 0;

        for (int k = 0; k < i; k++) {
            double dx = x[k+1] - x[k];
            sum += y[k]*dx + b[k]*dx*dx/2.0 + c[k]*dx*dx*dx/3.0;
        }

        double dx = z - x[i];
        sum += y[i]*dx + b[i]*dx*dx/2.0 + c[i]*dx*dx*dx/3.0;
        return sum;
    }
};

struct cspline {
    pp::vector x, y, b, c, d;

    void cinterp(const pp::vector& xs, const pp::vector& ys) {
        int n = xs.size();
        x = xs;
        y = ys;
        b.resize(n);
        c.resize(n - 1);
        d.resize(n - 1);

        pp::vector h(n - 1), p(n - 1);
        for (int i = 0; i < n - 1; i++) {
            h[i] = x[i+1] - x[i];
            p[i] = (y[i+1] - y[i]) / h[i];
        }

        // build tridiagonal system
        pp::vector D(n), Q(n - 1), B(n);
        D[0] = 2;
        for (int i = 0; i < n - 2; i++) {
            D[i+1] = 2*h[i]/h[i+1] + 2;
            D[n-1] = 2;
        }
        Q[0] = 1;
        for (int i = 0; i < n - 2; i++) {
            Q[i+1] = h[i]/h[i+1];
        }
        B[0] = 3*p[0];
        for (int i = 0; i < n - 2; i++) {
            B[i+1] = 3*(p[i] + p[i+1]*h[i]/h[i+1]);
            B[n-1] = 3*p[n-2];
        }

        // Gauss elimination
        for (int i = 1; i < n; i++) {
            D[i] -= Q[i-1] / D[i-1];
            B[i] -= B[i-1] / D[i-1];
        }

        // back-substitution
        b[n-1] = B[n-1] / D[n-1];
        for (int i = n - 2; i >= 0; i--) {
            b[i] = (B[i] - Q[i]*b[i+1]) / D[i];
        }

        // compute c and d coefficients
        for (int i = 0; i < n - 1; i++) {
            c[i] = (-2*b[i] - b[i+1] + 3*p[i]) / h[i];
            d[i] = ( b[i] + b[i+1] - 2*p[i]) / h[i] / h[i];
        }
    }

    double evaluate(double z) {
        /* evaluate the spline: y[i] + b[i]*(z-x[i]) + c[i]*(z-x[i])^2 + d[i]*(z-x[i])^3 */
        cinterp(x, y);
        int i = binsearch(x, z);
        double dx = z - x[i];
        return y[i] + dx*(b[i] + dx*(c[i] + dx*d[i]));
    }

    double derivative(double z) {
        /* derivative: b[i] + 2*c[i]*(z-x[i]) + 3*d[i]*(z-x[i])^2 */
        int i = binsearch(x, z);
        double dx = z - x[i];
        return b[i] + dx*(2*c[i] + dx*3*d[i]);
    }

    double integrate(double z) {
        /* analytic integral from x[0] to z */
        int i = binsearch(x, z);
        double sum = 0;

        for (int k = 0; k < i; k++) {
            double dx = x[k+1] - x[k];
            sum += dx*(y[k] + dx*(b[k]/2.0 + dx*(c[k]/3.0 + dx*d[k]/4.0)));
        }

        double dx = z - x[i];
        sum += dx*(y[i] + dx*(b[i]/2.0 + dx*(c[i]/3.0 + dx*d[i]/4.0)));
        return sum;
    }
};

int main() {
    // ---------------PART A-----------------
    // build the tabulated data
    pp::vector xs {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    pp::vector ys(xs.size());
    for (int i=0; i<xs.size(); i++) {
        ys[i] = std::cos(xs[i]);
    }

    lspline lspline;
    lspline.x = xs;
    lspline.y = ys;

    // save tabulated data
    std::ofstream table("tabulated.dat");
    for (int i=0; i<xs.size(); i++) {
        table << xs[i] << " " << ys[i] << std::endl;
    }
    table.close();

    // save interpolated data
    std::ofstream linterp("lspline.dat");
    pp::vector zs = linspace(xs[0], xs[xs.size() - 1], 100);
    for (int i=0; i<zs.size(); i++) {
        double z = lspline.evaluate(zs[i]);
        double s = lspline.integrate(zs[i]);

        linterp << zs[i] << " " << z << " " << s << std::endl;
    }
    linterp.close();

    //-----------PART B---------------
    qspline qspline;
    qspline.x = xs;
    qspline.y = ys;

    // save interpolated data
    std::ofstream qinterp("qspline.dat");
    for (int i=0; i<zs.size(); i++) {
        double z = qspline.evaluate(zs[i]);
        double s = qspline.integrate(zs[i]);
        double d = qspline.derivative(zs[i]);

        qinterp << zs[i] << " " << z << " " << s << " " << d << std::endl;
    }
    qinterp.close();

    //-----------PART C---------------
    cspline cspline;
    cspline.x = xs;
    cspline.y = ys;

    // save interpolated data
    std::ofstream cinterp("cspline.dat");
    for (int i=0; i<zs.size(); i++) {
        std::cout << zs[i] << std::endl;
        double z = cspline.evaluate(zs[i]);
        double s = cspline.integrate(zs[i]);
        double d = cspline.derivative(zs[i]);

        cinterp << zs[i] << " " << z << " " << s << " " << d << std::endl;
    }
    cinterp.close();

    return 0;
}

