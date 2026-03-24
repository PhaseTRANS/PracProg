#include"matrix.hpp" 
#include<vector> 
#include<iostream> 
#include<random>


pp::matrix random_matrix(int n, int m) {
    if (n < m) {
        throw std::invalid_argument("Require n > m");
    }

    // Random device for seeding
    std::random_device rd;

    // Mersenne Twister engine
    std::mt19937 gen(rd());

    // Uniform distribution between 0 and 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    pp::matrix A(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A.cols[j][i] = dist(gen);
        }
    }

    return A;
}

pp::vector random_vector(int n) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    pp::vector b(n);

    for (int j=0; j<n; j++) {
        b[j] = dist(gen);
    }
    return b;
}

struct QR{
    std::pair<pp::matrix, pp::matrix>  decomp(const pp::matrix& A) {
        pp::matrix Q = A.copy();
        pp::matrix R(A.size2(), A.size2());
        // orthogonalize Q and fill-in R
        int m = A.size2();
        for (int i = 0; i < m; ++i) {
            pp::vector v = A[i];
            for (int j = 0; j < i; ++j) {
                double r = v.dot(Q[j]);  // project updated v, not original A[i]
                R(j, i) = r;
                v = v - r * Q[j];
            }

            double norm_v = v.norm();
            R[i][i] = norm_v;
            Q[i] = v / norm_v;
        }
        return {Q, R};
    }

    pp::vector back_substitution(const pp::matrix& R, const pp::vector& y) {
        // helper function for back substitution, since it is needed seperatly and in solve
        pp::vector x(R.size2());
        for (int i = R.size2() - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < R.size2(); j++) {
                sum += R(i, j) * x[j];
            }
            x[i] = (y[i] - sum) / R(i, i);
        }
        return x;
    }   

    pp::vector solve(const pp::matrix& Q, const pp::matrix& R, const pp::vector& b) {
        pp::vector y(Q.size2());
        for (int i = 0; i < Q.size2(); i++) {
            y[i] = Q[i].dot(b);
        }
        return back_substitution(R, y);
    }

    double det(const pp::matrix& R) {
        // the determinant of a upper triangular matrix is the product 
        // of the diagonal elements
        double prod = 1;
        for (int i; i<R.size1(); i++) {
            prod *= R[i][i];
        }
        return prod;

    }

    pp::matrix inverse(const pp::matrix& Q, const pp::matrix& R) {
        // since A=RQ we have A^-1 = R^-1 * Q^T, so we need to find R^-1, which is upper triangular
        // thus we can find R^-1 column by column by solving R*c_k = e_k, where e_k is the k'th
        // basis vector which has 1 at the k'th position and 0 else
        pp::matrix R_inv(R);
        for (int i; i<R.size2(); i++) {
            pp::vector ek(R.size1());
            ek[i] = 1;
            pp::vector ck = back_substitution(R, ek);
            R_inv[i] = ck;
        }
        pp::matrix A_inv = R_inv * Q.transpose();
        return A_inv;
    }

    
};

int main(int argc, char** argv) {
    // PART 1
    std::cout << "----------PART A----------" << std::endl;
    int n = std::stoi(argv[1]);
    int m = n;

    // generate random matrix
    pp::matrix A = random_matrix(n, m);

    A.print("A=");
    
    // initialize QR
    QR qr;
    auto [Q, R] = qr.decomp(A);

    // check if it works
    std::cout << "------Checking decomp function---------" << std::endl;
    Q.print("Q=");

    // check if R is upper triangular
    R.print("R=");

    // check if QR = A
    pp::matrix QR = Q * R;
    QR.print("Q * R = A = ");

    // check if Q^T Q = I
    pp::matrix eye = Q.transpose() * Q;
    eye.print("Q^T * Q=");

    // generate random vector
    // check if solve works
    // !!!!!!!!!!!!!solve is not working A*x is not equal to b!!!!!!!!!!!!!!!!
    std::cout << "-----Checking solve function-----" << std::endl;
    pp::vector b = random_vector(A.size1());
    b.print("b= ");
    pp::vector x = qr.solve(Q, R, b);
    x.print("x= ");
    pp::vector Ax = A * x;
    Ax.print("A * x = ");
    std::cout << "approx(A*x, b)= " << pp::approx(b, Ax) << std::endl;
    double determinant = qr.det(R);
    std::cout << "det(R)= " << determinant << std::endl;


    // PART 2
    std::cout << "----------PART B----------" << std::endl;

    // check if inverse works
    std::cout << "-----Checking inverse function-----" << std::endl;
    pp::matrix A_inv = qr.inverse(Q, R);
    A_inv.print("A^-1= ");
    pp::matrix eye2 = A * A_inv;
    eye2.print("A * A^-1=");

    return 0;
}