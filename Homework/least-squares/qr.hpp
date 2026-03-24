#pragma once
#include"matrix.hpp" 

namespace pp{
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
        for (int i=0; i<R.size1(); i++) {
            prod *= R[i][i];
        }
        return prod;

    }

    pp::matrix inverse(const pp::matrix& Q, const pp::matrix& R) {
        // since A=RQ we have A^-1 = R^-1 * Q^T, so we need to find R^-1, which is upper triangular
        // thus we can find R^-1 column by column by solving R*c_k = e_k, where e_k is the k'th
        // basis vector which has 1 at the k'th position and 0 else
        pp::matrix R_inv(R);
        for (int i=0; i<R.size2(); i++) {
            pp::vector ek(R.size1());
            ek[i] = 1;
            pp::vector ck = back_substitution(R, ek);
            R_inv[i] = ck;
        }
        pp::matrix A_inv = R_inv * Q.transpose();
        return A_inv;
    }

    
};
}