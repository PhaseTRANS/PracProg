#include"matrix.hpp"
#include<iostream>

namespace pp{
    struct QR{
        std::pair<matrix, matrix> decomp(const matrix& A) {
        matrix Q = A.copy();
        matrix R(A.size2(), A.size2());
        // orthogonalize Q and fill-in R
        return {Q, R};
    }

    vector solve(const matrix& Q, const matrix& R, const vector& b) {
        // ...
    }

    double det(const matrix& R) {
        // ...
    }

    matrix inverse(const matrix& Q, const matrix& R) {
        // ...
    }
    };
}