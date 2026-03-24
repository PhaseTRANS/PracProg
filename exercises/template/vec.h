#pragma once
#include<iostream>
#include<cmath>
#include<complex>

bool approx(double a, double b, double acc=1e-9, double eps=1e-9) {
    double absolute = std::abs(a - b);
    double relative = absolute / std::max(std::abs(a), std::abs(b));
    if (absolute <= acc || relative <= eps) {
        return true;
    }
    else {
        return false;
    }
}

template <typename T>
struct vec{
    T x, y, z;

    // constructors
    vec(T x, T y, T z) : x(x), y(y), z(z) {}     // parameterized ctor
    vec() : vec(0, 0, 0) {}                                       // unparameterized ctor
    vec(const vec&) = default;                                  // copy
    vec(vec&&) = default;                                       // move
    ~vec() = default;                                           // destructor
    
    // assignments
    vec& operator=(const vec&) = default;                       // copy assignment
    vec& operator=(vec&&) = default;                            // move assignment

    // arithmetic
    vec& operator+=(const vec&);
    vec& operator-=(const vec&);
    vec& operator*=(const T& s);
    vec& operator/=(const T& s);

    // utility
    void set(T a, T b, T c) {
        x = a;
        y = b;
        z = c;
    }

    // debugging
    void print(const std::string& s = "") const {
        std::cout << s << x << " " << y << " " << z << std::endl;
    }

    bool approx(const double a, const double b) const;
    bool approx(const vec& a, const vec& b) const {
        if (!::approx(a.x, b.x)) {  
            return false;
        }
        if (!::approx(a.y, b.y)) {
            return false;
        }
        if (!::approx(a.z, b.z)) {
            return false;
        }
        return true;
    }

    // dot-product for complex and real
    auto dot(const vec& b) const {
        return x*b.x + y*b.y + z*b.z;
    }

    // cross-product for complex and real
    vec cross(const vec& b) const {
        T new_x = y*b.z - z*b.y;
        T new_y = z*b.x - x*b.z;
        T new_z = x*b.y - y*b.x;
        return vec(new_x, new_y, new_z);
    }


    // norm for complex and real
    auto norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

// non-member operators
// vec operator-(const vec&);
template <typename T>
vec<T> operator-(const vec<T>& a, const vec<T>& b) {
    return vec(a.x - b.x, a.y - b.y, a.z - b.z);
}
template <typename T>
vec<T> operator+(const vec<T>& a, const vec<T>& b) {
    return vec(a.x + b.x, a.y + b.y, a.z + b.z);
}
template <typename T>
vec<T> operator*(const vec<T>& a, const T& c) {
    return vec(a.x * c, a.y * c, a.z * c);
}
template <typename T>
vec<T> operator*(const T& c, const vec<T>& a) {
    return vec(c * a.x, c * a.y, c * a.z);
}
template <typename T>
vec<T> operator/(const vec<T>& a, const T& c) {
    return vec(a.x / c, a.y / c, a.z / c);
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const vec<T>& v){
    os << "{ " << v.x << ", " << v.y << ", " << v.z << " } ";
    return os;
}