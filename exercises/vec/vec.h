#pragma once
#include<iostream>
#include<cmath>

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

struct vec{
    double x, y, z;

    // constructors
    vec(double x, double y, double z) : x(x), y(y), z(z) {}     // parameterized ctor
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
    vec& operator*=(double);
    vec& operator/=(double);

    // utility
    void set(double a, double b, double c) {
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
        // Fixed logic: return false if ANY component is NOT approximately equal
        if (!::approx(a.x, b.x)) {  // using :: to call global approx function
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

    // dot-product
    double dot(const vec& b) {
        double sum = x * b.x + y * b.y + z * b.z;
        return sum;
    }

    // cross-product
    vec cross(const vec& b) {
        double new_x = y * b.z - z * b.y;
        double new_y = z * b.x - x * b.z;
        double new_z = x * b.y - y * b.x;
        return vec(new_x, new_y, new_z);
    }

    // norm
    double norm() {
        double length = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
        return length;
    }
};

// non-member operators
// vec operator-(const vec&);
vec operator-(const vec& a, const vec& b) {
    return vec(a.x - b.x, a.y - b.y, a.z - b.z);
}
vec operator+(const vec& a, const vec& b) {
    return vec(a.x + b.x, a.y + b.y, a.z + b.z);
}
vec operator*(const vec& a, double c) {
    return vec(a.x * c, a.y * c, a.z * c);
}
vec operator*(double c, const vec& a) {
    return vec(c * a.x, c * a.y, c * a.z);
}
vec operator/(const vec& a, double c) {
    return vec(a.x / c, a.y / c, a.z / c);
}
std::ostream& operator<<(std::ostream& os, const vec& v){
    os << "{ " << v.x << ", " << v.y << ", " << v.z << " } ";
    return os;
}