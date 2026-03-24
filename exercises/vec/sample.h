#pragma once
#include<string>
#include<iostream>

namespace pp{
    struct vec{
        double x, y, z;

        // ctors
        vec(double a, double b, double c) {  // parametrized ctor
            std::cout << "parameterized constructer called..." << std::endl;
            x=a; y=b; z=c;
        }

        vec() : vec(0, 0, 0) {
            std::cout << "default constructer called..." << std::endl;
        }
        vec(const vec&) = default; // copy ctor; vec a=b;
        vec(vec&&) = default;  // move ctor: vec a=b+c

        // dtor
        ~vec() {
            std::cout << "destructer called..." << std::endl;
        }

        // assignments
        vec& operator=(const vec&);  // copy assignment
        vec& operator=(vec&&);   // move assignment


        // member operators
        vec& operator+=(double);
        vec& operator-=(double);
        vec& operator*=(double);
        vec& operator/=(double);

        void print(const std::string& s="");

        // stream output
        friend std::ostream& operator<<(std::ostream&, const vec&);
    };
    // non-members
    vec operator+(const vec&, const vec&);
    vec operator-(const vec&, const vec&);
    vec operator*(const vec&, double);
    vec operator*(double, const vec&);
    vec operator/(const vec&, double);
    vec operator-(const vec&);
}