#pragma once
#include<iostream>
#include<cstdio>
#include<iomanip>
#include<cmath>
#include<string>
#include<vector>
#include<initializer_list>
#include<functional>
#include<ranges>
#include<cmath>

namespace pp{

struct vector {
	std::vector<double> data;

	vector() = default;
	vector(int n) : data(n) {}
	vector(std::initializer_list<double> list) : data(list) {}
	vector(const vector&) = default;
	vector(vector&&) noexcept = default;

	vector& operator=(const vector&) = default;
	vector& operator=(vector&&) noexcept = default;

	inline int size() const {return data.size();}
	//auto n(){return std::views::iota(0,size());}
	void resize(int n) {data.resize(n);}
	inline double& operator[](int i) {return data[i];}
	inline const double& operator[](int i) const {return data[i];}

	vector& operator+=(const vector& other){
		for(int i=0;i<size();i++)(*this)[i]+=other[i];
		//for(int i:n())(*this)[i]+=other[i];
		return (*this);
		}

	vector& operator-=(const vector& other){
		for(int i=0;i<size();i++)(*this)[i]-=other[i];
		return (*this);
		}

	vector& operator*=(double c){
		for(int i=0;i<size();i++)(*this)[i]*=c;
		return (*this);
		}

	vector& operator/=(double c){
		for(int i=0;i<size();i++)(*this)[i]/=c;
		return (*this);
		}

	double norm() const {
		double sum2=0;
		for(int i=0;i<size();i++)sum2+=(*this)[i]*(*this)[i];
		return std::sqrt(sum2);
	}

	double dot(const vector& other) const {
		double sum1 = 0;
		for (int i = 0; i < size(); i++) {
			sum1 += (*this)[i] * other[i];
		}
		return sum1;
	}

	void print(std::string s="") const {
		std::cout<<s<<" ";
		for(auto &x : data) std::cout<<x<<" ";
		std::cout<<"\n";
	}

	vector map(std::function<double(double)> f) const{
		vector r(size());
		for(int i=0;i<size();i++)r.data[i]=f(data[i]);
		return r;
	}

}; //vector

vector operator+(vector a, const vector& b)	{ a+=b ; return a; }
vector operator-(vector a)					{ a*=-1; return a; }
vector operator-(vector a, const vector& b)	{ a-=b ; return a; }
vector operator*(vector a, const double c)	{ a*=c ; return a; }
vector operator*(const double c, vector a)	{ a*=c ; return a; }
vector operator/(vector a, const double c)	{ a/=c ; return a; }

bool approx(double x, double y, double acc=1e-6, double eps=1e-6){
	if(std::abs(x-y) < acc) return true;
	if(std::abs(x-y) < eps*std::max(std::abs(x),std::abs(y))) return true;
	return false;
	}

bool approx(const vector& a, const vector& b, double acc=1e-6, double eps=1e-6){
	if(a.size() != b.size()) return false;
	for(int i=0;i<a.size();i++)
		if(!approx(a[i],b[i],acc,eps)) return false;
	return true;
	}

struct matrix {
	std::vector<pp::vector> cols;
	matrix()=default;
	matrix(int n,int m) : cols(m, pp::vector(n)) {}
	matrix(const matrix& other)=default;
	matrix(matrix&& other)=default;
	matrix& operator=(const matrix& other)=default;
	matrix& operator=(matrix&& other)=default;
	int size1() const {return cols.empty() ? 0 : cols[0].size(); }
	int size2() const {return cols.size();}
	double& operator()(int i, int j){return cols[j][i];}
	double& operator[](int i, int j){return cols[j][i];}
	const double& operator()(int i, int j)const{return cols[j][i];}
	const double& operator[](int i, int j)const{return cols[j][i];}
	vector& operator[](int i){return cols[i];}
	const vector& operator[](int i) const {return cols[i];}
//	void resize(int n, int m);
	void setid(){
		if(size1()!=size2())throw std::runtime_error("non-square matrix\n");
		for(int i=0;i<size1();i++){
			(*this)[i,i]=1;
			for(int j=i+1;j<size1();j++)(*this)[i,j]=(*this)[j,i]=0;
			}
		}
	matrix transpose() const{
		matrix R(size2(),size1());
		for(int i=0;i<size1();i++)
		for(int j=0;j<size2();j++)
			R[j,i]=(*this)[i,j];
		return R;
	}

//	matrix T() const;
	
	double get (int i, int j) {return cols[j][i];}
	void set(int i, int j, double value){cols[j][i] = value;}
//	vector get_col(int j);
//	void set_col(int j,vector& cj);

	matrix& operator+=(const matrix& B){
		for(int i=0;i<size2();i++)(*this)[i]+=B[i];
		return *this;
		}
	matrix& operator-=(const matrix& B){
		for(int i=0;i<size2();i++)(*this)[i]-=B[i];
		return *this;
		}
	matrix& operator*=(const double c){
		for(int i=0;i<size2();i++)(*this)[i]*=c;
		return *this;
		}
	matrix& operator/=(const double c){
		for(int i=0;i<size2();i++)(*this)[i]/=c;
		return *this;
		}
	matrix& operator*=(const matrix&);
	matrix  operator^(int);
	
	void print(std::string s="") const {
    std::cout << s << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    for(int i = 0; i < size1(); i++) {
        for(int j = 0; j < size2(); j++) {
            std::cout << std::setw(10) << (*this)(i, j);
        }
        std::cout << "\n";
    }
}
    matrix copy() const {
        return *this;
    }
};

matrix operator+(matrix A, const matrix& B){
	for(int i=0;i<A.size2();i++)A[i]+=B[i];
	return A;
	}

matrix operator-(matrix A, const matrix& B){
	for(int i=0;i<A.size2();i++)A[i]-=B[i];
	return A;
	}

matrix operator*(const matrix& A, const matrix& B){
	if(A.size2()!=B.size1()) throw std::invalid_argument("size mismatch");
	matrix R(A.size1(),B.size2());
	for(int k=0;k<A.size2();k++)
	for(int j=0;j<B.size2();j++) {
		double Bkj=B[k,j];
		for(int i=0;i<A.size1();i++)R[i,j]+=A[i,k]*Bkj;
		}
	return R;
	}

matrix operator*(matrix A, const double c){
	for(auto &col : A.cols) col*=c;
	return A;
	}
matrix operator*(const double c, matrix A){
	for(auto &col : A.cols) col*=c;
	return A;
	}
matrix operator/(matrix A, const double c){
	for(auto &col : A.cols) col/=c;
	return A;
	}
vector operator*(const matrix& A, const vector& v){
	vector r(A.size1());
	for(int j=0;j<A.size2();j++){
		double vj=v[j];
		for(int i=0;i<A.size1();i++) r[i]+=A[i,j]*vj;
		}
	return r;
	}
}