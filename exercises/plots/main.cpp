#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<iomanip>
#include<vector>
#include<numbers>

double erf(double x){
// single precision error function (Abramowitz and Stegun, from Wikipedia)
    if(x<0) return -erf(-x);
    std::vector<double> a {0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
    double t=1/(1+0.3275911*x);
    double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
    return 1-sum*std::exp(-x*x);
} 

double sgamma(double x){
    constexpr double PI = std::numbers::pi; // c++20
    if(x<0)return PI/std::sin(PI*x)/sgamma(1-x);
    if(x<9)return sgamma(x+1)/x;
    double lnsgamma=std::log(2*PI)/2+(x-0.5)*std::log(x)-x
        +(1.0/12)/x-(1.0/360)/(x*x*x)+(1.0/1260)/(x*x*x*x*x);
    return std::exp(lnsgamma);
}

double lngamma(double x){
    constexpr double PI = std::numbers::pi; // c++20
    if(x<=0) return NAN;
    if(x<9) return lngamma(x+1)-std::log(x);
    return x*std::log(x+1/(12*x-1/x/10))-x+std::log(2*PI/x)/2;
}


int main(int argc,char** argv){
	double xmin=0,xmax=10,dx=0.125,dx2=0.01;
	for(int i=0;i<argc;i++){
		std::string arg=argv[i];
		std::cerr<<"i= "<<i<<" arg="<<arg<<"\n";
		if(arg=="-xmin" && i+1<argc)xmin=std::stod(argv[i+1]);
		if(arg=="-xmax" && i+1<argc)xmax=std::stod(argv[i+1]);
		if(arg=="-dx" && i+1<argc)dx=std::stod(argv[i+1]);
		if(arg=="-dx2" && i+1<argc)dx2=std::stod(argv[i+1]);
	}
	std::cerr<<"xmin= "<<xmin<<"\n";
	std::cerr<<"xmax= "<<xmax<<"\n";
	std::cerr<<"dx= "<<dx<<"\n";
	std::cerr<<"dx2= "<<dx2<<"\n";
	
    // save data for error function
	std::ofstream erf_file("erf.dat");
	erf_file<<std::scientific;
	for(double x=xmin;x<=xmax;x+=dx){
		erf_file<<x<<" "<<erf(x)<<" "<< std::erf(x) <<"\n";
	}
	erf_file<<"\n\n";
	
	for(double x=xmin;x<=xmax;x+=dx2){
		erf_file<<x<<" "<<erf(x)<<" "<< std::erf(x) <<"\n";
	}
	erf_file.close();

    // save data for gamma
	std::ofstream gamma_file("gamma.dat");
	gamma_file<<std::scientific;
	for(double x=xmin;x<=xmax;x+=dx){
		gamma_file<<x<<" "<<sgamma(x)<<" "<< std::tgamma(x) <<"\n";
	}
	gamma_file<<"\n\n";
	
	for(double x=xmin;x<=xmax;x+=dx2){
		gamma_file<<x<<" "<<sgamma(x)<<" "<< std::tgamma(x) <<"\n";
	}
	gamma_file.close();

    // save data for log gamma
    std::ofstream lngamma_file("lngamma.dat");
    lngamma_file<<std::scientific;
    for(double x=xmin;x<=xmax;x+=dx){
        lngamma_file<<x<<" "<<lngamma(x)<<" "<<std::lgamma(x)<<"\n";
    }
    lngamma_file<<"\n\n";

	for(double x=xmin;x<=xmax;x+=dx2){
		lngamma_file<<x<<" "<<lgamma(x)<<" "<< std::lgamma(x) <<"\n";
	}
	lngamma_file.close();

	std::cerr<<"Data written to erf.dat, gamma.dat, and lngamma.dat\n";

	return 0;
}