#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <map>

#include <cmath>

#include "PhysConstants.hpp"

using namespace std;

PhysConstants consts;

//Some of the useful constants to be used
double c = consts.GetConstant("c").GetValue();
double hbar = consts.GetConstant("hbar").GetValue();
double eMass = consts.GetConstant("electronMass").GetValue() / (c*c);

//Gamma function for complex argument G(a+bi) returns |G|**2, a la Lanczos!
complex<double> Gamma(const double &a, const double &b) {
    double gg = 7.0;
    complex<double> z(a,b), cG, x, t;
    double c[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
		   771.32342877765313, -176.61502916214059, 12.507343278686905,
		   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

    if(real(z) < 0.5)
	z = 1. - z;
    z -= 1.0;
    x = c[0];
    
    for(int i = 1; i < 9; i++)
	x+=c[i]/(z+double(i));
    t = z + gg + 0.5;
    cG = sqrt(2*M_PI) * pow(t,z+0.5) * exp(-t) * x;
    return(cG);
}

//Courtesy of M. G. Bertolli - that fucking dick. 
//units on w = MeV and Fermi = unitless
double CalcFermiFunction(const int &z, const int &a, const double &w) {
    double energy = w;
    energy /= c*hbar;
    double r0 = 1.2e-15; // m
    double R = r0*pow(a, 1./3.);
    double p = sqrt(w*w-1.);
    double fineStructure = consts.GetConstant("fineStructure").GetValue();
    double gamma = sqrt(1.-pow(fineStructure*z,2));
    double y = fineStructure*z*w/p;

    return ( 2*(1.+gamma)*pow(2*p*R, -2*(1.-gamma))*
	     exp(M_PI*y)*norm(Gamma(gamma,y))/norm(Gamma(2*gamma+1.0,0)) ) ;
}

double CalcFermiIntegral(void) {
    //Some of the initial parameters
    double a = 76;
    double z = 30;
    double qbeta = 10.150;
    double step = 0.01;
    
    double integral = 0; 
    for(double i = 1+step; i < qbeta; i+=step) {
        double fermiFunc = CalcFermiFunction(z,a,-qbeta+i);
        double pe = sqrt(2*eMass*i);
        double diff = pow(qbeta-i,2.0);

        cout << i << " " << fermiFunc << " " << " " << pe << " " 
             << diff << " " << endl;
        integral += fermiFunc*pe*i*pow(qbeta-i, 2.0); 
    }

    cout << "Ingegral = " << integral << "  1/emass**5 = " << pow(eMass,5.) 
         << endl;
    return(integral/pow(eMass,5.));
}

int main() {
    double en = 0.269929;

    //for(double i = qbeta+step; i < 0; i += step)
    //    cout << i << " " << Fermi(z,a,qbeta-i) << endl;
    cout << en << " " << log10(CalcFermiIntegral()) << endl;
}
