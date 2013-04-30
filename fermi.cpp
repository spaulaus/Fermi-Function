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

//Gamma function for complex argument G(a+bi) returns |G|**2, a la Lanczos!
complex<double> Gamma(double a, double b) {
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
//Beta decay goes from Q to ~0
double Fermi(int z, int a, double w) {
    w /= consts.GetConstant("c").GetValue()*consts.GetConstant("hbar").GetValue();
    double r0 = 1.2e-15; // m
    double R = r0*pow(a, 1./3.);
    double p = sqrt(w*w-1.);
    double fineStructure = consts.GetConstant("fineStructure").GetValue();
    double gamma = sqrt(1.-pow(fineStructure*z,2));
    double y = fineStructure*z*w/p;

    return ( 2*(1.+gamma)*pow(2*p*R, -2*(1.-gamma))*
	     exp(M_PI*y)*norm(Gamma(gamma,y))/norm(Gamma(2*gamma+1.0,0)) ) ;
}

int main() {
    double a = 77;
    double z = 29;
    double qbeta = -10.150;
    double en = 0.269929;
    double step = 0.01;

    //for(double i = qbeta+step; i < 0; i += step)
    //    cout << i << " " << Fermi(z,a,qbeta-i) << endl;
    cout << en << " " << log10(Fermi(z,a,qbeta+en)) << endl;
}
