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
double c = consts.GetConstant("c").GetValue(); // m/s
double hbar = consts.GetConstant("hbar").GetValue()/1e6; // MeV*s
double me = consts.GetConstant("electronMass").GetValue()/c/c; // MeV/c/c

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

//Courtesy of M. G. Bertolli
//units on w = MeV and Fermi = unitless
double CalcFermiFunction(const int &z, const int &a, const double &w) {
    double R = (2.908e-3*pow(a, 1./3.) - 2.437e-3*pow(a, -1./3.))*(hbar/me/c);
    double p = sqrt(w*w-pow(me,2.))/c; 
    double fineStructure = consts.GetConstant("fineStructure").GetValue();
    double gamma = sqrt(1.-pow(fineStructure*z,2));
    double y = (fineStructure*z*w)/p/c;

    // cout << "Inside F(z,a,w) = " << R << " " << p << " " << gamma << " " << y 
    //      << " " << norm(Gamma(gamma,y)) << " " << norm(Gamma(2*gamma+1.0,0)) << endl;

    return ( 2*(1.+gamma)*pow(2*p*R*hbar/me/me/c/c/c/c, -2*(1.-gamma))*
	     exp(M_PI*y)*norm(Gamma(gamma,y))/norm(Gamma(2*gamma+1.0,0)) ) ;
}

double CalcFermiIntegral(void) {
    //The one I'm interested in
    //double a = 77, z = 29, qbeta = 10.28; //77cu; units of MeV
    //From Langanke
    //double a = 63, z = 27, qbeta = (3.661+0.087)/me/c/c; // 63co; units of MeV
    //From Gove and Martin
    double a = 120, z = 50, qbeta = 5; // units of MeV    
    double step = 5e-5;

    cout << "Qbeta = " << qbeta << endl;

    double integral = 0; 
    for(double w = 1+step; w < qbeta; w+=step) {
        double pe = sqrt(w*w-pow(me,2.))/c;  
        double diff = pow(qbeta-w,2.0);
        double fermiFunc = CalcFermiFunction(z+1,a,w);
        
        // cout << w << " " << fermiFunc << " " << " " << pe << " " 
        //      << diff << " " << endl << endl;
        integral += w*pe*diff*fermiFunc;
    }
    return(integral/me/me/me/me);
}

int main() {
    double f = CalcFermiIntegral();
    cout << "f = " << f << endl 
         << "log(f) = " << log10(f) << endl
         << "log(ft) = " << log10(f*27.4) << endl;
}
