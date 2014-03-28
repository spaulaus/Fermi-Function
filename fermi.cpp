#include <complex>
#include <iostream>
#include <iomanip>

#include <cmath>

#include <PhysConstants.hpp>

using namespace std;

PhysConstants consts;

//Some of the useful constants to be used
// double me = 510.999;
// double fineStructure = 1./137.;

double me = consts.GetConstant("electronMass").GetValue()*1000; // keV/c/c
double fineStructure = consts.GetConstant("fineStructure").GetValue();

//Gamma function for complex argument G(a+bi) returns the complex arguement
//from http://en.wikipedia.org/wiki/Lanczos_approximation
complex<double> Gamma(const double &a, const double &b) {
    double gg = 7.0;
    complex<double> z(a,b), cG, x, t;

    //Coefficients used by the GNU Scientific Library
    double c[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         	   771.32342877765313, -176.61502916214059, 12.507343278686905,
                   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
    //Reflection Formula
    if(z.real() < 0.5) {
        complex<double> tmp = 1.0 - z;
        return(M_PI/ (sin(M_PI*z)*Gamma(tmp.real(), tmp.imag())) );
    } else {
        z -= 1.0;
        x = c[0];

        for(int i = 1; i < gg+2; i++)
            x+=c[i]/(z+double(i));
        t = z + gg + 0.5;
        cG = sqrt(2*M_PI) * pow(t,z+0.5) * exp(-t) * x;
        return(cG);
    }
}

double CalcFermiFunction(const int &z, const int &a, const double &w) {
    //From Royer, NPA, 807, 105-118 (2008) 
    double R = 1.1818*pow(a,1/3) -0.089 + 1.5938 * pow(a,-1/3); 
    //double R = fineStructure*pow(a,1/3) / 2.;
    double p = sqrt(w*w-1); 
    double gamma = sqrt(1.-pow(fineStructure*z,2));
    double y = (fineStructure*z*w)/p;

    return ( 2*(1.+gamma)*pow(2*p*R, 2*gamma-2)*
             exp(M_PI*y)*
             ( norm(Gamma(gamma,y)) / norm(Gamma(2*gamma+1.0,0)) ) );
}

int main() {
    double z = 3, a = 6, q = 3505.21;
    for(double p = 0.0; p < 10.; p+= 0.1) {
        double W = sqrt(p*p + 1);
        double T = (W - 1)*me;
        if(T > q)
            break;
        else
            cout << setprecision(10) << p << " " << T << " " 
                 << CalcFermiFunction(z, a, W)*(p/W) << endl;
    }
}
