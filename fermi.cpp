#include <complex>
#include <iostream>
#include <iomanip>
#include <ostream>


#include <cmath>

using namespace std;

//Gamma function for complex argument G(a+bi) returns |G|**2, a la Lanczos!
complex<double> Gamma(double a, double b) {
    double gg = 7.0;
    complex<double> z(a,b), cG, x, t;
    double c[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
		   771.32342877765313, -176.61502916214059, 12.507343278686905,
		   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

    if(z.real() < 0.5) {
        return (M_PI / (sin(M_PI*z)*Gamma(1.0-z.real(), z.imag())) );
    } else {
        z -= 1.0;
        x = c[0];
	
	for(int i = 1; i < 9; i++)
	    x+=c[i]/(z+double(i));
	t = z + gg + 0.5;
	cG = sqrt(2*M_PI) * pow(t,z+0.5) * exp(-t) * x;
	return(cG);
    }
}

//Courtesy of M. G. Bertolli
//units on w and Fermi in MeV
//Beta decay goes from 1 to Q
double Fermi(int z, int a, double w) {
    double alpha = 0.0072973525698, r0 = 1.2;
    double R = r0*pow(a, 1./3.);
    double p = sqrt(w*w-1.);
    double gamma = sqrt(1.-pow(alpha*z,2));
    double y = alpha*z*w/p;
    
    return ( 2*(1.+gamma)*pow(2*p*R, -2*(1.-gamma))*
	     exp(M_PI*y)*norm(Gamma(gamma,y))/norm(Gamma(2*gamma+1.0,0)) ) ;
}

int main() {
    complex<double> a(0.5,0), b(1.5,0), c(2.5,0);
    cout << setprecision(12) << "G(2.5, 0) = " 
	 << Gamma(c.real(), c.imag()) << endl;

    cout << setprecision(12) << "Fermi(119,339,1.2) = " << Fermi(119,339,1.2) << endl;
	//<< "Fermi(29,73,1.2) = " << Fermi(29,73,1.2) << endl;
}
    
	
