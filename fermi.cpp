#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <map>
#include <vector>

#include <cmath>

#include "PhysConstants.hpp"

using namespace std;

PhysConstants consts;

//Some of the useful constants to be used
double c = consts.GetConstant("c").GetValue(); // m/s
double hbar = consts.GetConstant("hbar").GetValue()/1e6; // MeV*s
double me = consts.GetConstant("electronMass").GetValue(); // MeV/c/c

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

double CalcFermiFunction(const int &z, const int &a, const double &w) {
    double R = (1.123*pow(a, 1./3.) - 0.941*pow(a, -1./3.))*1e-15*me/hbar/c;
    double p = sqrt(w*w-1); 
    double fineStructure = consts.GetConstant("fineStructure").GetValue();
    double gamma = sqrt(1.-pow(fineStructure*z,2));
    double y = (fineStructure*z*w)/p;

    // cout << "Inside F(z,a,w) = " << R << " " << p << " " << gamma << " " << y 
    //      << " " << norm(Gamma(gamma,y)) << " " << norm(Gamma(2*gamma+1.0,0)) << endl;

    return ( 2*(1.+gamma)*pow(2*p*R, -2*(1.-gamma))*
	     exp(M_PI*y)*norm(Gamma(gamma,y))/norm(Gamma(2*gamma+1.0,0)) ) ;
}

double CalcFermiIntegral(void) {
    //The one I'm interested in
    //double a = 77, z = 29, qbeta = 10.28/me; //77cu; units of MeV
    //From Langanke
    //double a = 63, z = 27, qbeta = (3.661+0.087)/me; // 63co; units of MeV
    //From Gove and Martin
    double a = 120, z = 50; 
    double step = 5e-6;
    
    double q[] = {0.05, 0.5, 5.0};
    for(int i = 0; i < 3; i++) {
        double qbeta = q[i]/me;
        cout << "Qbeta = " << qbeta << endl;
        
        double integral = 0; 
        for(double w = 1+step; w < qbeta; w+=step) {
            double p = sqrt(w*w-1);  
            double diff = pow(qbeta-w,2.0);
            double fermiFunc = CalcFermiFunction(z+1,a,w);
            
            // cout << w << " " << fermiFunc << " " << " " << p << " " 
            //      << diff << " " << endl << endl;
            integral += w*p*diff*fermiFunc;
        }
        integral *= pow(me,4);
        // cout << "f = " << integral << endl 
        //      << "log(f) = " << log10(integral) << endl
        //      << "log(ft) = " << log10(integral*27.4) << endl << endl;
    }
    return(34973);
}

double KanteleFermi(const double &z, const double &e) {
    //Expects energy in kev
    //Some tabulated values for the function
    double m[16][11] = {
        //Values for Electron Decay
        {1, .3269, -.12948, .019016, -.00095348, 0, .0272646, -.0004201, -9.5474E-6, 0, 0}, 
        {10, 2.08738, -.57928, .0404785, .00305364, -.000334, .11416, .043251, -.0033661, 0, 0}, 
        {20, 2.80848, -.404105, -.068659, .019067, -.0010833, .127189, .128762, -.010038, 0, 0}, 
        {30, 3.33244, -.300744, -.11193, .023452, -.0012114, .083368, .244539, -.019246, 0, 0}, 
        {50, 4.38334, -.319865, -.091311, .016894, -.0007709, .047726, .515291, -.0424, 0, 0}, 
        {70, 5.51184, -.38239, -.058665, .0107613, -.000491881, .515219, .738253, -.066362, 0, 0}, 
        {90, 6.76374, -.41807, -.042962, .0086018, -.00051951, 2.52482, .66674, -.076984, 0, 0}, 
        {100, 7.4614, -.437171, -.035496, .0077881, -.00057477, 4.4908, .42503, -.071491, 0, 0}, 
        //Values for Positron Decay
        {1, -.313382, .117663, -.016218, .0007587, 0, -.024786, 5.7346E-5, 2.31865E-5, 0, 0}, 
        {10, -3.84456, 1.534, -.21992, .010662, 0, -.053804, -.041362, .002706, 0, 0}, 
        {20, -10.8276, 5.50759, -1.14606, .111243, -.00421536, .230165, -.15023, .00902366, 0, 0}, 
        {30, -12.1595, 4.7432, -.65041, .029935, 0, .37228, -.21268, .012221, 0, 0}, 
        {50, -22.5153, 8.92152, -1.22276, .056017, 0, -6.54927, 2.2527, -.28514, .011719, 0}, 
        {70, -33.01442, 13.1511, -1.789534, .0810547, 0, -11.6324, 4.19548, -.52875, .0214844, 0}, 
        {90, -42.6523, 16.9538, -2.2726, .10076, 0, -12.1595, 4.6966, -.61494, .025391, 0}, 
        {100, -47.5832, 18.9583, -2.53152, .111491, 0, -16.05257, 6.25352, -.81078, .033203, 0} 
    };

    int g;
    if(e > 2000)
        g = 5.;
    else
        g = 0.;

    int idx;
    for(int i = 0; i < 8; i++) {
        if( (m[i+1][0] - z) >= 0) {
            idx = i;
            break;
        }
    }
    
    double y[16];
    for(int j = idx; j <= idx+1; j++) {
        y[j] = exp(m[j][1 + g] + 
                   m[j][2 + g] * log(e) + 
                   m[j][3 + g] * pow(log(e),2) + 
                   m[j][4 + g] * pow(log(e),3) + 
                   m[j][5 + g] * pow(log(e),4) );
    }
    double f = exp(log(y[idx]) + 
                   ( z - m[idx][0] ) * 
                   ( log(y[idx + 1]) - log(y[idx]) ) / 
                   ( m[idx + 1][0] - m[idx][0] ) );
    return(f);
}

int main() {
    double z = 10;
    vector<double> vals = {2,4,6,8,10,20,40,60,80,100,200,400,600,800,1000,6000};
    for(const auto &i : vals)
        cout << i << " " << KanteleFermi(z, i) << endl;

}
