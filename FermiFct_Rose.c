#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double ALPHA = (1./137.);
const double mec2 = 510.999;				// in keV
const double Qb = 3505.21;

const double Z = 3;					// Atomic number
const double A = 6;					// Mass number
const double trTyp = -1.;				// -1 for electron; +1 for positron

const double PMIN = 0.0;
const double PMAX = 10.0;
const double PSTP = 0.1;

/**********************************************************************/
// Downloaded from http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html
#ifndef DCOMPLEX
struct dcomplex_ {
    double re;
    double im;
};
#define DCOMPLEX struct dcomplex_
#define DREAL(x) (x).re
#define DIMAG(x) (x).im
#define DCMPLX(x,y,z) (z).re = x, (z).im = y
#endif

/**********************************************************************/
// Downloaded from http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.htm

/* complex Gamma function in double precision */
DCOMPLEX cdgamma(DCOMPLEX x)
{
    DCOMPLEX y;
    double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;
    
    xr = DREAL(x);
    xi = DIMAG(x);
    if (xr < 0) {
        wr = 1 - xr;
        wi = -xi;
    } else {
        wr = xr;
        wi = xi;
    }
    ur = wr + 6.00009857740312429;
    vr = ur * (wr + 4.99999857982434025) - wi * wi;
    vi = wi * (wr + 4.99999857982434025) + ur * wi;
    yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 
        0.293729529320536228;
    yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
    ur = vr * (wr + 4.00000003016801681) - vi * wi;
    ui = vi * (wr + 4.00000003016801681) + vr * wi;
    vr = ur * (wr + 2.99999999944915534) - ui * wi;
    vi = ui * (wr + 2.99999999944915534) + ur * wi;
    yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
    yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
    ur = vr * (wr + 2.00000000000603851) - vi * wi;
    ui = vi * (wr + 2.00000000000603851) + vr * wi;
    vr = ur * (wr + 0.999999999999975753) - ui * wi;
    vi = ui * (wr + 0.999999999999975753) + ur * wi;
    yr += ur * 10.5400280458730808 + vr;
    yi += ui * 10.5400280458730808 + vi;
    ur = vr * wr - vi * wi;
    ui = vi * wr + vr * wi;
    t = ur * ur + ui * ui;
    vr = yr * ur + yi * ui + t * 0.0327673720261526849;
    vi = yi * ur - yr * ui;
    yr = wr + 7.31790632447016203;
    ur = log(yr * yr + wi * wi) * 0.5 - 1;
    ui = atan2(wi, yr);
    yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
    yi = ui * (wr - 0.5) + ur * wi;
    ur = yr * cos(yi);
    ui = yr * sin(yi);
    yr = ur * vr - ui * vi;
    yi = ui * vr + ur * vi;
    if (xr < 0) {
        wr = xr * 3.14159265358979324;
        wi = exp(xi * 3.14159265358979324);
        vi = 1 / wi;
        ur = (vi + wi) * sin(wr);
        ui = (vi - wi) * cos(wr);
        vr = ur * yr + ui * yi;
        vi = ui * yr - ur * yi;
        ur = 6.2831853071795862 / (vr * vr + vi * vi);
        yr = ur * vr;
        yi = ur * vi;
    }
    DCMPLX(yr, yi, y);
    return y;
}

/**********************************************************************/
double FermiF(double SZ, double W, double A)
{
    // for SZ<0 : repulsion (beta+)
    // for SZ>0 : attraction (beta-)
    
    DCOMPLEX z,Gam;
    double gr,gi;
    double p,gam,R,y;
    double s,t;
    
    p = sqrt(W*W - 1);
    gam = sqrt(1. - pow(ALPHA*SZ,2.));
    R = ALPHA * pow(A,1./3.) / 2.;			// from M.E. Rose
    y = ALPHA*SZ*W/p;
    
    // from M.E. Rose, Eq.(9c), p.280 in "Beta- and Gamma-ray Spectroscopy
    // K. Siegbahn, 1955
    
    s = 2.*(1. + gam)*pow(2.*p*R,2.*(gam-1.))*exp(M_PI*y);
    
    DCMPLX(gam,y,z);
    Gam = cdgamma(z);
    gr = DREAL(Gam);
    gi = DIMAG(Gam);
    
    //t = lgamma(gam) - lgamma(2.*gam + 1.);
    //t = exp(2.*t);
    
    t = gr*gr + gi*gi; 
    t /= exp(2.*lgamma(2.*gam + 1.));
    
    return(s*t);
}

/**********************************************************************/
int main(int argc, char** argv)
/***
    Driver for the calculation of the Fermi funcion F(+-Z,W).
    Prints the values of the "modified Fermi function" G=(p/W)F
    to compare with the results in the tables of M.E. Rose et al.
    (Appendix II, p.887, Siegbahn 1955)
    Uses the complex Gamma function cdgamma.
***/
{
    double p,W,SIGZ,T;
    
    for(p = PMIN; p <= PMAX; p += PSTP)
        {
            W = sqrt(p*p + 1);
            SIGZ = -trTyp * Z;
            T = (W - 1.)*mec2;
            if(T > Qb)
                break;
            printf("%5.2lf  %5.2lf  %e\n",p,T,FermiF(SIGZ,W,A)*(p/W));
        }
    return(0);
}
/**********************************************************************/
