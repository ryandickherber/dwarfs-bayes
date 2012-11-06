//continuum code written by Jim Buckley
//output modified by Ryan Dickherber to be compatible with dwarfs-bayes

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define UU  0
#define SS  1
#define TT  2
#define DD  3
#define CC  4
#define BB  5 //b/bbar?
#define WW  6
#define ZZ  7
#define gg  8
#define TAUTAU  9 //tau+/tau-?

double m_chi=1000.0;   /* Mass of neutralino in GeV */
int ichan=BB;
#define NE 31
double Erange[NE];
double Eincrement=0.15;
double Estart=100.0;
#define NM 30
double Mlist[NM];
double Mincrement=0.15;
double Mstart=115.0;

/*
double alpha1[] = {0.95, 0.0, 1.1, 0.0, 0.0, 1.0, 0.73, 0.73, 0.0, 0.0} ;
double alpha2[] = {6.5, 0.0, 15.1, 0.0, 0.0, 10.7, 7.76, 7.76, 0.0, 0.0} ;
double  a[] = {   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,   -1.31}; 
double  b[] = { 0.0063,   0.04,  -0.45, 0.0063,   0.17,   0.37,  -0.95,  -0.83,   0.48,    6.94};
double  c[] = {  -8.62,  -8.84, -19.05,  -8.62, -10.23, -16.05,  -9.86,-11.175, -20.51,   -4.93};
double  d[] = {   8.53,   2.77,  21.96,   8.53,   2.13,  18.01,   6.25, 6.5902,  24.42,   -0.51};
double  e[] = {  -9.73,  -7.71, -15.18,  -9.73,  -7.00, -19.50,  -4.37,-3.6468, -19.56,   -4.53};
double  eta[]={    1.0,    1.0,    2.0,    1.0,    1.0,    1.0,    2.0,    2.0,   1.0} ;
*/
// 110704 JB I continue to use the same fit functions, but changed parameters to
// better fit values and mass scaling found by M. Vivier through refitting Fornengo's functions
// to Pythia simulations.
double alpha1[]=
 {0.95,  0.0,  1.1,  0.0,  0.0,   1.0,   0.73,   0.73,   0.0,   0.0} ;
double alpha2[]=
 {6.5,   0.0, 15.1,  0.0,  0.0,   10.7,  7.76,   7.76,   0.0,   0.0} ;
//                  uu,     ss,     tt,     dd,     cc,     bb,     ww,     zz,     gg,   tautau
double  a[] = {   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,   -1.5,  -1.36,   -1.4,   -1.5,   -1.29}; 
double  b[] = { 0.0063,   0.04,  -0.45, 0.0063,   0.17,   1.05,   0.50,   0.47,   0.48,    7.83};
double  c[] = {  -8.62,  -8.84, -19.05,  -8.62, -10.23, -17.80,  -12.6,  -14.4, -20.51,   -2.70};
double  d[] = {   8.53,   2.77,  21.96,   8.53,   2.13,   12.3,   12.1,  16.30,  24.42,  -0.009};
double  e[] = {  -9.73,  -7.71, -15.18,  -9.73,  -7.00,  -1.86,  -9.86, -13.60, -19.56,   -5.06};
double  eta[]={    1.0,    1.0,    2.0,    1.0,    1.0,    1.0,    2.0,    2.0,    1.0,     1.0} ;
    

/*
    In the following functions for the continuum spectrum, we 
    always return spectra as e dn/de.  
 
    9/21/09 Note that in the past I changed from natural log of e (ln(e)) to
    log_10(e).  This introduced some errors since dn/dln(e)=e dn/de but
    dlog_10(e) = 1/(e ln(10)) and dn/dlog_10(e) = e ln(10) dn/de
        so, correcting an earlier error, this function and the other functions
    that calculate the continuum spectrum actually return...

        dn/dx = m_chi*dn/de, so 
        dn/dlog_10(e) = e ln(10) dn/de
         = (e/m_chi) ln(10) dn/dx =  x ln(10) dn/dx
*/

/*
    Function: continuum

    Return the differential spectrum of continuum photons from
    dark matter annihilation as x dn/dx or e dn/de.  Note that
    this functional form is taken from Fornengo, Pieri and Scopel,
    Phys. Rev. D., 70, 103529 where the formula is given by
    dn/dx (with x = E/m_chi) and where the gamma-ray multiplicity
    (e.g. 2 or 1 for different channels) is given in a separate factor
    eta.
*/
double continuum(double loge, int ichan) { /* loge is log_10 of energy in GeV */
    double  x ;
    double  dndx ;
    double  cont ;
    double  av,bv,cv,dv,ev ;
    double  logmchi ;

    x = pow(10.0,loge)/m_chi ;
    if(x>1.0) return(0.0) ;     // Force a cutoff for E> m_chi
    av=a[ichan]; bv=b[ichan]; cv=c[ichan]; dv=d[ichan]; ev=e[ichan] ;
    if(ichan == TAUTAU) {
      logmchi = log(m_chi) ;
      bv= 6.08 + 0.27*logmchi ;
          cv= -6.9 + 5.2e-4*m_chi-4.2e-08*m_chi*m_chi+1.6e-12*m_chi*m_chi*m_chi;
          ev= -5.25 + 0.36*logmchi - 0.040*logmchi*logmchi;
      dndx = pow(x,av)*(bv*x + cv*x*x + dv*pow(x,3))
                 * exp(ev*x) ;
    }
    else {
          if(ichan == BB) {
        logmchi = log(m_chi) ;  // CHECK IF LOG means ln or log_10 in MV's NOTES
        av = -1.46 + (4.26e-2)*logmchi - (4.37e-3)*logmchi*logmchi ;
          }
      dndx = eta[ichan]*pow(x,av)*exp(bv+cv*x +dv*x*x + ev*pow(x,3)) ;
    }
    cont = x*dndx ; 
    if (cont < 1.0e-20) cont=1.0e-20 ;  
    return(cont) ;
}

double continuum_div_e(double e) {
    return (1.0/log(10))*continuum(log10(e), ichan)/e;
}

double f(double x, void * params) {
    return continuum_div_e(x);
}

//integrated continuum
double int_cont(double e_start, double e_end) {
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function F;
    F.function = &f;

    gsl_integration_qags (&F, e_start, e_end,
            0, 1e-7, 1000, w, &result, &error);

    gsl_integration_workspace_free(w);
    return result;
}

int main()
{
    printf("ichan=%i\n", ichan);
    printf("NM=%i\n", NM);
    printf("Mstart=%f #GeV\n", Mstart);
    printf("Mint=%f\n", Mincrement);
    printf("NE=%i\n", NE);
    printf("Estart=%f #GeV\n", Estart);
    printf("Eint=%f\n", Eincrement);
    printf("integrated=[\\\n");
    for (int i=0; i<NE; i++) {
        if (i==0)
            Erange[i]=Estart;
        else
            Erange[i]=Erange[i-1]+Erange[i-1]*Eincrement;
    }
    for (int i=0; i<NM; i++) {
        if (i==0)
            Mlist[i]=Mstart;
        else
            Mlist[i]=Mlist[i-1]+Mlist[i-1]*Mincrement;
    }
    for (int j=0; j<NM; j++) {
        printf("    [\\\n");
        for (int i=0; i<NE-1; i++) {
                m_chi=Mlist[j];
                double r;
                r=int_cont(Erange[i], Erange[i+1]);
                //printf("%i\t%f\t%i\t%f\t%e\n",
                //        i, Erange[i], j, m_chi, r);
                printf("        [%i,%f,%i,%f,%f,%e],\\\n",
                        j, m_chi, i, Erange[i], Erange[i+1], r);
        }
        printf("    ],\\\n");
    }
    printf("]\n");

    return 0;
}
