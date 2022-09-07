#pragma once
#include <cmath>
#include "typedefs.hpp"
static inline real heaviside(real x){
	if(x <= 0) return 0;
	return 1.0; 
}
static inline real h_x(real x, real xmin, real xmax, real sigma, real xc){
    real eval = 0.0;
    real k, k1, k2, invk;
    real sigma2 = sigma*sigma;
    k = sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma2));
    invk = 1/k;
    k1 = invk*sqrt(0.5*M_PI*sigma2);
    k2 = invk*exp(-0.5*(xc*xc)/sigma2);
    eval = (k1 * erf((xmax-x)/(sqrt(2)*sigma)) - k2*(xmax-x) - 0.5)*heaviside(xc - fabs(xmax-x))
    + (k1 * erf((x-xmin)/(sqrt(2)*sigma)) - k2*(x-xmin) - 0.5)*heaviside(xc - fabs(x-xmin))
    + heaviside(xc + 0.5*(xmax-xmin) - fabs(x - 0.5*(xmin+xmax)));
    return eval;
}