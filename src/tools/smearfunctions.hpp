#pragma once
#include <cmath>
#include "typedefs.hpp"
static inline real heaviside(real x){
	if(x <= 0) return 0;
	return 1.0; 
}
static inline real phi_x(real x, real sigma, real xc){
    real sigma2 = sigma*sigma;
    real k = sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma2));
    real invk = 1.0/k;
    return invk*( exp(-x*x/(2*sigma2) - exp(xc*xc/(2*sigma2)) ))*heaviside(xc - fabs(x));
}
static inline real indus_k(real sigma, real xc){
    return sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma*sigma));
}
static inline real indus_k1(real k, real sigma){
    return (1.0/k)*sqrt(0.5*M_PI*sigma*sigma);
}
static inline real indus_k2(real k, real sigma, real xc){
    return (1.0/k)*exp(-0.5*(xc*xc)/(sigma*sigma));
}
static inline real h_x(real x, real xmin, real xmax, real sigma, real xc){
    real eval = 0.0;
    real k, k1, k2, invk;
    real sigma2 = sigma*sigma;
    k = indus_k(sigma, xc);
    k1 = indus_k1(k, sigma);
    k2 = indus_k2(k, sigma, xc);
    eval = (k1 * erf((xmax-x)/(sqrt(2)*sigma)) - k2*(xmax-x) - 0.5)*heaviside(xc - fabs(xmax-x))
    + (k1 * erf((x-xmin)/(sqrt(2)*sigma)) - k2*(x-xmin) - 0.5)*heaviside(xc - fabs(x-xmin))
    + heaviside(xc + 0.5*(xmax-xmin) - fabs(x - 0.5*(xmin+xmax)));
    return eval;
}
static inline real dh_x(real x, real xmin, real xmax, real sigma, real xc){
    return phi_x(xmin-x, sigma, xc)-phi_x(xmax-x, sigma, xc);
}
static inline real h_r(real x, real xmax, real sigma, real xc){
    real eval = 0.0;
    real k, k1, k2, invk;
    real sigma2 = sigma*sigma;
    k = sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma2));
    invk = 1/k;
    k1 = invk*sqrt(0.5*M_PI*sigma2);
    k2 = invk*exp(-0.5*(xc*xc)/sigma2);
    eval = (k1 * erf((xmax-x)/(sqrt(2)*sigma)) - k2*(xmax-x) - 0.5)*heaviside(xc - fabs(xmax-x))
    + heaviside(xc + xmax-x);
    return eval;
}

static inline real dh_r(real x, real xmin, real xmax, real sigma, real xc){
    return -phi_x(xmax-x, sigma, xc);
}
