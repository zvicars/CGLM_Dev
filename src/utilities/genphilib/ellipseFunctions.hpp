#pragma once
#include "../../tools/typedefs.hpp"
#include "../../tools/stlmath.hpp"
#include <cmath>
#include <limits>

namespace ep{
    const double EPSILON = 1e-4;
    const int MAX_ITER = 1000;

    //for integrating within a single cell, will just place a plane tangent to the ellipse after the first fit and use that for subsequent evaluations
    //need to set/unset integrate flags appropriately
    Vec3<real> LAST_NORMAL;
    Vec3<real> LAST_POS;
    bool INTEGRATE_FLAG = 0;
    bool LAST_NORMAL_SET = 0;

    real getPlanarDist(Vec3<real> x){
        return dot(LAST_NORMAL, x - LAST_POS);
    }

}

// Function to evaluate the ellipsoid equation
static inline real ellipsoidEquation(const Vec3<real>& x, const Vec3<real>& a) {
    return ((x[0] * x[0]) / (a[0] * a[0])) + ((x[1] * x[1]) / (a[1] * a[1])) + ((x[2] * x[2]) / (a[2] * a[2])) - 1.0;
}

static inline real gradEll(const Vec3<real>& x, const Vec3<real>& a, std::size_t i){
    return 2.0*x[i] / (a[i] * a[i]);
}

static inline real gradDist(const Vec3<real>& x, std::size_t i){
    return 2.0 * x[i];
}

static inline Vec3<real> guessPoint(const Vec3<real>& x, const Vec3<real>& a) {
    Vec3<real> retVal;
    real xnorm = norm2(x);
    for(int i = 0; i < 3; i++){
        retVal[i] = x[i]*a[i] / xnorm;
    }
    return retVal;
}

// Function to find the nearest point on the ellipsoid
static inline Vec3<real> nearestPointOnEllipsoid(const Vec3<real>& p, const Vec3<real>& a) {
    real lambda = 1;
    real dLambda;
    real LR = 0.1;
    real c = 0.1;
    Vec3<real> x, dx;
    //std::cout << p[0] << "  " << p[1] << "  " << p[2] << std::endl;
    x = guessPoint(p, a);
    //x.fill(0.0);
    int it;
    for (it = 0; it < ep::MAX_ITER; ++it) {
        //apply gradient descent to the lagrangian
        Vec3<real> gradL;
        real gradL_Lambda = ellipsoidEquation(x, a);
        auto rvec = x - p;
        //real L = lambda * gradL_Lambda + dot(rvec,rvec);
        for(int i = 0; i < 3; i++){
            gradL[i] = gradDist(rvec, i) + (lambda*gradEll(x, a, i)) + (c*gradL_Lambda*gradEll(x, a, i)) ;
        }
        for(int i = 0; i < 3; i++){
            dx[i] = -LR*gradL[i];
        }
        dLambda =  gradL_Lambda;
        x = x + dx;
        lambda += dLambda;

        //std::cout << x[0] << "  " << x[1] << "  " << x[2] << std::endl;
        //std::cout << dLambda + dx[0] + dx[1] + dx[2]<< std::endl;
        //std::cin.get();

        if (std::abs(dLambda) + std::abs(dx[0]) + std::abs(dx[1]) + std::abs(dx[2]) < ep::EPSILON) {
            break;
        }
    }
    if(ep::INTEGRATE_FLAG == 1 && ep::LAST_NORMAL_SET == 0){
        ep::LAST_NORMAL = {gradEll(x, a, 0), gradEll(x, a, 1),  gradEll(x, a, 2)};
        ep::LAST_POS    = x;
        ep::LAST_NORMAL_SET = 0;
    }
    return x;
}