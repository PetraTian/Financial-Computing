/*
 *  black_scholes.h
 *  Monte Carlo European
 *
 *  Created by Ramavarapu Sreenivas on 10/24/14.
 *  Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved.
 *
 */
#include "normdist.h"
#include <random>
#include <chrono>

double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};

double option_price_delta_call_black_scholes(const double& S,     // spot price
                                             const double& K,     // Strike (exercise) price,
                                             const double& r,     // interest rate
                                             const double& sigma, // volatility
                                             const double& time){  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double delta = N(d1);
    return delta;
};

double option_price_delta_put_black_scholes(const double& S, // spot price
                                            const double& K, // Strike (exercise) price,
                                            const double& r,  // interest rate
                                            const double& sigma,
                                            const double& time) {
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double delta = -N(-d1);
    return delta;
}

double max(double a, double b) {
    return (b < a )? a:b;
}

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
// unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
// default_random_engine generator (seed);
// initially just use default_random_engine generator;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);

// u.i.i.d. generator
double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

// unit-normal i.i.d. generator
double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}
