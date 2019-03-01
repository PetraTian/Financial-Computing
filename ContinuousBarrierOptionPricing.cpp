//
//  ContinuousBarrierOptionPricing.cpp
//  BarrierOptionPricing
//
//  Created by PC on 2018/12/1.
//  Copyright © 2018 PC. All rights reserved.
//

# include <iostream>
# include <cmath>
# include "black_scholes.h"
using namespace std;


double T, rf, sigma, S0, K, B;
int num_trials, num_divisions;
double **path;

double theoratical_down_and_out_call_price()
{
    // get Black-Scholes price for vanilla call option
    double time_sqrt = sqrt(T);
    double d1 = (log(S0/K)+rf*T)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    double c = S0*N(d1) - K*exp(-rf*T)*N(d2);
    
    // get theoratical value of barrier call option
    double barrier_call_price;
    double lamda = (rf+0.5*pow(sigma,2.0)) / pow(sigma,2.0);
    if (B<=K){
        double y = log(pow(B,2.0)/(S0*K)) / (sigma*time_sqrt) + lamda*sigma*time_sqrt;
        double cdi = S0*pow(B/S0,2.0*lamda)*N(y) - K*exp(-rf*T)*pow(B/S0,2.0*lamda-2.0)*N(y-sigma*time_sqrt);
        barrier_call_price = c-cdi;
    }
    else{
        double x1 = log(S0/B)/(sigma*time_sqrt) + lamda*sigma*time_sqrt;
        double y1 = log(B/S0)/(sigma*time_sqrt) + lamda*sigma*time_sqrt;
        barrier_call_price = S0*N(x1)-K*exp(-rf*T)*N(x1-sigma*time_sqrt)-S0*pow(B/S0,2.0*lamda)*N(y1)+K*exp(-rf*T)*pow(K/S0,2.0*lamda-2.0)*N(y1-sigma*time_sqrt);
    }
    return barrier_call_price;
}

double theoratical_down_and_out_put_price()
{
    // get Black-Scholes price for vanilla put option
    double time_sqrt = sqrt(T);
    double d1 = (log(S0/K)+rf*T)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    double p = K*exp(-rf*time_sqrt)*N(-d2) - S0*N(-d1);
    
    double barrier_put_price;
    double lamda = (rf+0.5*pow(sigma,2.0)) / pow(sigma,2.0);
    double y = log(pow(B,2.0)/(S0*K)) / (sigma*time_sqrt) + lamda*sigma*time_sqrt;
    double x1 = log(S0/B)/(sigma*time_sqrt) + lamda*sigma*time_sqrt;
    double y1 = log(B/S0)/(sigma*time_sqrt) + lamda*sigma*time_sqrt;
    // get theoratical value for barrier put option
    if (B<=K){
        double pdi = -S0*N(-x1)+K*exp(-rf*T)*N(-x1+sigma*time_sqrt)+S0*pow(B/S0,2.0*lamda)*(N(y)-N(y1))-K*exp(-rf*T)*pow(B/S0,2.0*lamda-2.0)*(N(y-sigma*time_sqrt)-N(y1-sigma*time_sqrt));
        barrier_put_price = p-pdi;
    }
    else
        barrier_put_price = 0;
    
    return barrier_put_price;
}

int main(int argc, const char * argv[]) {
    sscanf(argv[1], "%lf", &T);
    sscanf(argv[2], "%lf", &rf);
    sscanf(argv[3], "%lf", &sigma);
    sscanf(argv[4], "%lf", &S0);
    sscanf(argv[5], "%lf", &K);
    sscanf(argv[6], "%d", &num_trials);
    sscanf(argv[7], "%d", &num_divisions);
    sscanf(argv[8], "%lf", &B);
    
    // initialize path array
    path = new double*[num_trials];
    for (int i=0; i<num_trials; i++){
        path[i] = new double[num_divisions+1];
        path[i][0] = S0;
    }
    
    // calculate parameter for motion in each segment
    double delta_T = T / num_divisions;
    double delta_R = (rf - 0.5*pow(sigma,2.0)) * delta_T;
    double delta_SD = sigma * sqrt(delta_T);
    
    // do monte-carlo simulation
    int num_samplings = num_trials / 4;
    for (int i=0; i<num_samplings; i++){
        for (int j=0; j<num_divisions; j++){
            // generate unit-normal random numbers using Box-Muller method
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            // simulate one step forward
            path[i*4][j+1] = path[i*4][j] * exp(delta_R + delta_SD * a);
            path[i*4+1][j+1] = path[i*4+1][j] * exp(delta_R - delta_SD * a);
            path[i*4+2][j+1] = path[i*4+2][j] * exp(delta_R + delta_SD * b);
            path[i*4+3][j+1] = path[i*4+3][j] * exp(delta_R - delta_SD * b);
        }
    }
    
    
    // get call option price sum at time T
    double call_option_price_sum = 0;
    double put_option_price_sum = 0;
    for (int i=0; i<num_trials; i++){
        for (int j=1; j<num_divisions+1; j++){
            if (path[i][j] <=B)
                break;
            if (j==num_divisions){
                call_option_price_sum += max(0.0, path[i][j] - K);
                put_option_price_sum += max(0.0, K - path[i][j]);
            }
        }
    }
    
    // get average adjusted price
    double adjusted_call_option_price_sum = 0;
    double adjusted_put_option_price_sum = 0;
    for (int i=0; i<num_trials; i++){
        double ST = path[i][num_divisions];
        if (ST > B){
            double pc = exp(-(2.0*log(S0/B)*log(ST/B))/(pow(sigma,2.0)*T));
            adjusted_call_option_price_sum += max(0, ST-K)*(1-pc);
            adjusted_put_option_price_sum += max(0, K-ST)*(1-pc);
        }
    }
    
    // get theoratical value
    double theoratical_call_option_price = theoratical_down_and_out_call_price();
    double theoratical_put_option_price = theoratical_down_and_out_put_price();
    
    // output
    cout << "----------------------------------------------------" << endl;
    cout << "European Down-and-Out Continuous Barrier Options Pricing via Monte-Carlo Simulation" << endl;
    cout << "Expiration Time(Years) = " << T << endl;
    cout << "Risk Free Interest Rate = " << rf << endl;
    cout << "Volatility (%age of stock value) = " << sigma*100 << endl;
    cout << "Initial Stock Price = " << S0 << endl;
    cout << "Strick Price = " << K << endl;
    cout << "Barrier Price = " << B << endl;
    cout << "Number of Trials = " << num_trials << endl;
    cout << "Number of Divisions = " << num_divisions << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "Average Call Price by Explicit Simulation = " << exp(-rf*T)*(call_option_price_sum / double(num_trials)) << endl;
    cout << "The Call Price Using the (1-pc)-Adjustment Term = " << exp(-rf*T)*(adjusted_call_option_price_sum / double(num_trials)) << endl;
    cout << "Theoratical Call Price = " << theoratical_call_option_price << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "Average Put Price by Explicit Simulation = " << exp(-rf*T)*(put_option_price_sum / double(num_trials)) << endl;
    cout << "The Put Price Using the (1-pc)-Adjustment Term = " << exp(-rf*T)*(adjusted_put_option_price_sum / double(num_trials)) << endl;
    cout << "Theoratical Put Price = " << theoratical_put_option_price << endl;
    
    return 0;
}
