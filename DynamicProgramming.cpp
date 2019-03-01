//
//  DynamicProgramming.cpp
//  AmericanOptionPricing
//
//  Created by PC on 2018/11/12.
//  Copyright © 2018 PC. All rights reserved.
//

# include <iostream>
# include <cmath>
# include "newmat.h"
using namespace std;

double T, n, rf, sigma, S0, K;
double R, u, pu, pd;

double max(double a, double b)
{
    return (a>b)? a:b;
}

double american_call_option_trinomial_dyn_prog()
{
    double state_num_at_T = 2*n+1;
    
    // create transition matrix
    Matrix transition_matrix(state_num_at_T,state_num_at_T);
    for (int i=1; i<=state_num_at_T; i++)
        transition_matrix = 0.0;
    transition_matrix(1,1) = 1-pu;
    transition_matrix(1,2) = pu;
    transition_matrix(state_num_at_T, state_num_at_T-1) = pd;
    transition_matrix(state_num_at_T, state_num_at_T) = 1-pd;
    for (int i=2; i<state_num_at_T; i++){
        transition_matrix(i,i-1) = pd;
        transition_matrix(i,i+1) = pu;
        transition_matrix(i,i) = 1-pu-pd;
    }
    
    // option price matrix at time T
    Matrix V_t(state_num_at_T, 1);
    for (int i=1; i<=state_num_at_T; i++)
        V_t(i,1) = max(0.0, S0*pow(u, -n+i-1)-K);
    
    // use transition matrix to calculate initial option price
    Matrix one_step_backward_value(state_num_at_T, 1);
    Matrix option_exercise_now_value(state_num_at_T, 1);
    for (int i=n; i>0; i--){
        one_step_backward_value = transition_matrix * V_t / R;
        for (int j=1; j<=state_num_at_T; j++)
            option_exercise_now_value(j,1) = max(0.0, S0*pow(u, -n+j-1)-K);
        for (int j=1; j<=state_num_at_T; j++)
            V_t(j,1) = max(one_step_backward_value(j,1),option_exercise_now_value(j,1));
    }
    
    return V_t(n+1,1);
}

double american_put_option_trinomial_dyn_prog()
{
    double state_num_at_T = 2*n+1;
    
    // create transition matrix
    Matrix transition_matrix(state_num_at_T,state_num_at_T);
    for (int i=1; i<=state_num_at_T; i++)
        transition_matrix = 0.0;
    transition_matrix(1,1) = 1-pu;
    transition_matrix(1,2) = pu;
    transition_matrix(state_num_at_T, state_num_at_T-1) = pd;
    transition_matrix(state_num_at_T, state_num_at_T) = 1-pd;
    for (int i=2; i<state_num_at_T; i++){
        transition_matrix(i,i-1) = pd;
        transition_matrix(i,i+1) = pu;
        transition_matrix(i,i) = 1-pu-pd;
    }
    
    // option price matrix at time T
    Matrix V_t(state_num_at_T, 1);
    for (int i=1; i<=state_num_at_T; i++)
        V_t(i,1) = max(0.0, K-S0*pow(u, -n+i-1));
    
    // use transition matrix to calculate initial option price
    Matrix one_step_backward_value(state_num_at_T, 1);
    Matrix option_exercise_now_value(state_num_at_T, 1);
    for (int i=n; i>0; i--){
        one_step_backward_value = transition_matrix * V_t / R;
        for (int j=1; j<=state_num_at_T; j++)
            option_exercise_now_value(j,1) = max(0.0, K-S0*pow(u, -n+j-1));
        for (int j=1; j<=state_num_at_T; j++)
            V_t(j,1) = max(one_step_backward_value(j,1),option_exercise_now_value(j,1));
    }
    
    return V_t(n+1,1);
}


int main(int argc, char * argv[])
{
    // read data from command line
    sscanf(argv[1], "%lf", &T);
    sscanf(argv[2], "%lf", &n);
    sscanf(argv[3], "%lf", &rf);
    sscanf(argv[4], "%lf", &sigma);
    sscanf(argv[5], "%lf", &S0);
    sscanf(argv[6], "%lf", &K);
    
    // calculate R, u, pu, pd, pn
    R = exp(rf*T/n);
    u = exp(sigma*sqrt(2*T/n));
    pu = pow((sqrt(R)-1/sqrt(u))/(sqrt(u)-1/sqrt(u)),2);
    pd = pow((sqrt(u)-sqrt(R))/(sqrt(u)-1/sqrt(u)),2);
    
    // print data
    cout << "American Option Pricing by Trinomial-Model-Inspired Dynamic Programming" << endl;
    cout << "Expiration Time(Years) = " << T << endl;
    cout << "Number of Divisions = " << n << endl;
    cout << "Risk Free Interest Rate = " << rf << endl;
    cout << "Volatility (%age of stock value) = " << sigma*100 << endl;
    cout << "Initial Stock Price = " << S0 << endl;
    cout << "Strick Price = " << K << endl;
    cout << "---------------------------------" << endl;
    cout << "R = " << R << endl;
    cout << "Up factor = " << u << endl;
    cout << "Uptick Probability = " << pu << endl;
    cout << "Downtick Probability = " << pd << endl;
    cout << "Notick Probability = " << 1-pu-pd << endl;
    cout << "---------------------------------" << endl;
    
    // pricing option with recursion
    double call_price_dyn, put_price_dyn;
    call_price_dyn = american_call_option_trinomial_dyn_prog();
    put_price_dyn = american_put_option_trinomial_dyn_prog();
    cout << "Trinomial Price of an American Call Option = " << call_price_dyn << endl;
    cout << "Trinomial Price of an American Put Option = " << put_price_dyn << endl;
    cout << "---------------------------------" << endl;
    
}
