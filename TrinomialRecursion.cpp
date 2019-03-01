//
//  main.cpp
//  AmericanOptionPricing
//
//  Created by PC on 2018/11/12.
//  Copyright © 2018 PC. All rights reserved.
//

# include <iostream>
# include <cmath>
using namespace std;

// global variables
double T, n, rf, sigma, S0, K;
double R, u, pu, pd;

// this function return the larger value between two passed-in value
double max(double a, double b)
{
    return (b<a)? a:b;
}

double american_call_option_trinomial_recursion(int k, double current_price)
{
    if (k==n)
        return max(0.0, current_price-K);
    else
        return max(current_price-K, (pu*american_call_option_trinomial_recursion(k+1, current_price*u) + pd*american_call_option_trinomial_recursion(k+1, current_price/u) + (1-pu-pd)*american_call_option_trinomial_recursion(k+1,current_price)) / R);
}

double american_put_option_trinomial_recursion(int k, double current_price)
{
    if (k==n)
        return max(0.0, K-current_price);
    else
        return max(K-current_price, (pu*american_put_option_trinomial_recursion(k+1, current_price*u) + pd*american_put_option_trinomial_recursion(k+1, current_price/u) + (1-pu-pd)*american_put_option_trinomial_recursion(k+1, current_price)) / R);
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
    cout << "Recursive Trinomial American Option Pricing" << endl;
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
    double call_price_recursion, put_price_recursion;
    call_price_recursion = american_call_option_trinomial_recursion(0, S0);
    put_price_recursion = american_put_option_trinomial_recursion(0, S0);
    cout << "Trinomial Price of an American Call Option = " << call_price_recursion << endl;
    cout << "Trinomial Price of an American Put Option = " << put_price_recursion << endl;
    cout << "---------------------------------" << endl;
    
    
    
}
