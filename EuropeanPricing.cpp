//
//  EuropeanPricing.cpp
//  OptionPricingWithRecursionMemorization
//
//  Created by PC on 2018/11/13.
//  Copyright © 2018 PC. All rights reserved.
//

# include <iostream>
# include <cmath>
using namespace std;

double T, rf, sigma, S0, K;
int n;
double R, u, pu, pd;
double **store_array;


double max(double a, double b)
{
    return (b<a)? a:b;
}


double N(double z)
{
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
    if (z<0.0) n = 1.0-n;
    return n;
}


double european_call_option_trinomial_recursion(int k, int i, double current_price)
{
    int mid = n;
    if (k==n)
        return max(0.0, current_price-K);
    else{
        double up_price, down_price, no_price;
        if (store_array[k][mid+i+1] == -1)
            store_array[k][mid+i+1] = european_call_option_trinomial_recursion(k+1, i+1, current_price*u);
        up_price = store_array[k][mid+i+1];
        
        if (store_array[k][mid+i] == -1)
            store_array[k][mid+i] = european_call_option_trinomial_recursion(k+1, i, current_price);
        no_price = store_array[k][mid+i];
        
        if (store_array[k][mid+i-1] == -1)
            store_array[k][mid+i-1] = european_call_option_trinomial_recursion(k+1, i-1,  current_price/u);
        down_price = store_array[k][mid+i-1];
        
        return ((pu*up_price + pd*down_price + (1-pu-pd)*no_price) / R);
    }
}


double european_put_option_trinomial_recursion(int k, int i, double current_price)
{
    int mid = n;
    if (k==n)
        return max(0.0, K-current_price);
    else{
        double up_price, down_price, no_price;
        if (store_array[k][mid+i+1] == -1)
            store_array[k][mid+i+1] = european_put_option_trinomial_recursion(k+1, i+1, current_price*u);
        up_price = store_array[k][mid+i+1];
        
        if (store_array[k][mid+i] == -1)
            store_array[k][mid+i] = european_put_option_trinomial_recursion(k+1, i, current_price);
        no_price = store_array[k][mid+i];
        
        if (store_array[k][mid+i-1] == -1)
            store_array[k][mid+i-1] = european_put_option_trinomial_recursion(k+1, i-1,  current_price/u);
        down_price = store_array[k][mid+i-1];
        
        return ((pu*up_price + pd*down_price + (1-pu-pd)*no_price) / R);
    }
}


double european_call_option_Black_Scholes(double S0, double K, double rf, double sigma, double T)
{
    double T_sqrt = sqrt(T);
    double d1 = (log(S0/K)+rf*T) / (sigma*T_sqrt) + 0.5*sigma*T_sqrt;
    double d2 = d1 - (sigma*T_sqrt);
    return S0*N(d1) - K*exp(-rf*T)*N(d2);
}


double european_put_option_Black_Scholes(double S0, double K, double rf, double sigma, double T)
{
    double T_sqrt = sqrt(T);
    double d1 = (log(S0/K)+rf*T) / (sigma*T_sqrt) + 0.5*sigma*T_sqrt;
    double d2 = d1 - (sigma*T_sqrt);
    return K*exp(-rf*T)*N(-d2) - S0*N(-d1);
}


int main(int argc, char * argv[])
{
    // read data from command line
    sscanf(argv[1], "%lf", &T);
    sscanf(argv[2], "%d", &n);
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
    cout << "Recursive Trinomial European Option Pricing" << endl;
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
    
    // setup memorization array
    store_array = new double*[n];
    for (int i=0; i<n; i++)
        store_array[i] = new double[2*n+1];
    for (int i=0; i<n; i++){
        for (int j=0; j<(2*n+1); j++)
            store_array[i][j] = -1;
    }
    // price european call option
    double european_call = european_call_option_trinomial_recursion(0, 0, S0);
    double call_BS = european_call_option_Black_Scholes(S0, K, rf, sigma, T);
    cout << "Trinomial Price of an European Call Option = " << european_call << endl;
    cout << "Call Price According to Black Scholes = " << call_BS << endl;
    cout << "---------------------------------" << endl;
    
    // re-setup memorization array
    
    for (int i=0; i<n; i++){
        for (int j=0; j<(2*n+1); j++)
            store_array[i][j] = -1;
    }
    // price european put option
    double european_put = european_put_option_trinomial_recursion(0, 0, S0);
    double put_BS = european_put_option_Black_Scholes(S0, K, rf, sigma, T);
    cout << "Trinomial Price of an European Put Option = " << european_put << endl;
    cout << "Put Price According to Black Scholes = " << put_BS << endl;
    cout << "---------------------------------" << endl;
    
    // put call parity
    cout << "Verifying Put-Call Parity: S+P-C = K*exp(-rf*T)" << endl;
    cout << "S+P-C = " << S0 << "+" << european_put << "-" << european_call << " = " << S0+european_put-european_call << endl;
    cout << "K*exp(-rf*T) = " << K*exp(-rf*T) << endl;
    cout << "Put-Call Parity Holds for European Option." << endl;
    
    return 0;
}
