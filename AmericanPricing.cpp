//
//  AmericanPricing.cpp
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


double american_call_option_trinomial_recursion(int k, int i, double current_price)
{
    int mid = n;
    if (k==n)
        return max(0.0, current_price-K);
    else{
        double up_price, down_price, no_price;
        if (store_array[k][mid+i+1] == -1)
            store_array[k][mid+i+1] = american_call_option_trinomial_recursion(k+1, i+1, current_price*u);
        up_price = store_array[k][mid+i+1];
        
        if (store_array[k][mid+i] == -1)
            store_array[k][mid+i] = american_call_option_trinomial_recursion(k+1, i, current_price);
        no_price = store_array[k][mid+i];
        
        if (store_array[k][mid+i-1] == -1)
            store_array[k][mid+i-1] = american_call_option_trinomial_recursion(k+1, i-1,  current_price/u);
        down_price = store_array[k][mid+i-1];
        
        return max(current_price - K, (pu*up_price + pd*down_price + (1-pu-pd)*no_price) / R);
    }
}


double american_put_option_trinomial_recursion(int k, int i, double current_price)
{
    int mid = n;
    if (k==n)
        return max(0.0, K-current_price);
    else{
        double up_price, down_price, no_price;
        if (store_array[k][mid+i+1] == -1)
            store_array[k][mid+i+1] = american_put_option_trinomial_recursion(k+1, i+1, current_price*u);
        up_price = store_array[k][mid+i+1];
        
        if (store_array[k][mid+i] == -1)
            store_array[k][mid+i] = american_put_option_trinomial_recursion(k+1, i, current_price);
        no_price = store_array[k][mid+i];
        
        if (store_array[k][mid+i-1] == -1)
            store_array[k][mid+i-1] = american_put_option_trinomial_recursion(k+1, i-1,  current_price/u);
        down_price = store_array[k][mid+i-1];
        
        return max(K-current_price, (pu*up_price + pd*down_price + (1-pu-pd)*no_price) / R);
    }
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
    
    // setup memorization array
    store_array = new double*[n];
    for (int i=0; i<n; i++)
        store_array[i] = new double[int(2*n+1)];
    for (int i=0; i<n; i++){
        for (int j=0; j<(2*n+1); j++)
            store_array[i][j] = -1;
    }
    // price american call option
    double american_call = american_call_option_trinomial_recursion(0, 0, S0);
    cout << "Trinomial Price of an American Call Option = " << american_call << endl;
    
    // re-setup memorization array
    for (int i=0; i<n; i++){
        for (int j=0; j<(2*n+1); j++)
            store_array[i][j] = -1;
    }
    // price american put option
    double american_put = american_put_option_trinomial_recursion(0, 0, S0);
    cout << "Trinomial Price of an American Put Option = " << american_put << endl;
    cout << "---------------------------------" << endl;
    
    // put call parity
    cout << "Verifying Put-Call Parity: S+P-C = K*exp(-rf*T)" << endl;
    cout << "S+P-C = " << S0 << "+" << american_put << "-" << american_call << " = " << S0+american_put-american_call << endl;
    cout << "K*exp(-rf*T) = " << K*exp(-rf*T) << endl;
    cout << "Put-Call Parity dose not hold." << endl;
    
    return 0;
}
