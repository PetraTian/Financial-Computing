/*
 *  alice_and_bob.h
 *  Loosing as little as possible
 *
 *  Created by Ramavarapu Sreenivas on 9/2/12.
 *  Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved.
 *
 */

#ifndef ALICE_AND_BOB
#define ALICE_AND_BOB

#include <cmath>
#include <random>

using namespace std;

class I_have_nothing_apropos_for_this_class
{
private:
    double alice_probability, bob_probability;
    
    // private member function: uniform RV generator
    double get_uniform()
    {
        // write the appropriate code here
        return (((float) random()) / (pow(2.0, 31.0)-1.0));
    }
    
    // private member function: nCi (i.e. n-take-i)
    int take(int n, int i)
    {
        // write a **RECURSIVE** implementation of n-take-i.
        // If you made it non-recurisive (i.e. n!/((n-i)!i!)) -- it
        // will take too long for large sizes
        if (i == 0)
            return 1;
        if (i == 1)
            return n;
        return (take(n,i-1) * (n-(i-1)) / i);
    }
    
    // this routine implements the probability that Alice has more
    // heads than Bob after n-many coin tosses
    double theoretical_value(double q, double p, int n)
    {
        // implement equation 1.1 of Addona-Wagon-Wilf paper
        // r represents the number of heads that Bob gets
        // s represents the number of heads that Alice gets
        
        // store take number to save computational time
        double take_array[n];
        for (int i=0; i<=n; i++)
            take_array[i] = -1;
        // calculate theoretical result
        double func_result = 0;
        for (int r=0; r<=n; r++){
            double temp = 0;
            for (int s=r+1; s<=n; s++){
                if(take_array[s] == -1)
                    take_array[s] = take(n,s);
                temp += take_array[s]*pow(q,s)*pow(1-q,n-s);
            }
            if (take_array[r] == -1)
                take_array[r] = take(n,r);
            func_result += take_array[r]*pow(p,r)*pow(1-p,n-r)*temp;
        }
        return func_result;
    }
    
public:
    // public function:
    void set_probability(double alice_p, double bob_p)
    {
        alice_probability = alice_p;
        bob_probability = bob_p;
    }
    
    // get theoratical probability for N(p,q)
    double get_theoratical_value(int n){
        return theoretical_value(alice_probability, bob_probability, n);
    }
    
    // probability of Alice winning the game.
    double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
    {
        int no_of_wins_for_alice = 0;
        for (int i = 0; i < no_of_trials; i++)
        {
            int number_of_heads_for_alice = 0;
            int number_of_heads_for_bob = 0;
            for (int j = 0; j < number_of_coin_tosses_in_each_game; j++)
            {
                if (get_uniform() < alice_probability)
                    number_of_heads_for_alice++;
                if (get_uniform() < bob_probability)
                    number_of_heads_for_bob++;
            }
            if (number_of_heads_for_alice > number_of_heads_for_bob)
                no_of_wins_for_alice++;
        }
        return (((double) no_of_wins_for_alice)/((double) no_of_trials));
    }
    
    int search_result()
    {
        // implememt a discrete-search procedure for the optimal n-value.
        // start with n = 1 and find the discrete-value of n that has
        // the largest probability for Alice winning.  Why would this work?
        // See Theorem 2.2 of the paper for the reason!
        double take_array[2];
        take_array[0] = take(1,0);
        take_array[1] = take(1,1);
        double last_f = theoretical_value(alice_probability, bob_probability, 1);
        double diff;
        for (int n=2; n<=100; n++){
            double take_array[n+1];
            for (int i=0; i<=n; i++)
                take_array[i] = take(n,i);
            double current_f = theoretical_value(alice_probability, bob_probability, n);
            diff = current_f - last_f;
            if (diff > 0)
                last_f = current_f;
            else
                return n-1;
        }
        return -1;
    }
};
#endif
