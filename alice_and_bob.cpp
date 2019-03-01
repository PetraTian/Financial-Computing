// IE523: Financial Computation
// "How to lose as little as possible" by Addona, Wagon and Wilf
// Written by Prof. Sreenivas
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "alice_and_bob.h"
using namespace std;

int main (int argc, char* argv[])
{
    I_have_nothing_apropos_for_this_class x;
    double alice_success_prob, bob_success_prob;
    sscanf (argv[1], "%lf", &alice_success_prob);
    sscanf (argv[2], "%lf", &bob_success_prob);
    
    cout << "Probability of success for Alice = " << alice_success_prob << endl;
    cout << "Probability of success for Bob = " << bob_success_prob << endl;
    
    x.set_probability(alice_success_prob, bob_success_prob);
    
    int optimal = x.search_result();
    if (optimal > 0)
        cout << "The optimal number of coin tosses in each game is " << optimal << endl;
    else {
        cout << "The optimal number of coin tosses in each game exceeds 100... Quitting" << endl;
    }
    
    // simulate games for n from 1 to 50
    int no_of_trials = 1000000;
    int max_n = 30;
    double prob_theoratical[max_n];
    double prob_experiment[max_n];
    for (int n=1; n<=max_n; n++){
        prob_theoratical[n-1] = x.get_theoratical_value(n);
        prob_experiment[n-1] = x.simulated_value(n, no_of_trials);
    }
    
    ofstream outfile("simulated_value.csv");
    for (int i=0; i<sizeof(prob_theoratical)/sizeof(prob_theoratical[0]); i++)
        outfile << prob_theoratical[i] << "," << prob_experiment[i] << endl;
}



