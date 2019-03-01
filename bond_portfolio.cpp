// Written by Prof. Sreenivas for IE523: Financial Computing
// Modified by Peichen Tian

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include "lp_lib.h"

using namespace std;

const double ERROR = 1e-10;
int number_of_cash_flows;  //number of bonds used for debt immune
vector <double> price_list; //PV of bonds
vector <int> maturity_list; //bond maturity
vector <vector <double> > cash_flow_list; //cash flow of all bonds
double debt_obligation_amount;
double time_when_debt_is_due;
double debt_PV; //present value of debt
vector <double> yield_to_maturity; //bonds yield to maturity
double ytm_mean; //mean of all bonds ytm
vector <double> duration;
vector <double> duration_lp; //duration used for lp
vector <double> convexity;
vector <double> convexity_lp; //convexity used for lp
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;
double max_convexity; //max convexity to get from immunization
double *optimal_weight = new double[number_of_cash_flows]; //optimal weight of bonds(result from lp)


double my_function(vector <double> cash_flow, double price, int maturity, double rate){
    // This function computes f(r) in page 2 of lesson 3 lecture notes
    double func_r_second_part = 0;
    for (int t=1; t<=maturity; t++){
        func_r_second_part += cash_flow[t-1] * pow(1+rate, maturity-t);
    }
    double func = price * pow(1+rate, maturity) - func_r_second_part;
    return func;
}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate){
    // This function computes f'(r) in the bottom of page 2 of lesson 3 lecture notes
    double der_func_second_part = 0;
    for (int t=1; t<=maturity; t++)
        der_func_second_part += cash_flow[t-1] * (maturity-t) * pow(1+rate, maturity-t-1);
    double der_func = maturity * price * pow(1+rate, maturity) - der_func_second_part;
    return der_func;
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate){
    // This function finds the (only) +ve root of f(r) of page 2 of
    // lesson 3 using Newton-Raphson method.
    // When difference between rate and updated rate is smaller than the predefined ERROR,
    // we have got the solution of the function.
    double rate1 = rate;
    do{
        rate = rate1;
        double func = my_function(cash_flow, price, maturity, rate);
        double der_func = derivative_function(cash_flow, price, maturity, rate);
        rate1 = rate - func / der_func;
    }while(abs(rate - rate1) > ERROR);
    return rate1;
}


void get_average_ytm(){
    // This function computes the average of ytm
    double ytm_sum = 0;
    for (int i=0; i<number_of_cash_flows; i++){
        ytm_sum += yield_to_maturity[i];
    }
    ytm_mean = ytm_sum / number_of_cash_flows;
}

void get_yield_to_maturity_list(){
    // This function computes yield to maturity list for all bonds
    // and also call the function to compute the average of all ytm to get ytm_mean.
    // Yield to maturity is the discount that match future cash flows of bond
    // with its present value.
    for (int i=0; i<number_of_cash_flows; i++){
        double price = price_list[i];
        double maturity = maturity_list[i];
        vector <double> cash_flow = cash_flow_list[i];
        double initial_rate = 0.05;
        double ytm = Newton_Raphson(cash_flow, price, maturity, initial_rate);
        yield_to_maturity.push_back(ytm);
    }
    get_average_ytm();
}

void present_value_of_debt(){
    // This function compute PV of future debt obligation using the average-value-of-the-YTMs
    debt_PV = debt_obligation_amount / pow((1+ytm_mean), time_when_debt_is_due);
}

void get_percentage_of_cash_flow_to_meet_debt_obligation(){
    // This function calls the function to compute present value of debt obligation,
    // match the first term of bond cash flow with debt obligation.
    // Percentage_of_cash_flow_to_meet_debt_obligation is the discount rate to match
    // present value of bond with present value of debt obligation (What percentage
    // of bonds to buy to immune future debt if interest doesn't change).
    present_value_of_debt();
    for (int i=0; i<number_of_cash_flows; i++){
        double percentage = debt_PV / price_list[i];
        percentage_of_cash_flow_to_meet_debt_obligation.push_back(percentage);
    }
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate){
    // This function computes the duration of a single bonds(cash flow).
    // duration = sum(t*Pt/(1+r)^t, t=1:maturity) / bond_present_value, r is ytm of the bond
    double numerator = 0;
    for (int t=1; t<=maturity; t++)
        numerator += t * cash_flow[t-1] / pow(1+rate, t);
    double dur = numerator / price;
    return dur;
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate){
    // This function computes the convexity of a single bonds(cash flow).
    // convexity = sum(t*(t+1)*Pt/(1+r)^(t+2), t=1:maturity) / bond_present_value
    // r is ytm of the bond.
    double numerator = 0;
    for (int t=1; t<=maturity; t++)
        numerator += t * (t+1) * cash_flow[t-1] / pow(1+rate, t+2);
    double conv = numerator / price;
    return conv;
}

void get_duration_convexity_list(){
    // This function computes the duration and convexity of all bonds,
    // and also the duration and convexity used for lp_solve.
    // Use ytm for particular bond as the discount rate to compute duration and conv.
    for (int i=0; i<number_of_cash_flows; i++){
        double dur = get_duration(cash_flow_list[i],price_list[i],maturity_list[i],
                                  yield_to_maturity[i]);
        duration.push_back(dur);
        double dur_lp = dur / percentage_of_cash_flow_to_meet_debt_obligation[i];
        duration_lp.push_back(dur_lp);
        double conv = get_convexity(cash_flow_list[i],price_list[i],
                                    maturity_list[i], yield_to_maturity[i]);
        convexity.push_back(conv);
        double conv_lp = conv / percentage_of_cash_flow_to_meet_debt_obligation[i];
        convexity_lp.push_back(conv_lp);
    }
}


void get_data(char* argv[]){
    // This function reads the data from the file identified on the command-line,
    // and call several functions to compute global variables.
    // Modified from: https://blog.csdn.net/sunshineacm/article/details/78068987
    ifstream input_file(argv[1]);
    if (input_file.is_open()){
        input_file >> number_of_cash_flows;
        input_file.get();
        for (int i=0; i<number_of_cash_flows; i++){
            vector <double> temp;
            cash_flow_list.push_back(temp);
            
            string line;
            stringstream ss;
            getline(input_file, line);
            ss.clear();
            ss.str(line);
            int count = 0;
            while(true){
                double value;
                ss >> value;
                if (ss.fail())
                    break;
                count++;
                if (count == 1)
                    price_list.push_back(value);
                else if (count == 2)
                    maturity_list.push_back(value);
                else
                    cash_flow_list[i].push_back(value);
            }
        }
        input_file >> debt_obligation_amount;
        input_file >> time_when_debt_is_due;
    }
    get_yield_to_maturity_list();
    get_percentage_of_cash_flow_to_meet_debt_obligation();
    get_duration_convexity_list();
}

void print_data(char *filename){
    // This function prints out the data
    cout << "Input File: " << filename << endl;
    cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
    cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
    for (int i = 0; i < number_of_cash_flows; i++)
    {
        cout << "---------------------------" << endl;
        cout << "Cash Flow #" << i+1 << endl;
        cout << "Price = " << price_list[i] << endl;
        cout << "Maturity = " << maturity_list[i] << endl;
        cout << "Percentage of Face Value that would meet the obligation = " <<
        percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
        cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
        cout << "Duration = " << duration[i] << endl;
        cout << "Duration (to be used in LP formulation below = " << duration_lp[i] << endl;
        cout << "(Note) " << duration[i] << " = " << duration_lp[i] << " * "  << percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
        cout << "Convexity = " << convexity[i] << endl;
        cout << "Convexity (to be used in LP formulation below = " << convexity_lp[i] << endl;
        cout << "(Note) " << convexity[i] << " = " << convexity_lp[i] << " * " << percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
    }
    cout << "***************************" << endl;
    cout << "Average YTM = " << ytm_mean << endl;
    cout << "Present value of debt " << debt_PV << endl;
}

void print_optimal_portfolio(double ret, lprec *lp){
    // This function print out the fraction of PV of each bonds to buy,
    // and the amount of each bonds' cash flow to buy
    cout << "***************************" << endl;
    print_lp(lp);
    if(ret == 0){
        cout << "Largest Convexity we can get is:" << max_convexity << endl;
        cout << "Optimal portfolio" << endl;
        // This print out the fraction of bond Present Value to buy for each bond
        for (int i=0; i<number_of_cash_flows; i++){
            double lamda = optimal_weight[i];
            cout << "%Cash Flow:" << i+1 << "  " << lamda << endl;
        }
        cout << "That is, buy" << endl;
        // This print out the real amount of cash flow to buy for each bond
        for (int i=0; i<number_of_cash_flows; i++){
            if (optimal_weight[i] != 0){
                double value = price_list[i] * optimal_weight[i];
                cout << "$" << value << " of Cash Flow#" << i+1 << endl;
            }
        }
    }
    else //if the value of ret is not 0, then no solution
        cout << "There is no portfolio that meets the duration constraint of 10 years." << endl;
}

void get_optimal_portfolio(){
    // This function use lp_solve to solve fraction of bonds cash flow to buy in order to
    // immune future debt.
    double price_array[number_of_cash_flows+1];
    price_array[0] = 0;
    for (int i=1; i<=number_of_cash_flows; i++)
        price_array[i] = price_list[i-1];
    
    double duration_array[number_of_cash_flows+1];
    duration_array[0] = 0;
    for (int i=1; i<=number_of_cash_flows; i++)
        duration_array[i] = duration_lp[i-1];
    
    // Because lp only minimize objective function, in order to get the maximize
    // objective function, we let lp solve the minimize of the opposite objective function
    double convexity_array[number_of_cash_flows+1];
    convexity_array[0] = 0;
    for (int i=1; i<=number_of_cash_flows; i++)
        convexity_array[i] = -convexity_lp[i-1];
    
    lprec *lp;
    lp = make_lp(0,number_of_cash_flows); //five variables lp
    set_verbose(lp,3);
    
    // Match the first term of the bond(present value)
    add_constraint(lp, price_array, EQ, debt_PV);
    
    // Match the second term of cash flow (weighted duration of bonds should add up to
    // "maturity" of debt.
    add_constraint(lp, duration_array, EQ, time_when_debt_is_due);
    
    // Objective is to maximize weighted convexity of bonds portfolio.
    set_obj_fn(lp, convexity_array);
    
    // Positive solution is inherently satisfied in lp_solve.
    // ret is the return value from lp_solve(whether lp solve have found solution or not)
    double ret = solve(lp);
    
    // The max convexity we want to get should be the opposite number of the value
    // given by lp_solve.
    max_convexity = -get_objective(lp);
    // Get the optimize weight.
    get_variables(lp, optimal_weight);
    print_optimal_portfolio(ret, lp);
}


int main (int argc, char* argv[])
{
    if (argc == 1)
        cout << "Input filename missing" << endl;
    else {
        get_data(argv);
        print_data(argv[1]);
        get_optimal_portfolio();
    }
    return (0);
}


