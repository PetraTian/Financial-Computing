#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include"discrete.h"

using namespace std;

int main(int argc, char* argv[])
{

	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%f", &risk_free_rate);
	sscanf(argv[3], "%f", &volatility);
	sscanf(argv[4], "%f", &initial_stock_price);
	sscanf(argv[5], "%f", &strike_price);
	sscanf(argv[6], "%d", &number_of_trails);
	sscanf(argv[7], "%d", &number_of_divisions);
	sscanf(argv[8], "%f", &barrier_price);

	cout << "--------------------------------------" << endl;
	cout << "European Down-and-Out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trails = " << number_of_trails << endl;
	cout << "Number of Divisions = " << number_of_divisions << endl;
	cout << "--------------------------------------" << endl;

	cout << "--------------------------------------" << endl;
	down_and_out_discrete_simulation();
	cout << "--------------------------------------" << endl;
}
