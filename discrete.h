#include <time.h>
#include <chrono>
#include <random>
#include <cstdlib>

int number_of_trails, number_of_divisions;
float up_factor, uptick_prob, risk_free_rate, strike_price, barrier_price;
float initial_stock_price, expiration_time, volatility, R;


#pragma once
#pragma once
double get_uniform()
{
	return (((double)random()) / (pow(2.0, 31.0) - 1.0));
}

double max(double a, double b)
{
	return (b < a) ? a : b;
}

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
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
}

double probability_not_breach(double s_t)
{
	if (initial_stock_price <= barrier_price || s_t <= barrier_price)
		return (0);
	else
	{
		double p_b = 1;
		for (int i = 1; i < number_of_divisions; i++)
		{
			double mu, sigma_square;
			mu = initial_stock_price + (expiration_time*i / number_of_divisions)*(s_t - initial_stock_price) / expiration_time;
			sigma_square = (expiration_time*i / number_of_divisions)*(expiration_time - (expiration_time*i / number_of_divisions)) / expiration_time;
			double p = 1 - N((barrier_price - mu) / sqrt(sigma_square));
			p_b = p_b * p;
		}
		return (p_b);
	}
}

double down_and_out_discrete_simulation()
{
	double delta_T = expiration_time / ((float)number_of_divisions);//why we need to convert it into float form here?
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility * sqrt(delta_T);
	double call_price = 0;
	double put_price = 0;
	double call_price_1 = 0;
	double put_price_1 = 0;

	for (int i = 0; i < number_of_trails; i++)
	{
		// by sharing random variables we create 4 paths 
		double current_stock_price1 = initial_stock_price;
		double current_stock_price2 = initial_stock_price;
		double current_stock_price3 = initial_stock_price;
		double current_stock_price4 = initial_stock_price;

		for (int j = 0; j < number_of_divisions; j++)
		{
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			if (current_stock_price1 <= barrier_price)
				current_stock_price1 = 0;
			else
				current_stock_price1 = current_stock_price1 * exp(delta_R + delta_SD * a);

			if (current_stock_price2 <= barrier_price)
				current_stock_price2 = 0;
			else
				current_stock_price2 = current_stock_price2 * exp(delta_R - delta_SD * a);

			if (current_stock_price3 <= barrier_price)
				current_stock_price3 = 0;
			else
				current_stock_price3 = current_stock_price3 * exp(delta_R + delta_SD * b);

			if (current_stock_price4 <= barrier_price)
				current_stock_price4 = 0;
			else
				current_stock_price4 = current_stock_price4 * exp(delta_R - delta_SD * b);


		}
		//std::cout << current_stock_price1 << std::endl << current_stock_price2 << std::endl << current_stock_price3 << std::endl << current_stock_price4 << std::endl;

		call_price = call_price + (max(0, current_stock_price1 - strike_price) + max(0, current_stock_price2 - strike_price) +
			max(0, current_stock_price3 - strike_price) + max(0, current_stock_price4 - strike_price)) / 4;
		put_price = put_price + (max(0, strike_price - current_stock_price1) + max(0, strike_price - current_stock_price2) +
			max(0, strike_price - current_stock_price3) + max(0, strike_price - current_stock_price4)) / 4;

		//std::cout << p_b_1 << std::endl << p_b_2 << std::endl << p_b_3 << std::endl << p_b_4 << std::endl;
		//if (p_b_1 > 1 || p_b_1 < 0 || p_b_2 > 1 || p_b_2 < 0 || p_b_3 > 1 || p_b_3 < 0 || p_b_4 > 1 || p_b_4 < 0)
		//std::cout << "1" << std::endl;

		double R = (risk_free_rate - 0.5*pow(volatility, 2))*expiration_time;
		double SD = volatility * sqrt(expiration_time);
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		double S1 = initial_stock_price * exp(R + SD * a);
		double S2 = initial_stock_price * exp(R - SD * a);
		double S3 = initial_stock_price * exp(R + SD * b);
		double S4 = initial_stock_price * exp(R - SD * b);

		double p_b_1 = 1, p_b_2 = 1, p_b_3 = 1, p_b_4 = 1;
		p_b_1 = probability_not_breach(S1);
		p_b_2 = probability_not_breach(S2);
		p_b_3 = probability_not_breach(S3);
		p_b_4 = probability_not_breach(S4);

		call_price_1 = call_price_1 + (max(0, S1 - strike_price)*p_b_1 +
			max(0, S2 - strike_price)*p_b_2 +
			max(0, S3 - strike_price)*p_b_3 +
			max(0, S4 - strike_price)*p_b_4) / 4;
		put_price_1 = put_price_1 + (max(0, strike_price - S1)*p_b_1 +
			max(0, strike_price - S2)*p_b_2 +
			max(0, strike_price - S3)*p_b_3 +
			max(0, strike_price - S4)*p_b_4) / 4;
		//std::cout << i << std::endl;
	}
	std::cout << "Price of an European Down and Out Discrete Barrier Call Option from Simulation = " << call_price * exp(-risk_free_rate * expiration_time) / number_of_trails << std::endl;
	std::cout << "Price of an European Down and Out Discrete Barrier Call Option with Brownian-Bridge correction = " << call_price_1 * exp(-risk_free_rate * expiration_time) / number_of_trails << std::endl;
	
	std::cout << "Price of an European Down and Out Discrete Barrier Put Option from Simulation = " << put_price * exp(-risk_free_rate * expiration_time) / number_of_trails << std::endl;
	std::cout << "Price of an European Down and Out Discrete Barrier Put Option with Brownian-Bridge correction = " << put_price_1 * exp(-risk_free_rate * expiration_time) / number_of_trails << std::endl;

	return (0);
}
