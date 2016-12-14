//#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <string>
#include "string.h"

#define DEBUG_SPACING 10000

#define DELTAF false
#define LASTX false
#define ITERX true

#define DELTAF_PRECISION 0.0000001
#define PRECISION 0.0000001
#define MAX_ITERATIONS 1000000
#define GAMMA 0.01

using namespace std;

double f1(double x, double y)
{
	return log(1 + log(x) * log(x) + log(y) * log(y));
}

double f2(double x, double y)
{
	return log(1 + x * x + (x * x - y) * (x * x - y));
}

void df1(double x, double y, double* out) {
	// x   (2 log(x))/(x (log^2(x) + log^2(y) + 1))
	out[0] = (2 * log(x)) / (x * (1 + log(x) * log(x) + log(y) * log(y)));
	// y   (2 log(y))/(y (log^2(x) + log^2(y) + 1))
	out[1] = (2 * log(y)) / (y * (1 + log(y) * log(y) + log(x) * log(x)));
}

void df2(double x, double y, double* out) {
	// x   (2 (2 x^3 - 2 x y + x))/((x^2 - y)^2 + x^2 + 1)
	out[0] = (2 * (2 * x*x*x - 2 * x*y + x)) / (pow(x*x - y, 2) + x*x + 1);
	// y  -(2 (x^2 - y))/((x^2 - y)^2 + x^2 + 1)
	out[1] = -1 * (2 * (x*x - y)) / (pow(x*x - y, 2) + x*x + 1);
}

double goldenSectionSearch(double precision, double in_a, double in_b, double *x, double *dx, double(*f) (double, double))
{
	const double theta1 = ((3-sqrt(5)) / 2);
	const double theta2 = 1 - theta1;

	// initialize bounds
	double a = in_a;
	double c = (in_a + in_b) / 2; //center
	double b = in_b;

	while (f(x[0] + c*dx[0], x[1] + c*dx[1]) < f(x[0] + c*dx[0], x[1] + c*dx[1]))
	{
		a = c;
		c = b;
		b = 2*b;
	}

	double u = a + theta1 * (b - a);
	double v = a + theta2 * (b - a);

	while (b - a > precision)
	{
		if (f(x[0] + u*dx[0], x[1] + u*dx[1]) < f(x[0] + v*dx[0], x[1] + v*dx[1]))
		{
			b = v;
			v = u;
			u = a + theta1 * (b - a);
		}
		else
		{
			a = u;
			u = v;
			v = a + theta1 * (b - a);
		}
	}
	return (u + v) / 2;
}

void gradient(string name, double precision, double gamma, int max_iterations, double *out, void(*d_function)(double, double, double*), double(*function) (double, double))
{
	int current_iteration = 0;
	double grad[2] = { 0, 0 };
	double tmp[2] = { 0, 0 };
	double last_out[2] = { 0, 0 };

	cout << "Starting Gradient descent for "<< name <<".\nParameters:" << endl;
	cout << "Gamma: " << gamma << endl;
	cout << "Precision: " << precision << endl;
	cout << "Max number of iterations: " << max_iterations << endl;
	cout << "Starting..." << endl;
	//main loop
	while (true)
	{
		if (current_iteration % DEBUG_SPACING == 0)
		{
			cout << "Iteration: " << current_iteration << ". Result: (" << out[0] << ", " << out[1] << ")" << endl;
		}
		// calculate new derivate
		d_function(out[0], out[1], grad);
		tmp[0] += -gamma * grad[0];
		tmp[1] += -gamma * grad[1];

		// golden section search
		double t;
		if (name == "F1")
		{
			// we cant use GSS on f1 because it is not unimodal
			t = 0.01;
		}
		else
		{
			t = goldenSectionSearch(precision, 0, 10, out, tmp, function);
		}

		last_out[0] = out[0];
		last_out[1] = out[1];

		out[0] = out[0] + t * tmp[0];
		out[1] = out[1] + t * tmp[1];

		current_iteration += 1;
		if ((current_iteration > max_iterations) && ITERX)
		{
			break;
		}
		if ((abs(function(out[0], out[1]) - function(last_out[0], last_out[1])) < DELTAF_PRECISION) && DELTAF)
		{
			break;
		}
		if (((abs(out[0] - last_out[0]) < precision) && (abs(out[1] - last_out[1]) < precision)) && LASTX)
		{
			break;
		}
	}
	cout << "Proccess finalized for " << name << ". Result: (" << out[0] << ", " << out[1] << ") Value: " << function(out[0], out[1]) << endl;
	cout << "Total iterations: " << current_iteration << endl;
}

int main()
{
	// starting point
	double new_arr[2];
	// step size
	double gamma = GAMMA;
	// precision needed for stop
	double precision = PRECISION ;
	// max iterations
	int max_iterations = MAX_ITERATIONS;

	new_arr[0] = 1.5;
	new_arr[1] = 1.5;
	gradient("F1", precision, gamma, max_iterations, new_arr, df1, f1);

	new_arr[0] = 0.5;
	new_arr[1] = 0.5;
	gradient("F2", precision, gamma, max_iterations, new_arr, df2, f2);

	cout << "Exiting..." << endl;
	cin.ignore();

	return 0;
}
