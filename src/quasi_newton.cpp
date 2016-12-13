#include <iostream>
#include <cmath>
#include <cstring>
#include "string.h"

#define DEBUG false

using namespace std;

double f1 (double x, double y)
{
  return log(1 + log(x) * log(x) + log(y) * log(y));
}

double f2 (double x, double y)
{
  return log(1 + x * x + (x * x - y) * (x * x - y));
}

void df1(double x, double y, double* out) {
  // x   (2 log(x))/(x (log^2(x) + log^2(y) + 1))
  out[0] = (2 * log (x)) / (x * (1 + log(x) * log (x) + log(y) * log(y)));
  // y   (2 log(y))/(y (log^2(x) + log^2(y) + 1))
  out[1] = (2 * log (y)) / (y * (1 + log(y) * log (y) + log(x) * log(x)));
}

void df2(double x, double y, double* out) {
  // x   (2 (2 x^3 - 2 x y + x))/((x^2 - y)^2 + x^2 + 1)
  out[0] = (2*(2*x*x*x - 2*x*y + x))/(pow(x*x - y, 2) + x*x + 1);
  // y  -(2 (x^2 - y))/((x^2 - y)^2 + x^2 + 1)
  out[1] = -1 * (2*(x*x - y))/(pow(x*x - y, 2) + x*x + 1);
}

void Hf1 (double x, double y, double* out)
{
  //xx  -(2 (log(x) (log^2(y) + 1) + log^3(x) + log^2(x) - log^2(y) - 1))/(x^2 (log^2(x) + log^2(y) + 1)^2)
  out[0] = -1*(2*(log(x)*(pow(log(y),2) + 1) + pow(log(x),3) + pow(log(x),2) - pow(log(y),2) - 1))/(x*x*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //xy  -(4 log(x) log(y))/(x y (log^2(x) + log^2(y) + 1)^2)
  out[1] = -1*(4*log(x)*log(y))/(x*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //yx  -(4 log(x) log(y))/(x y (log^2(x) + log^2(y) + 1)^2)
  out[2] = -1*(4*log(x)*log(y))/(x*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //yy  -(2 (log^2(x) (log(y) - 1) + log^3(y) + log^2(y) + log(y) - 1))/(y^2 (log^2(x) + log^2(y) + 1)^2)
  out[3] = -1*(2*(pow(log(x),2)*(log(y) - 1) + pow(log(y),3) + pow(log(y),2) + log(y) - 1))/(y*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
}

void Hf2 (double x, double y, double* out)
{
  //xx    (2 (-2 x^6 + x^4 (2 y - 1) + x^2 (2 y^2 + 4 y + 5) - 2 y^3 + y^2 - 2 y + 1))/(x^4 + x^2 (1 - 2 y) + y^2 + 1)^2
  out[0] = (2*(-2*pow(x,6) + pow(x,4)*(2*y - 1) + x*x*(2*y*y + 4*y + 5) - 2*y*y*y + y*y - 2*y + 1))/pow(pow(x,4) + x*x*(1 - 2*y) + y*y + 1, 2);
  //xy    (4 x (x^4 - 2 x^2 y + y^2 - y - 1))/(x^4 + x^2 (1 - 2 y) + y^2 + 1)^2
  out[1] = (4*x*(pow(x,4) - 2*x*x*y + y*y - y - 1))/pow(pow(x,4) + x*x*(1 - 2*y) + y*y + 1, 2);
  //yx    (4 x (x^4 - 2 x^2 y + y^2 - y - 1))/(x^4 + x^2 (1 - 2 y) + y^2 + 1)^2
  out[2] = (4*x*(pow(x,4) - 2*x*x*y + y*y - y - 1))/pow(pow(x,4) + x*x*(1 - 2*y) + y*y + 1, 2);
  //yy   -(2 (x^4 - x^2 (2 y + 1) + y^2 - 1))/(x^4 + x^2 (1 - 2 y) + y^2 + 1)^2
  out[3] = -1 * (2*(pow(x,4) + x*x*(2*y - 1) + y*y - 1))/pow(pow(x,4) + x*x*(2*y + 1) + y*y + 1, 2);
}

void invert2dmatrix (double* original, double* inverted)
{
  double idet = 1 / (original[0]*original[3]-original[1]*original[2]);
  inverted[0] = original[3] * idet;
  inverted[1] = original[1] * idet * -1;
  inverted[2] = original[2] * idet * -1;
  inverted[3] = original[0] * idet;
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

void quasiNewton (string name, double precision, int max_iterations, double *out, double (*function) (double, double), void (*d_function)(double, double, double*), void (*h_function)(double, double, double*))
{
  int current_iteration = 0;
  double grad [2] = {0, 0};
  double est_inv_hess [4] = {0, 0, 0, 0};
  double dk [2] = {0, 0};
  double tmp [2] = {0, 0};
  double p [2] = {0, 0};
  double q [2] = {0, 0};
  out[0] = out[1] = 1.3;

  cout << "Starting Quasi-Newtonn method.\nParameters:" << endl;
  cout << "Precision: " << precision << endl;
  cout << "Max number of iterations: " << max_iterations << endl;
  cout << "Starting..." << endl;

  //initialize H
  h_function(out[0], out[1], est_inv_hess);

  //main loop
  while (current_iteration < max_iterations)
  {
    //tmp = x[k-1]

    // p = x[k] - x[k-1]
    p[0] = out[0] - tmp[0];
    p[1] = out[1] - tmp[1];

    d_function(tmp[0], tmp[1], grad);

    // q = grad(x[k]) - grad(x[k-1])
    q[0] = -1 * grad[0];
    q[1] = -1 * grad[1];

    memcpy(&tmp, out, 2 * sizeof(double));
    //tmp = x[k];
    d_function(tmp[0], tmp[1], grad);

    q[0] += grad[0];
    q[1] += grad[1];

    //BFGS?
    //h[k] = h[k-1] + (1+(qk^T*Hk*qk)*(pk^T*qk)^-1)*((pk*pk^-T)*(pk^T*qk)^-1) - ((pk*qk^T*Hk + Hk*qk*pk^T)*(pk^T*qk)^-1)
    // update est_inv_hess

    if (grad[0] == 0 && grad[1] == 0)
    {
      cout << "Grad exit " << grad[0] << " " << grad[1] << endl;
      break;
    }

    if (current_iteration%100 == 0 || DEBUG)
    {
      cout << "i_HESS: "<<  est_inv_hess[0]<<", "<<est_inv_hess[1]<<", "<<est_inv_hess[2]<<", "<<est_inv_hess[3]<<endl;
      cout << "GRAD: "<<  grad[0]<<", "<<grad[1]<<endl;
    }

    dk[0] = -1 * (est_inv_hess[0]*grad[0] + est_inv_hess[1]*grad[1]);
    dk[1] = -1 * (est_inv_hess[2]*grad[0] + est_inv_hess[3]*grad[1]);

    double t = 1;//goldenSectionSearch(precision, 0, 10, out, dk, function);

    if (current_iteration%100 == 0  || DEBUG)
    {
      cout << "Dk: "<<dk[0]<<" "<<dk[1]<< " t: "<<t<<endl;
    }

    out[0] += t * dk[0];
    out[1] += t * dk[1];

    if (current_iteration%100 == 0 || DEBUG)
    {
      cout <<"Iteration: " << current_iteration <<" Result: (" << out[0] << ", " << out[1] << ")" << endl;
    }

    current_iteration += 1;

    if ((abs(out[0] - tmp[0]) < precision) && (abs(out[1] - tmp[1]) < precision))
    {
      cout << "Precision Exit" << endl;
      break;
    }
  }
  cout << "Proccess finalized for " << name <<". Result: (" << out[0] << ", " << out[1] << ")" << endl;
  cout << "Total iterations: " << current_iteration << endl;
}

int main()
{
  // starting point
  double new_arr [2] = {10, 10};
  // step size
  double gamma = 0.01;
  // precision neede for stop
  double precision = 0.00001;
  // max iterations
  int max_iterations = 100;

  quasiNewton("F1",precision, max_iterations, new_arr, f1, df1, Hf1);

  quasiNewton("F2",precision, max_iterations, new_arr, f2, df2, Hf2);

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}