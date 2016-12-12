#include <iostream>
#include <cmath>
#include <cstring>
#include "string.h"

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
  out[0] = (2*(2*x*x*x - 2*x*y + x))/(pow(x*x - y, 2) + x*x + 1)
  // y  -(2 (x^2 - y))/((x^2 - y)^2 + x^2 + 1)
  out[1] = -1 * (2*(x*x - y))/(pow(x*x - y, 2) + x*x + 1)
}

double goldenSectionSearch (double precision, double in_a, double in_b, double *x, double *dx, double (*f) (double, double))
{
  const double gr = ((sqrt(5)-1) / 2);

  // initialize bounds
  double a = in_a;
  double c = (in_a + in_b) / 2; //center
  double b = in_b;

  while (f(x[0]+c*dx[0], x[1]+c*dx[1]) < f(x[0]+c*dx[0], x[1]+c*dx[1]))
  {
    a = c;
    c = b;
    b = b + (c - a);
  }

  double d = a + (1 - gr) * (b - a);
  double e = a + gr * (b - a);

  double aux_a = a;
  double aux_b = b;

  while (b - a > precision)
  {
    if (f(x[0]+d*dx[0], x[1]+d*dx[1]) < f(x[0]+e*dx[0], x[1]+e*dx[1]))
    {
      b = e;
      e = d;
      d = a + (1 - gr) * (b - a);
    }
    else
    {
      a = d;
      d = e;
      e = a + gr * (b - a);
    }

    if ((aux_a == a) && (aux_b == b))
    {
      break;
    }
    aux_a = a;
    aux_b = b;

  }
  return (a + b) / 2;
}

void gradient (string name, double precision, double gamma, int max_iterations, double *out, void (*d_function)(double, double, double*), double (*function) (double, double))
{
  int current_iteration = 0;
  double grad [2] = {0, 0};
  double tmp [2] = {0, 0};
  double last_out [2] = {0, 0};
  out[0] = out[1] = 10;

  cout << "Starting Gradient descent.\nParameters:" << endl;
  cout << "Gamma: " << gamma << endl;
  cout << "Precision: " << precision << endl;
  cout << "Max number of iterations: " << max_iterations << endl;
  cout << "Starting..." << endl;
  //main loop
  while (current_iteration < max_iterations)
  {
    if (current_iteration%5000 == 0)
    {
      cout << "Iteration: " << current_iteration << endl;
    }
    // calculate new derivate
    d_function(out[0], out[1], grad);
    tmp[0] += -gamma * grad[0];
    tmp[1] += -gamma * grad[1];

    // golden section search
    double t = goldenSectionSearch(precision, 0, 10, out, tmp, function);

    last_out[0] = out[0];
    last_out[1] = out[1];

    out[0] = out[0] + t * tmp[0];
    out[1] = out[1] + t * tmp[1];

    current_iteration += 1;

    if ((abs(out[0] - last_out[0]) < precision) && (abs(out[1] - last_out[1]) < precision))
    {
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
  // precision needed for stop
  double precision = 0.0000001;
  // max iterations
  int max_iterations = 1000000;

  gradient("F1",precision, gamma, max_iterations, new_arr, df1, f1);

  // is F2 unimodal?
  gradient("F2",precision, gamma, max_iterations, new_arr, df2, f2);

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}
