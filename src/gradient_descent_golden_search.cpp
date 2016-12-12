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
  return log(1 + x * x + (x * x + y) * (x * x + y));
}

void df1(double x, double y, double* out) {
  out[0] = (2 * log (x)) / (x * (1 + log(x) * log (x) + log(y) * log(y)));
  out[1] = (2 * log (y)) / (y * (1 + log(y) * log (y) + log(x) * log(x)));
}

void df2(double x, double y, double* out) {
  out[0] = (4*x*(x*x + y) + 2*x) / ((x*x + y)*(x*x + y) + x*x + 1);
  out[1] = (2*(y + x*x)) / ((y + x*x)*(y + x*x) + x*x + 1);
}

double goldenSectionSearch (double precision, double in_a, double in_b, double *x, double *dx, double (*f) (double, double))
{
  const double gr = ((sqrt(5)-1) / 2);

  // initialize bounds
  double a = in_a;
  double c = (in_a + in_b) / 2; //center
  double b = in_b;

  //cout << "GSS in (a, b, c): (" << a << ", "<<b<<", "<<c<<")   "<< endl;

  while (f(x[0]+c*dx[0], x[1]+c*dx[1]) < f(x[0]+c*dx[0], x[1]+c*dx[1]))
  {
    cout << f(x[0]+c*dx[0], x[1]+c*dx[1]) << " < " << f(x[0]+c*dx[0], x[1]+c*dx[1]) << endl;
    a = c;
    c = b;
    b = b + (c - a);
  }

  double d = a + (1 - gr) * (b - a);
  double e = a + gr * (b - a);

  //cout << "GSS reduced (a, b, c, d, e): ("<<a<<", "<<b<<", "<<c<<", "<<d<<", "<<e<<")"<< endl;
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
  //cout << "GSS final (a, b, c, d, e): ("<<a<<", "<<b<<", "<<c<<", "<<d<<", "<<e<<")"<< endl;
  return (a + b) / 2;
}

void gradient (string name, double precision, double gamma, int max_iterations, double *out, void (*d_function)(double, double, double*), double (*function) (double, double))
{
  int current_iteration = 0;
  double grad [2] = {0, 0};
  double tmp [2] = {0, 0};
  out[0] = out[1] = 10;

  cout << "Starting Gradient descent.\nParameters:" << endl;
  cout << "Gamma: " << gamma << endl;
  cout << "Precision: " << precision << endl;
  cout << "Max number of iterations: " << max_iterations << endl;
  cout << "Starting..." << endl;
  //main loop
  while (current_iteration < max_iterations)
  {
    if (current_iteration%1000 == 0)
    {
      //cout << "Iteration: " << current_iteration << endl;
    }
    //tmp = out;
    memcpy(&tmp, out, 2 * sizeof(double));
    // calculate new derivate
    d_function(tmp[0], tmp[1], grad);
    out[0] += -gamma * grad[0];
    out[1] += -gamma * grad[1];

    //cout << "Iteration: " << current_iteration << " Golden Search" << endl;
    // golden section search
    double t = goldenSectionSearch(precision, 0, 10, tmp, out, function);

    out[0] = out[0] + t * out[0];
    out[1] = out[1] + t * out[1];

    current_iteration += 1;

    if ((abs(out[0] - tmp[0]) < precision) && (abs(out[1] - tmp[1]) < precision))
    {
      break;
    }
  }
  cout << "Proccess finalized for " << name <<". Result: (" << out[0] << ", " << out[1] << ")" << endl;
  cout << "Total iterations: " << current_iteration << endl;
  //delete [] grad;
  //delete [] tmp;
}

int main()
{
  // starting point
  double new_arr [2] = {10, 10};
  // step size
  double gamma = 0.01;
  // precision neede for stop
  double precision = 0.0000001;
  // max iterations
  int max_iterations = 100000;

  gradient("F1",precision, gamma, max_iterations, new_arr, df1, f1);

  gradient("F2",precision, gamma, max_iterations, new_arr, df2, f2);

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}
