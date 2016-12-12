#include <iostream>
#include <cmath>
#include <cstring>
#include "string.h"

using namespace std;

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

void gradient (string name, double precision, double gamma, int max_iterations, double *out, void (*d_function)(double, double, double*))
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
    //tmp = out;
    memcpy(&tmp, out, 2 * sizeof(double));
    // calculate new derivate
    d_function(tmp[0], tmp[1], grad);
    out[0] += -gamma * grad[0];
    out[1] += -gamma * grad[1];
    current_iteration += 1;

    if ((abs(out[0] - tmp[0]) < precision) && (abs(out[1] - tmp[1]) < precision))
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
  // precision neede for stop
  double precision = 0.00001;
  // max iterations
  int max_iterations = 100000;

  gradient("F1",precision, gamma, max_iterations, new_arr, df1);

  gradient("F2",precision, gamma, max_iterations, new_arr, df2);

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}
