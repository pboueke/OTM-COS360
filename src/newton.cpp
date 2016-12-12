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

void df1(double x, double y, double* out)
{
  out[0] = (2 * log (x)) / (x * (1 + log(x) * log (x) + log(y) * log(y)));
  out[1] = (2 * log (y)) / (y * (1 + log(y) * log (y) + log(x) * log(x)));
}

void df2(double x, double y, double* out)
{
  out[0] = (4*x*(x*x + y) + 2*x) / ((x*x + y)*(x*x + y) + x*x + 1);
  out[1] = (2*(y + x*x)) / ((y + x*x)*(y + x*x) + x*x + 1);
}

void Hf1 (double x, double y, double* out)
{
  //xx  -(2 (log(x) (log^2(y) + 1) + log^3(x) + log^2(x) - log^2(y) - 1))/(x^2 (log^2(x) + log^2(y) + 1)^2)
  out[0] = -1*(2*(log(x)*(pow(log(y),2) + 1) + pow(log(x),3) + pow(log(x),2) - pow(log(y),2) - 1))/(x*x*pow((pow(log(x),2) + pow(log(y),2) + 1), 2));
  //xy  -(4 log(x) log(y))/(x y (log^2(x) + log^2(y) + 1)^2)
  out[1] = -1*(4*log(x)*log(y))/(x*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //yx  -(4 log(x) log(y))/(x y (log^2(x) + log^2(y) + 1)^2)
  out[2] = -1*(4*log(x)*log(y))/(x*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //yy  -(2 (log^2(x) (log(y) - 1) + log^3(y) + log^2(y) + log(y) - 1))/(y^2 (log^2(x) + log^2(y) + 1)^2)
  out[3] = -1*(2*(pow(log(x),2)*(log(y) - 1) + pow(log(y),3) + pow(log(y),2) + log(y) - 1))/(y*y*pow(pow(log(x)2,) + pow(log(y),2) + 1, 2));
}

void Hf2 (double x, double y, double* out)
{
  //xx    (4 (x^2 + y) + 8 x^2 + 2)/((x^2 + y)^2 + x^2 + 1) - (4 x (x^2 + y) + 2 x)^2/((x^2 + y)^2 + x^2 + 1)^2
  out[0] = (4*(x*x + y) + 8*x*x + 2)/(pow(x*x + y,2) + x*x + 1) - (4*x (x*x + y) + 2*x)^2/pow(pow(x*x + y,2) + x*x + 1, 2);
  //xy   -(4 x (x^4 + 2 x^2 y + y^2 + y - 1))/(x^4 + x^2 (2 y + 1) + y^2 + 1)^2
  out[1] = -1 * (4*x (pow(x,4) + 2*x*x y + y*y + y - 1))/pow(pow(x,4) + x*x (2*y + 1) + y*2 + 1, 2);
  //yx   -(4 x (x^4 + 2 x^2 y + y^2 + y - 1))/(x^4 + x^2 (2 y + 1) + y^2 + 1)^2
  out[2] =  -1 *(4*x (pow(x, 4) + 2 x*x*y + y*y + y - 1))/pow(pow(x,4) + x*x*(2*y + 1) + y*y + 1, 2);
  //yy   -(2 (x^4 + x^2 (2 y - 1) + y^2 - 1))/(x^4 + x^2 (2 y + 1) + y^2 + 1)^2
  out[3] = -1 * (2*(pow(x,4) + x*x*(2*y - 1) + y*y - 1))/pow(x,4) + x*x (2*y + 1) + y*y + 1, 2);
}

void newton (string name, double precision, int max_iterations, double *out, void (*d_function)(double, double, double*), double (*function)(double, double))
{
  int current_iteration = 0;
  double grad [2] = {0, 0};
  double tmp [2] = {0, 0};
  out[0] = out[1] = 10;

  cout << "Starting Gradient descent.\nParameters:" << endl;
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
    double f_value = function(out[0], out[1]);

    cout << "out[0]: " <<out[0]<<" func: " << f_value << " der[0]: " << grad[0] << " ratio: " << f_value / grad[0] << " final: " << out[0] - f_value / grad[0] << endl;

    out[0] = tmp[0] - f_value / grad[0];
    out[1] = tmp[1] - f_value / grad[1];
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

  newton("F1",precision, 3, new_arr, df1, f1);

  //newton("F2",precision, max_iterations, new_arr, df2, f2);

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}
