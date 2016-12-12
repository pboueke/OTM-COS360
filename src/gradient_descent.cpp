#include <iostream>
#include <cmath>
#include <cstring>

using namespace std;

void df1(double x, double y, double* out) {
  out[0] = (2 * log (x)) / (x * (1 + log(x) * log (x) + log(y) * log(y)));
  out[1] = (2 * log (y)) / (y * (1 + log(y) * log (y) + log(x) * log(x)));
}

void df2(double x, double y, double* out) {
  out[0] = (4*x*(x*x + y) + 2*x) / ((x*x + y)*(x*x + y) + x*x + 1);
  out[1] = (2*(y + x*x)) / ((y + x*x)*(y + x*x) + x*x + 1);
}

int main()
{
  double old_arr [2] = {0, 0};
  // starting point
  double new_arr [2] = {10, 10};
  // step size
  double gamma = 0.01;
  // precision neede for stop
  double precision = 0.00001;
  // max iterations
  int max_iterations = 100000;
  // counter
  int current_iteration = 0;
  // gradient storage
  double grad [2] = {0, 0};

  cout << "Starting Gradient descent.\nParameters:" << endl;
  cout << "Gamma: " << gamma << endl;
  cout << "Precision: " << precision << endl;
  cout << "Max number of iterations: " << max_iterations << endl;
  cout << "[f1] Starting..." << endl;

  //main loop
  while (
          (    (new_arr[0] - old_arr[0] > precision)
            && (new_arr[1] - old_arr[1] > precision)
          )
          || current_iteration < max_iterations
        )
  {
    //old_arr = new_arr;
    memcpy(&old_arr, &new_arr, 2 * sizeof(double));
    // calculate new derivate
    df1(old_arr[0], old_arr[1], grad);
    new_arr[0] += -gamma * grad[0];
    new_arr[1] += -gamma * grad[1];

    if (current_iteration%10000 == 0)
    {
      cout << "[f1] Iteration: " << current_iteration << endl;
    }
    current_iteration += 1;
  }

  cout << "Proccess finalized for f1. Result: (" << new_arr[0] << ", " << new_arr[1] << ")" << endl;
  cout << "Total iterations for f1: " << current_iteration << endl;
  cout << "[f2] Starting..." << endl;

  old_arr[0] = old_arr[1] = 0;
  // starting point
  new_arr[0] = new_arr[1] = 10;
  //main loop
  current_iteration = 0;
  // gradient storage
  grad[0] = grad[1] = 0;

  while (
          (    (new_arr[0] - old_arr[0] > precision)
            && (new_arr[1] - old_arr[1] > precision)
          )
          || current_iteration < max_iterations
        )
  {
    //old_arr = new_arr;
    memcpy(&old_arr, &new_arr, 2 * sizeof(double));
    // calculate new derivate
    df2(old_arr[0], old_arr[1], grad);
    new_arr[0] += -gamma * grad[0];
    new_arr[1] += -gamma * grad[1];

    if (current_iteration%10000 == 0)
    {
      cout << "[f2] Iteration: " << current_iteration << endl;
    }
    current_iteration += 1;
  }

  cout << "Proccess finalized for f2. Result: (" << new_arr[0] << ", " << new_arr[1] << ")" << endl;
  cout << "Total iterations for f2: " << current_iteration << endl;

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}
