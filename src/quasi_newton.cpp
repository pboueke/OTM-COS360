#include <iostream>
#include <cmath>
#include <cstring>
#include "stdlib.h"
#include "string.h"

#define DEBUG false
#define DEBUG_SPACING 100000

#define DELTAF false
#define LASTX false
#define ITERX true

#define DELTAF_PRECISION 0.000000001
#define PRECISION 0.00001
#define MAX_ITERATIONS 100000
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
  out[0] = (-(4*pow(log(x),2))/(x*x*pow((pow(log(x), 2) + pow(log(y), 2) + 1), 2)) - (2*log(x))/(x*x*(pow(log(x), 2) + pow(log(y), 2) + 1)) + 2/(x*x*(pow(log(x), 2) + pow(log(y), 2) + 1)));
  //xy  -(4 log(x) log(y))/(x y (log^2(x) + log^2(y) + 1)^2)
  out[1] = -1*(4*log(x)*log(y))/(x*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //yx  -(4 log(x) log(y))/(x y (log^2(x) + log^2(y) + 1)^2)
  out[2] = -1*(4*log(x)*log(y))/(x*y*pow(pow(log(x),2) + pow(log(y),2) + 1, 2));
  //yy  -(2 (log^2(x) (log(y) - 1) + log^3(y) + log^2(y) + log(y) - 1))/(y^2 (log^2(x) + log^2(y) + 1)^2)
  out[3] = (-(4*pow(log(y),2))/(y*y*pow((pow(log(x), 2) + pow(log(y), 2) + 1), 2)) - (2*log(y))/(y*y*(pow(log(x), 2) + pow(log(y), 2) + 1)) + 2/(y*y*(pow(log(x), 2) + pow(log(y), 2) + 1)));;
}

void Hf2 (double x, double y, double* out)
{
  //xx    (4 (x^2 - y) + 8 x^2 + 2)/((x^2 - y)^2 + x^2 + 1) - (4 x (x^2 - y) + 2 x)^2/((x^2 - y)^2 + x^2 + 1)^2
  out[0] = (4*(pow(x,2) - y) + 8*pow(x,2) + 2)/(pow(pow(x,2) - y, 2) + pow(x,2) + 1) - pow(4*x*(pow(x,2) - y) + 2*x, 2)/pow(pow(pow(x,2) - y, 2) + pow(x,2) + 1, 2);
  //xy    (2 (4 x (x^2 - y) + 2 x) (x^2 - y))/((x^2 - y)^2 + x^2 + 1)^2 - (4 x)/((x^2 - y)^2 + x^2 + 1)
  out[1] = (2*(4*x*(pow(x,2) - y) + 2*x)*(pow(x,2) - y))/(pow(pow(pow(x,2) - y, 2) + pow(x,2) + 1, 2)) - (4*x)/(pow(pow(x,2) - y, 2) + pow(x,2) + 1);
  //yx    (2 (4 x (x^2 - y) + 2 x) (x^2 - y))/((x^2 - y)^2 + x^2 + 1)^2 - (4 x)/((x^2 - y)^2 + x^2 + 1)
  out[2] = (2*(4*x*(pow(x,2) - y) + 2*x)*(pow(x,2) - y))/(pow(pow(pow(x,2) - y, 2) + pow(x,2) + 1, 2)) - (4*x)/(pow(pow(x,2) - y, 2) + pow(x,2) + 1);
  //yy   -(2 (x^4 - x^2 (2 y + 1) + y^2 - 1))/(x^4 + x^2 (1 - 2 y) + y^2 + 1)^2
  out[3] = 2/(pow(pow(x,2) - y, 2) + pow(x,2) + 1) - (4*pow(pow(x,2) - y, 2))/pow(pow(pow(x,2) - y, 2) + pow(x,2) + 1, 2);
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

void update_est_inv_hess(double *hess, double *p, double *q) {
  double old_hess[4];
  for (int i = 0; i < 4; i++) {
    //cout << "OLDHESS: " << hess[i] << endl;
    old_hess[i] = hess[i];
  }
  double common_divider = p[0]*q[0] + p[1]*q[1];

  double first_factor = 1;
  double first_ratio = q[0] * (q[0]*old_hess[0] + q[1]*old_hess[2]);
  first_ratio += q[1] * (q[0]*old_hess[1] + q[1]*old_hess[3]);
  first_ratio /= common_divider;
  first_factor += first_ratio;

  double first_matrix[4];
  first_matrix[0] = p[0]*p[0];
  first_matrix[1] = p[0]*p[1];
  first_matrix[2] = p[1]*p[0];
  first_matrix[3] = p[1]*p[1];
  for (int i = 0; i < 4; i++) {
    first_matrix[i] *= first_ratio;
  }

  double aux_matrix_a[4];
  aux_matrix_a[0] = old_hess[0]*p[0]*q[0] + old_hess[2]*p[0]*q[1];
  aux_matrix_a[1] = old_hess[1]*p[0]*q[0] + old_hess[3]*p[0]*q[1];
  aux_matrix_a[2] = old_hess[0]*p[1]*q[0] + old_hess[2]*p[1]*q[1];
  aux_matrix_a[3] = old_hess[1]*p[1]*q[0] + old_hess[3]*p[1]*q[1];

  double aux_matrix_b[4];
  aux_matrix_b[0] = old_hess[0]*q[0]*p[0] + old_hess[1]*q[1]*p[0];
  aux_matrix_b[1] = old_hess[0]*q[0]*p[1] + old_hess[1]*q[1]*p[1];
  aux_matrix_b[2] = old_hess[2]*q[0]*p[0] + old_hess[3]*q[1]*p[0];
  aux_matrix_b[3] = old_hess[2]*q[0]*p[1] + old_hess[3]*q[1]*p[1];

  for (int i = 0; i < 4; i++) {
    aux_matrix_a[i] += aux_matrix_b[i];
    aux_matrix_a[i] /= common_divider;
  }

  for (int i = 0; i < 4; i++) {
    hess[i] = old_hess[i] + first_matrix[i] - aux_matrix_a[i];
    //cout << "HESS: " << hess[i] << endl;
  }
  // exit(0);
}

void quasiNewton (string name, double precision, int max_iterations, double *out, double (*function) (double, double), void (*d_function)(double, double, double*), void (*h_function)(double, double, double*))
{
  int current_iteration = 0;
  double grad [2] = {0, 0};
  double est_inv_hess [4] = {1, 0, 0, 1};//I
  double tmp_hess [4];
  double dk [2] = {0, 0};
  double tmp [2] = {0, 0};
  double p [2] = {0, 0};
  double q [2] = {0, 0};

  cout << "Starting Quasi-Newtonn method.\nParameters:" << endl;
  cout << "Precision: " << precision << endl;
  cout << "Max number of iterations: " << max_iterations << endl;
  cout << "Starting... at " << out[0] << " " << out[1] << endl;
  cout << "---------------------------------------------" << endl;

  //init grad
  d_function(out[0], out[1], grad);

  //main loop
  while (!(grad[0] == 0 && grad[1] == 0))
  {
    memcpy(&tmp, out, 2 * sizeof(double));
    //tmp = x[k];

    if (current_iteration%DEBUG_SPACING == 0 || DEBUG)
    {
      cout << "i_HESS: "<<  est_inv_hess[0]<<", "<<est_inv_hess[1]<<", "<<est_inv_hess[2]<<", "<<est_inv_hess[3]<<endl;
      cout << "GRAD: "<<  grad[0]<<", "<<grad[1]<<endl;
      cout << "p: "<<p[0] << " "<<p[1]<<" q: "<<q[0] << " "<<q[1]<<endl;
    }

    dk[0] = -1 * (est_inv_hess[0]*grad[0] + est_inv_hess[1]*grad[1]);
    dk[1] = -1 * (est_inv_hess[2]*grad[0] + est_inv_hess[3]*grad[1]);

    // golden section search
		double t;
		if (name == "F1")
		{
			// we cant use GSS on f1 because it is not unimodal
			t = 0.01;
		}
		else
		{
			t = goldenSectionSearch(precision, 0, 10, out, dk, function);
		}

    if (current_iteration%DEBUG_SPACING == 0  || DEBUG)
    {
      cout << "Dk: "<<dk[0]<<" "<<dk[1]<< " t: "<<t<<endl;
    }

    out[0] += t * dk[0];
    out[1] += t * dk[1];

    if ((current_iteration > max_iterations) && ITERX)
		{
			break;
		}
		if ((abs(function(out[0], out[1]) - function(tmp[0], tmp[1])) < DELTAF_PRECISION) && DELTAF)
		{
			break;
		}
		if (((abs(out[0] - tmp[0]) < precision) && (abs(out[1] - tmp[1]) < precision)) && LASTX)
		{
			break;
		}

    q[0] = -1*grad[0];
    q[1] = -1*grad[1];
    d_function(out[0], out[1], grad);
    q[0] += grad[0];
    q[0] += grad[1];

    p[0] = out[0] - tmp[0];
    p[1] = out[1] - tmp[1];

    update_est_inv_hess(est_inv_hess, p, q);
    //invert2dmatrix(est_inv_hess, tmp_hess);
    //memcpy(&est_inv_hess, tmp_hess, 2 * sizeof(double));

    if (current_iteration%DEBUG_SPACING == 0 || DEBUG)
    {
      cout << "new GRAD: "<<  grad[0]<<", "<<grad[1]<<endl;
      cout << "p: "<<p[0] << " "<<p[1]<<" q: "<<q[0] << " "<<q[1]<<endl;
      cout <<"Iteration: " << current_iteration <<" Result: (" << out[0] << ", " << out[1] << ") Value: " << function(out[0], out[1]) << endl;
      cout << "---------------------------------------------" << endl;
    }

    current_iteration += 1;
  }
  cout << "Proccess finalized for " << name <<". Result: (" << out[0] << ", " << out[1] << ") Value: " << function(out[0], out[1]) << endl;
  cout << "Total iterations: " << current_iteration << endl;
  cout << "\n\n=========================================\n\n";
}

int main()
{
  // starting point
  double new_arr [2];
  // step size
  // precision neede for stop
  double precision = PRECISION;
  // max iterations
  int max_iterations = MAX_ITERATIONS;

  // new_arr[0] = 0.9;
  // new_arr[1] = 0.9;
  // quasiNewton("F1",precision, max_iterations, new_arr, f1, df1, Hf1);

  new_arr[0] = 0;
  new_arr[1] = 0;
  quasiNewton("F2",precision, max_iterations, new_arr, f2, df2, Hf2);

  cout << "Exiting..." << endl;
  //cin.ignore();

  return 0;
}
