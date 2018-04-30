#include <iostream>
#include "Grid.h"
#include "initialize.h"
#include "fast_march.h"

using namespace std;

// boundary for the domain used in the thesis
double C = 3.507715078;
double boundary(double x, double y, double z, double q) {
    // compute both eigenvalues
    double tmp1 = 1/2. * (-(pow(x, 2) + pow(y,2)) + C);
    double tmp2 = 1/2. * sqrt(pow(x,4) + 2*pow(x,2)*pow(y,2) - 4*pow(x,2)*z - 8*x*y*q
                        + pow(y,4) + 4*pow(y,2)*z + 4*pow(z,2) + 4*pow(q,2));
    double lambda_1 = tmp1 + tmp2;
    double lambda_2 = tmp1 - tmp2;
    if(lambda_1 > 0 && lambda_2 > 0) {
        return -lambda_1*lambda_2;
    }
    return 0;
}

// sphere boundary
//double boundary(double x, double y, double z, double q) {
//    return pow(x,2)+ pow(y,2) + pow(z,2) + pow(q,2) - 0.5*0.5;
//}

int main(int argc, char* argv[]) {

    int n = 50;
    int number[4] = {n, n, n, n};
    double xmin = -sqrt(3.507715078);
    double xmax = sqrt(3.507715078);

    Grid *g = new Grid(0, number);
    fast_march(g, xmin, xmax,
                  xmin, xmax,
                  xmin, xmax,
                  xmin, xmax,
                  boundary);

    // write solution to file
    g->writeToFile("distance_50p.txt", n, xmin, xmax);
    return 0;
}
