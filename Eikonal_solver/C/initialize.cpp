#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Grid.h"
#include "initialize.h"
using namespace std;

Initialized initialize(Grid *grid, double xmin, double xmax,
                                   double ymin, double ymax,
                                   double zmin, double zmax,
                                   double qmin, double qmax,
                                   double boundary(double, double, double, double)
                                   ) {
    std::vector<ValueIndex> narrowband = {};
    int n1 = grid->getShape().a;
    int n2 = grid->getShape().b;
    int n3 = grid->getShape().c;
    int n4 = grid->getShape().d;

    double dx = (xmax - xmin) / (n1 - 1);
    double dy = (ymax - ymin) / (n2 - 1);
    double dz = (zmax - zmin) / (n3 - 1);
    double dq = (qmax - qmin) / (n4 - 1);

    // initialize inner points
    long unknown = n1*n2*n3*n4;
    int l;
    for(int i = 0; i < n1; i++) {
    for(int j = 0; j < n2; j++) {
    for(int k = 0; k < n3; k++) {
        for(l = 0; l < n4; l++) {
            if (boundary(i*dx + xmin, j*dy + ymin,
                         k*dz + zmin, l*dq + qmin) < 0) {
                if(!add_boundary(i, j, k, l, dx, dy, dz, dq,
                                 xmin, ymin, zmin, qmin, grid, boundary, &narrowband)) {
                    break;
                }
            } else {
                unknown--;
            }
        }
        for(int m = n4-1; m > l; m--) {
            if (boundary(i*dx + xmin, j*dy + ymin,
                         k*dz + zmin, m*dq + qmin) < 0) {
                if(!add_boundary(i, j, k, m, dx, dy, dz, dq,
                                 xmin, ymin, zmin, qmin, grid, boundary, &narrowband)) {
                    break;
                }
            } else {
                unknown--;
            }
        }
    }
    }
    }
    return (struct Initialized) {*grid, narrowband, unknown};
}

double dist_p0_boundary(double p0[4], double p1[4],
                        double boundary(double, double, double, double)) {
    double t = boundary(p0[0], p0[1], p0[2], p0[3]) /
              (boundary(p1[0], p1[1], p1[2], p1[3]) -
               boundary(p0[0], p0[1], p0[2], p0[3]));
    return sqrt(pow(t*(p0[0] - p1[0]), 2) + pow(t*(p0[1] - p1[1]), 2)
              + pow(t*(p0[2] - p1[2]), 2) + pow(t*(p0[3] - p1[3]), 2));
}

bool add_boundary(int i, int j, int k, int l, 
                  double dx, double dy, double dz, double dq,
                  double xmin, double ymin, double zmin, double qmin, Grid *grid,
                  double boundary(double, double, double, double),
                  std::vector<ValueIndex> *narrowband) {
    double s1 = INFINITY; double s2 = INFINITY; double s3 = INFINITY;
    double s4 = INFINITY; double s5 = INFINITY; double s6 = INFINITY;
    double s7 = INFINITY; double s8 = INFINITY;

    double p0[4] = {0, 0, 0, 0};
    double p1[4] = {0, 0, 0, 0};

    double zero = 0;
    p0[0] = xmin + i*dx;     p0[1] = ymin + j*dy;
    p0[2] = zmin + k*dz;     p0[3] = qmin + l*dq;
    p1[0] = xmin + (i+1)*dx; p1[1] = p0[1];
    p1[2] = p0[2];           p1[3] = p0[3];
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s1 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i+1, j, k, l}, zero);
    }
    p1[0] = xmin + (i-1)*dx; p1[1] = p0[1];
    p1[2] = p0[2];           p1[3] = p0[3];
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s5 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i-1, j, k, l}, zero);
    }
    p1[0] = p0[0];           p1[1] = ymin + (j+1)*dy;
    p1[2] = p0[2];           p1[3] = p0[3];
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s2 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i, j+1, k, l}, zero);
    }
    p1[0] = p0[0];           p1[1] = ymin + (j-1)*dy;
    p1[2] = p0[2];           p1[3] = p0[3];
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s6 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i, j-1, k, l}, zero);
    }
    p1[0] = p0[0];           p1[1] = p0[1];
    p1[2] = zmin + (k+1)*dz; p1[3] = p0[3];
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s3 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i, j, k+1, l}, zero);
    }
    p1[0] = p0[0];           p1[1] = p0[1];
    p1[2] = zmin + (k-1)*dz; p1[3] = p0[3];
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s7 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i, j, k-1, l}, zero);
    }
    p1[0] = p0[0];           p1[1] = p0[1];
    p1[2] = p0[2];           p1[3] = qmin + (l+1)*dq;
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s4 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i, j, k, l+1}, zero);
    }
    p1[0] = p0[0];           p1[1] = p0[1];
    p1[2] = p0[2];           p1[3] = qmin + (l-1)*dq;
    if (boundary(p1[0], p1[1], p1[2], p1[3]) >= 0) {
        s8 = dist_p0_boundary(p0, p1, boundary);
        grid->setItem((int[4]) {i, j, k, l-1}, zero);
    }
    double t1 = dx*min(s1, s5);
    double t2 = dy*min(s2, s6);
    double t3 = dz*min(s3, s7);
    double t4 = dq*min(s4, s8);

    double value = sqrt(1/(1/pow(t1, 2)+ 1/pow(t2, 2) + 1/pow(t3, 2) + 1/pow(t4, 2)));
    if(value == INFINITY) {
        return false;
    }
    grid->setItem((int[4]) {i, j, k, l}, value);
    ValueIndex v = {value, i, j, k, l};
    narrowband->push_back(v);
    return true;
}
