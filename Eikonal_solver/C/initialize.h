#include <vector>
#ifndef INITIALIZE_H
#define INITIALIZE_H

struct ValueIndex {
    double value;
    int i;
    int j;
    int k;
    int l;
};

struct Initialized {
    Grid grid;
    std::vector<ValueIndex> narrowband;
    long not_known;
};

Initialized initialize(Grid *grid, double xmin, double xmax,
                                  double ymin, double ymax,
                                  double zmin, double zmax,
                                  double qmin, double qmax,
                                  double boundary(double, double, double, double));

double dist_p0_boundary(double p0[4], double p1[4],
                        double boundary(double, double, double, double));


bool add_boundary(int i, int j, int k, int l, 
                  double dx, double dy, double dz, double dq,
                  double xmin, double ymin, double zmin, double qmin, Grid *grid,
                  double boundary(double, double, double, double),
                  std::vector<ValueIndex> *narrowband);
#endif
