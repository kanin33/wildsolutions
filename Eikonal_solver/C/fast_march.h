#ifndef FAST_MARCH_H
#define FAST_MARCH_H

// F = 1
void fast_march(Grid *grid, double xmin, double xmax,
                            double ymin, double ymax,
                            double zmin, double zmax,
                            double qmin, double qmax,
                            double boundary(double, double, double, double));

void march_forward(Grid *grid, double xmin, double xmax,
                               double ymin, double ymax,
                               double zmin, double zmax,
                               double qmin, double qmax,
                               double boundary(double, double, double, double),
                               std::vector<ValueIndex> narrowband,
                               long not_known);

bool compare_ValueIndex(ValueIndex v1, ValueIndex v2);

double advance(Grid *grid, int i, int j, int k, int l,
               double dx, double dy, double dz, double dq,
               double xmin, double ymin, double zmin, double qmin);

double sum_squared(std::vector<double> T, int j);

#endif
