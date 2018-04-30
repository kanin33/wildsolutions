#include <iostream>
#include <algorithm> // for heap
#include "Grid.h"
#include "initialize.h"
#include "fast_march.h"
using namespace std;

void fast_march(Grid *grid, double xmin, double xmax,
                            double ymin, double ymax,
                            double zmin, double zmax,
                            double qmin, double qmax,
                            double boundary(double, double, double, double)) {
    // initialize
    Initialized init = initialize(grid, xmin, xmax,
                                        ymin, ymax,
                                        zmin, zmax,
                                        qmin, qmax,
                                        boundary);
    cout << "Done with initialize" << endl;

    // March forward
    march_forward(grid, xmin, xmax,
                        ymin, ymax,
                        zmin, zmax,
                        qmin, qmax,
                        boundary,
                        init.narrowband,
                        init.not_known);
    cout << "Done with fast march" << endl;
}

void march_forward(Grid *grid, double xmin, double xmax,
                               double ymin, double ymax,
                               double zmin, double zmax,
                               double qmin, double qmax,
                               double boundary(double, double, double, double),
                               std::vector<ValueIndex> narrowband,
                               long not_known) {
    int n1 = grid->getShape().a;
    int n2 = grid->getShape().b;
    int n3 = grid->getShape().c;
    int n4 = grid->getShape().d;

    double dx = (xmax - xmin) / (n1 - 1);
    double dy = (ymax - ymin) / (n2 - 1);
    double dz = (zmax - zmin) / (n3 - 1);
    double dq = (qmax - qmin) / (n4 - 1);

    make_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);

    ValueIndex A;
    int i; int j; int k; int l;
    int index[4];
    double gip1;
    double value;
    while(not_known > 0) {
        // Step 1
        A = narrowband.front();
        pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
        // remove A from vector
        narrowband.pop_back();

        // Step 2
        i = A.i; j = A.j;
        k = A.k; l = A.l;
        grid->setItem((int[4]) {i,j,k,l}, -A.value);

        // Step 3
        gip1 = grid->getItem((int[4]) {i+1, j, k, l});
        // If Trial of Far, update value
        if(gip1 > 0) {
            value = advance(grid, i+1, j, k, l, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            grid->setItem((int[4]) {i+1, j, k, l}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i+1, j, k, l};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i-1, j, k, l});
        if(gip1 > 0) {
            value = advance(grid, i-1, j, k, l, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            grid->setItem((int[4]) {i-1, j, k, l}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i-1, j, k, l};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i, j+1, k, l});
        if(gip1 > 0) {
            value = advance(grid, i, j+1, k, l, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            grid->setItem((int[4]) {i, j+1, k, l}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i, j+1, k, l};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i, j-1, k, l});
        if(gip1 > 0) {
            value = advance(grid, i, j-1, k, l, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            grid->setItem((int[4]) {i, j-1, k, l}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i, j-1, k, l};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i, j, k+1, l});
        if(gip1 > 0) {
            value = advance(grid, i, j, k+1, l, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            //grid->setItem(index, value);
            grid->setItem((int[4]) {i, j, k+1, l}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i, j, k+1, l};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i, j, k-1, l});
        if(gip1 > 0) {
            value = advance(grid, i, j, k-1, l, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            //grid->setItem(index, value);
            grid->setItem((int[4]) {i, j, k-1, l}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i, j, k-1, l};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i, j, k, l+1});
        if(gip1 > 0) {
            value = advance(grid, i, j, k, l+1, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            //grid->setItem(index, value);
            grid->setItem((int[4]) {i, j, k, l+1}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i, j, k, l+1};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        gip1 = grid->getItem((int[4]) {i, j, k, l-1});
        if(gip1 > 0) {
            value = advance(grid, i, j, k, l-1, dx, dy, dz, dq,
                            xmin, ymin, zmin, qmin);
            grid->setItem((int[4]) {i, j, k, l-1}, value);
            // If Far, add to narrow band
            if(gip1 == INFINITY) {
                ValueIndex v = {value, i, j, k, l-1};
                // add to vector
                narrowband.push_back(v);
                // push to heap
                pop_heap (narrowband.begin(), narrowband.end(), compare_ValueIndex);
            }
        }

        not_known--;
    }

}

bool compare_ValueIndex(ValueIndex v1, ValueIndex v2) {
    return v1.value > v2.value;
}

double advance(Grid *grid, int i, int j, int k, int l,
               double dx, double dy, double dz, double dq,
               double xmin, double ymin, double zmin, double qmin) {
    std::vector<double> T;
    int index1[4] = {i+1, j, k, l}; int index2[4] = {i-1, j, k, l};
    double T_H = min(abs(grid->getItem(index1)), abs(grid->getItem(index2)));
    if(T_H != INFINITY) {
        T.push_back(T_H);
    }
    index1[0] = i; index1[1] = j+1; index2[0] = i; index2[1] = j-1;
    double T_V = min(abs(grid->getItem(index1)), abs(grid->getItem(index2)));
    if(T_V != INFINITY) {
        T.push_back(T_V);
    }
    index1[1] = j; index1[2] = k+1; index2[1] = j; index2[2] = k-1;
    double T_Z = min(abs(grid->getItem(index1)), abs(grid->getItem(index2)));
    if(T_Z != INFINITY) {
        T.push_back(T_Z);
    }
    index1[2] = k; index1[3] = l+1; index2[2] = k; index2[3] = l-1;
    double T_Q = min(abs(grid->getItem(index1)), abs(grid->getItem(index2)));
    if(T_Q != INFINITY) {
        T.push_back(T_Q);
    }

    std::sort (T.begin(), T.end());
    // NB: n here is j in python code
    int n = T.size();
    double dx_square = pow(dx, 2);
    while(sum_squared(T, n) >= dx_square) {
        n--;
    }
    double v1 = 0;
    double v2 = 0;
    for(int i = 0; i < n; i++) {
        v1 += T[i];
        v2 += pow(T[i], 2);
    }
    double new_value = v1 + sqrt(pow(v1, 2) - n*(v2 - dx_square));
    new_value *= 1./n;
    return new_value;
}

double sum_squared(std::vector<double> T, int j) {
    double sum = 0;
    for(int i = 0; i < j-1; i++) {
        sum += pow(T[j-1] - T[i], 2);
    }
    return sum;
}

