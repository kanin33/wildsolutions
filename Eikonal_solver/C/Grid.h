#ifndef GRID_H
#define GRID_H

#include "sparsepp/spp.h"

struct Numbers {
    int a;
    int b;
    int c;
    int d;
};

class Grid {
    private:
        spp::sparse_hash_map<unsigned long long, double> elements;
        double standard_value;
        struct Numbers shape;
    public:
        Grid(double standard_value, int array[4]);
        void setItem(int array[4], double value);
        double getItem(int array[4]);
        struct Numbers getShape();
        void writeToFile(char* filename, int n, double xmin, double xmax);
};
#endif
