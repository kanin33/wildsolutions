#include "Grid.h"
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

Grid::Grid(double standard_value, int shape[4]) {
    this->standard_value = standard_value;
    this->shape = (struct Numbers) {shape[0], shape[1], shape[2], shape[3]};
}

void Grid::setItem(int array[4], double value) {
    unsigned long long ind = array[3] + (array[2] * 1000)
                          + (array[1] * 1000000) + (array[0] * 1000000000);
    elements[ind] = value;
}

double Grid::getItem(int array[4]) {
    unsigned long long ind = array[3] + (array[2] * 1000) + (array[1] * 1000000) + (array[0] * 1000000000);
    return elements[ind];
}

struct Numbers Grid::getShape() {
    return this->shape;
}

void Grid::writeToFile(char* filename, int n, double xmin, double xmax) {
    ofstream file;
    file.open(filename);
    file << n << " " << xmin << " " << xmax << "\n";
    for (const auto& e : elements) {
        file.fill('0');
        file.width(12);
        file << e.first << " " << e.second << "\n";
    }
    file.close();
}
