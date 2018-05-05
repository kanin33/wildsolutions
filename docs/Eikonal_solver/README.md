# Eikonal solver

Python 2.7 and C++ programs to solve the Eikonal equation ``|grad(u(x))|=1/f(x),  x in Omega `` when ``Omega`` is a 4 dimensional bounded domain.

### Prerequisites

Dependencies for the Python program are
- Numpy 1.13.3 or newer
- Matplotlib 2.1.0 or newer

Dependencies for the C++ program are
- Sparseepp ([GitHub repo](https://github.com/greg7mdp/sparsepp))
- Clang++ for compiling

### Python
To compute the solution to the Eikonal equation when ``Omega`` is a sphere and ``f`` is constantly equal to 1, simply type:
```python
solution_grid = Grid(standard_value=0, shape=[100]*4)
def boundary_function(x, y, z, q):
    return x**2 + y**2 + z**2 + q**2 - 1
fast_march(grid=solution_grid,
           xmin=-2, xmax=2,
           ymin=-2, ymax=2,
           zmin=-2, zmax=2,
           qmin=-2, qmax=2,
           boundary=boundary_function,
           F=lambda x, y, z, q: 1)
print solution_grid[50, 50, 50, 50]
```
The solution will be stored in the Grid-object. The above code will print the solution in the mid point of the grid, which corresponds to the point ``(0, 0, 0, 0)``.

The boundary function should return a negative number if the given point is inside the domain, and a positive number if it is outside.

Note that ``Omega`` should be completely contained in the box defined by the minimum and maximum values.

### C++
##### Compiling
The program should be compiled with Clang++ using the included Makefile. From the Terminal run
```
>> make
```
when in the C++ directory. The output file is named *tale*.

Note that the Sparseepp library must be contained in the same directory as compiling from.

##### Running
To compute the solution to the Eikonal equation when ``Omega`` is a sphere and ``f`` is constantly equal to 1, this should be written in the main file:
```C++
int n = 50;
int number[4] = {n, n, n, n};
double xmin = -2;
double xmax = 2;

Grid *g = new Grid(0, number);

fast_march(g, xmin, xmax,
              xmin, xmax,
              xmin, xmax,
              xmin, xmax,
              boundary);
```
where `boundary` is a function of four variables, given by
```C++
double boundary(double x, double y, double z, double q) {
    return pow(x,2)+ pow(y,2) + pow(z,2) + pow(q,2) - 1;
}
```

Note that in the C++ program the ``f`` function is always 1.

## Authors

**Tale Ulfsby**

## Acknowledgments

Thanks to Halvard Kværna for help with the C++ program.