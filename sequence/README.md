# Wild solutions
Python 2.7 and Matlab programs used to numerically construct an approximation to a wild solution.

### Prerequisites

Dependencies for the Python program are
- Numpy 1.13.3 or newer
- Matplotlib 2.1.0 or newer
- Matlab engine for Python
- ??

Dependencies for the Matlab programs are
- Matlab r2016a or newer
- Chebfun 5.7.0 or newer

# Usage
## Setup
The Matlab files included in *matlab\_programs* must be in the Matlab path on your local machine.
You must have a distance-file with data to check the distance to the boundary of U. This can be created using the Eikonal solver.

## Example
The following example shows how to continue computing a sequence using existing
data from a previous simulation. If there is no data in the given directory,
the file starts a new simulation.

```python
from sequence_examples import *

max_iterations = 10
dir_name       = "test_dir"
dist_filename  = "dist.txt"
r_values       = [0.5 - 0.01*i for i in range(11)]

subsol_values = example_5_6()
mu_0     = subsol_values[2][0]
mu_1     = subsol_values[2][1]
v_tilde  = subsol_values[5]
u_tilde  = subsol_values[6]
C        = subsol_values[7]

t_min   = 1;  t_max = 2
x_1_min = -2; x_1_max = 2

dist_obj = Dist_dU(dist_filename) # dist_obj(x, y, z, q) = dist((x, y, z, q), dU)
P        = [t_min, t_max, mu_0, mu_1, x_1_min, x_1_max]
# create object to save the solution
z_k      = Sequence(C, P, dist_obj,
                    v_tilde, u_tilde,
                    dirname=dir_name)
# read earlier simulations from file
z_k.read_data_from_file()

seq, z_k = compute_gen_ex_5_6(max_iterations,
                              t_min, t_max,
                              x_1_min, x_1_max,
                              dist_filename=dist_filename,
                              r_values=r_values,
                              dirname=dir_name,
                              z_k=z_k)
```


## Authors

**Tale Ulfsby**
