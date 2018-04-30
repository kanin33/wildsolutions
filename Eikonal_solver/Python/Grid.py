import numpy as np

class Grid:
    """
    Class for sparse storage of numpy_float64 values
    """
    def __init__(self, standard_value, shape):
        self.elements = {}
        self.standard_value = standard_value
        self.shape = shape

    def __setitem__(self, index, value):
        self.elements[index] = value

    def __getitem__(self, index):
        if index in self.elements:
            return self.elements[index]
        return self.standard_value

    def to_array(self):
        array = np.zeros(self.shape)
        array[:] = self.standard_value
        for index in self.elements:
            array[index] = self.elements[index]
        return array

if __name__ == '__main__':
    grid = Grid(0)
    grid[0, 0, 0, 0] = 1
    print grid[0, 0, 0, 0]
    print grid[1, 1, 1, 1]
