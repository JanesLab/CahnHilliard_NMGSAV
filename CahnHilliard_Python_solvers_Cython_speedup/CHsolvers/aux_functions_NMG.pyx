import numpy as np
import tracemalloc
import time


def dmatrix(nrows, ncols):
    """
    Create an empty matrix of size nrows x ncols
    :param nrows: Number of rows
    :param ncols: Number of columns
    :return: Matrix of size nrows x ncols
    """
    return np.empty((nrows, ncols), dtype=float)


def print_data(filename, a):
    """
    Write data to file in space-delimited format
    :param filename: Name of file to write data to
    :param a: Data that is written to file
    :return: None
    """
    np.savetxt(filename, a, fmt="%16.15f", delimiter=" ")


def calculate_mass(phi, h2, nx, ny):
    ave_mass = np.sum(phi) / (h2 * nx * ny)
    return ave_mass


def f(phi):
    fphi = (1 / 4) * ((phi ** 2) - 1) ** 2
    return fphi


def calculate_discrete_energy(phi, h2, epsilon2):
    a = h2 * np.sum(f(phi))
    sum_i = np.sum((phi[2:-1, :] - phi[1:-2, :]) ** 2)

    b = (epsilon2 / 2) * sum_i
    sum_j = np.sum((phi[:, 2:-1] - phi[:, 1:-2]) ** 2)

    c = (epsilon2 / 2) * sum_j
    E = a + b + c
    return E


# def calculate_discrete_norm_energy(phi, phi0, h2, epsilon2):
#     E0 = calculate_discrete_energy(phi0, h2, epsilon2)
#     E = calculate_discrete_energy(phi, h2, epsilon2)
#     return E / E0


def time_and_mem(func):
    """
    Decorator to measure execution time and memory usage of a function.
    :param func: Function to be decorated
    :return: results of the function with time and memory measurement in a dictionary called result_dict. The results from func will be have key values of results1, results2, etc.
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        # Start tracing memory allocations
        tracemalloc.start()
        results = func(*args, **kwargs)
        # End measuring time
        end_time = time.time()
        # Get current memory usage
        current, peak = tracemalloc.get_traced_memory()
        # Stop tracing memory allocations
        tracemalloc.stop()

        # Ensure results is always a tuple (for consistency)
        if not isinstance(results, tuple):
            results = (results,)

        # Calculate execution time
        execution_time = end_time - start_time
        # print(f"Execution time: {execution_time:.4f} seconds")
        # print(f"Current memory usage: {current / 1024 / 1024:.2f} MB")
        # print(f"Peak memory usage: {peak / 1024 / 1024:.2f} MB")
        result_dict = {f"results{i+1}": result for i,
                       result in enumerate(results)}
        # Add computation time
        result_dict["computation_time_Sec"] = execution_time
        # Add memory usage
        result_dict["memory_usage_MB"] = peak / 1024 / 1024
        return result_dict
    return wrapper
