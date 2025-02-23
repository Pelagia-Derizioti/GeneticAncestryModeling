from typing import List, Dict
import numpy as np
import json
import os
import logging


# Sum of a list of integers
def sum_list_of_int(a: List[int]) -> int:
    s = 0
    for i in range(len(a)):
        s += a[i]
    return s

# Sum of a list of floats
def sum_list_of_float(a: List[float]) -> float:
    s = 0.0
    for i in range(len(a)):
        s += a[i]
    return s

# Sum of values in a specific column (vertical sum) in a 2D list
def sum_vertical(a: List[List[int]], j: int) -> int:
    s = 0
    for i in range(len(a)):
        s += a[i][j]
    return s

# Mean of a list of integers
def mean_list_of_int(a: List[int]) -> float:
    s = sum_list_of_int(a)
    return float(s) / float(len(a))

# Mean of a list of floats
def mean_list_of_float(a: List[float]) -> float:
    s = sum_list_of_float(a)
    return float(s) / float(len(a))

# Mean ignoring zeros in the list
def mean_non_zero(a: List[int]) -> float:
    s = sum_list_of_int(a)
    n = 0
    for i in range(len(a)):
        if a[i] != 0:
            n += 1
    if n==0:
        return None
    return float(s) / float(n)

# Create a matrix (n x m) filled with zeros
def new_matrix(n: int, m: int) -> List[List[int]]:
    matrix = [[] for i in range(n)]
    for i in range(n):
        for j in range(m):
            matrix[i].append(0)
    return matrix

# Find the minimum in a list of integers
def min_list_of_int(a: List[int]) -> int:
    return min(a)

# Find the maximum in a list of integers
def max_list_of_int(a: List[int]) -> int:
    return max(a)

# Read integers from a file
def read_file(file_path: str) -> List[int]:
    try:
        with open(file_path, 'r') as file:
            data = file.read()
    except Exception as e:
        print(f"Error reading file: {e}")
        raise ValueError("Error reading file")

    integers = []
    parts = data.split()

    for part in parts:
        try:
            num = int(part)
        except Exception as e:
            print(f"Error inside Auxiliary_functions.py in function read_file(file_path) converting string to integer: {part}")
            print(f"Error: {e}")
            continue

        integers.append(num)

    return integers

# Print dictionary as formatted JSON
def print_json(data: Dict[str, object]) -> None:
    json_data = json.dumps(data, sort_keys=True, indent=2)
    print(json_data)

    