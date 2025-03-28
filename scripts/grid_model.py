# Make a random single-epoch model that follows a square grid.
#   Usage: python grid_model.py <rows> <columns>
import random
import sys
import numpy as np

RATE_LB = 1.0e-5
RATE_UB = 0.01


def print_matrix_row(row, is_first, is_last, indent=2):
    if is_first:
        sys.stdout.write((" " * indent) + "- [ ")
    else:
        sys.stdout.write((" " * indent) + (" " * 4))
    sys.stdout.write("[")
    for i, value in enumerate(row):
        sys.stdout.write(str(value))
        if i != len(row) - 1:
            sys.stdout.write(", ")
    sys.stdout.write("]")
    if is_last:
        sys.stdout.write(" ]")
    else:
        sys.stdout.write(",")
    sys.stdout.write("\n")


def print_matrix(matrix):
    for i, row in enumerate(matrix):
        is_first = i == 0
        is_last = i == len(matrix) - 1
        print_matrix_row(row, is_first, is_last)


HEADER = """ploidy: 2
epochTimeSplit: null
populationConversion: null
populationCoordinates: null"""

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: grid_model.py <rows> <columns>")
        exit(1)
    rows = int(sys.argv[1])
    cols = int(sys.argv[2])
    demes = rows * cols

    print(HEADER)
    coal_rates = []
    for d in range(demes):
        row = [0 for _ in range(demes)]
        row[d] = d + 1
        coal_rates.append(row)
    print("coalescence:")
    print("  matrices:")
    print_matrix(coal_rates)
    print("  parameters:")
    for d in range(demes):
        rate = random.uniform(RATE_LB, RATE_UB)
        print(f"  - ground_truth: {rate}")
        print(f"    lb: {RATE_LB}")
        print(f"    ub: {RATE_UB}")

    parameters = 0
    mig_rates = np.zeros((demes, demes), dtype=int)

    def row_col_to_d(i, j):
        return (i * cols) + j

    for i in range(rows):
        for j in range(cols):
            src = row_col_to_d(i, j)
            left = j - 1
            right = j + 1
            up = i - 1
            down = i + 1
            if left >= 0:
                tgt = row_col_to_d(i, left)
                parameters += 1
                mig_rates[src, tgt] = parameters
            if right < cols:
                tgt = row_col_to_d(i, right)
                parameters += 1
                mig_rates[src, tgt] = parameters
            if up >= 0:
                tgt = row_col_to_d(up, j)
                parameters += 1
                mig_rates[src, tgt] = parameters
            if down < rows:
                tgt = row_col_to_d(down, j)
                parameters += 1
                mig_rates[src, tgt] = parameters
    print("migration:")
    print("  matrices:")
    print_matrix(mig_rates)
    print("  parameters:")
    for _ in range(parameters):
        rate = random.uniform(RATE_LB, RATE_UB)
        print(f"  - ground_truth: {rate}")
        print(f"    lb: {RATE_LB}")
        print(f"    ub: {RATE_UB}")
