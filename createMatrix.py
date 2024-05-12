import argparse
import numpy as np


def create_sparse_matrix(filename, p, n):
    # Generate sparse matrix with probability p for each entry
    sparse_matrix = np.random.choice([0, 1], size=(n, n), p=[1-p, p])
    # Convert to float for non-zero entries
    sparse_matrix = sparse_matrix.astype(float)

    # Generate random floating-point numbers for non-zero entries
    non_zero_indices = np.argwhere(sparse_matrix != 0)
    for i, j in non_zero_indices:
        sparse_matrix[i][j] = np.random.rand()

    # Print the dense matrix
    print("Dense Matrix:")
    print(sparse_matrix)
    print()

    # Count non-zero elements
    nnz = np.count_nonzero(sparse_matrix)

    # Write matrix to file
    with open(filename, 'w') as file:
        # Write n and number of nonzeros as the first line
        file.write(f"{n} {nnz}\n")

        # Write each non-zero entry of the matrix
        for i in range(n):
            for j in range(n):
                val = sparse_matrix[i][j]
                if val != 0:
                    file.write(f"{i} {j} {val}\n")


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Create a sparse matrix with probability p')
    parser.add_argument('filename', type=str, help='output filename')
    parser.add_argument(
        'p', type=float, help='probability of non-zero elements')
    parser.add_argument('n', type=int, help='dimension of the matrix')
    args = parser.parse_args()

    # Create sparse matrix and write to file
    create_sparse_matrix(args.filename, args.p, args.n)


if __name__ == "__main__":
    main()
