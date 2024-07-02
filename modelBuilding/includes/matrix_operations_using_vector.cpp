// matrix_operations_using_vector.cpp
// 注意这里还有不少函数没有作isNaN或者isEmpty的判断，先放掉。20231122
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "matrix_operations_using_vector.h"

const double eps = 1e-10;

std::vector<std::vector<double>> madd(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    std::vector<std::vector<double>> matrix3;
    try {
        if (m1 == m2 && n1 == n2) {
            int m3 = m1;
            int n3 = n1;
            matrix3.resize(m3);
            for (int i=0; i<m3; i++) {
                matrix3[i].resize(n3);
            }
            for (int i=0; i<m3; i++) {
                for (int j=0; j<n3; j++) {
                    matrix3[i][j] = matrix1[i][j] + matrix2[i][j];
                }
            }
        } else {
            std::string dimensionError = "DimensionError: Number of Cols and Rows of Matrix1 and Matrix2 do not match! Current Value: m1 = " 
                + std::to_string(m1) + ", m2 = " + std::to_string(m2) + "; n1 = " + std::to_string(n1) + ", n2 = " + std::to_string(n2) + "\n";
            throw std::invalid_argument(dimensionError);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return matrix3;
}

std::vector<std::vector<double>> msubtract(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    std::vector<std::vector<double>> matrix3;
    try {
        if (m1 == m2 && n1 == n2) {
            int m3 = m1;
            int n3 = n1;
            matrix3.resize(m3);
            for (int i=0; i<m3; i++) {
                matrix3[i].resize(n3);
            }
            for (int i=0; i<m3; i++) {
                for (int j=0; j<n3; j++) {
                    matrix3[i][j] = matrix1[i][j] - matrix2[i][j]; //only this sign was changed from madd
                }
            }
        } else {
            std::string dimensionError = "DimensionError: Number of Cols and Rows of Matrix1 and Matrix2 do not match! Current Value: m1 = " 
                + std::to_string(n1) + ", m2 = " + std::to_string(m2) + "; n1 = " + std::to_string(n1) + ", n2 = " + std::to_string(n2) + "\n";
            throw std::invalid_argument(dimensionError);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return matrix3;
}

std::vector<std::vector<double>> mscalarmult(double scalar, std::vector<std::vector<double>>& matrix1) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    std::vector<std::vector<double>> matrix3;
    int m3 = m1, n3 = n1;
    matrix3.resize(m3);
    for (int i=0; i<m3; i++) {
        matrix3[i].resize(n3);
    }
    for (int i=0; i<m3; i++) {
        for (int j=0; j<n3; j++) {
            matrix3[i][j] = matrix1[i][j]*1.0*scalar;
        }
    }
    return matrix3;
}

std::vector<std::vector<double>> mmult(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    std::vector<std::vector<double>> matrix3;
    try {
        if (n1 == m2) {
            int m3 = m1;
            int n3 = n2;
            matrix3.resize(m3);
            for (int i=0; i<m3; i++) {
                matrix3[i].resize(n3);
            }
            double tmp = 0.0;
            for (int i=0; i<m3; i++) {
                for (int j=0; j<n3; j++) {
                    tmp = 0.0;
                    for (int k=0; k<n1; k++) {
                        tmp += matrix1[i][k] * matrix2[k][j];
                    }
                    matrix3[i][j] = tmp;
                }
            }

        } else {
            std::string dimensionError = "DimensionError: Number of Cols of Matrix1 and number of Rows of Matrix2 do not match! Current Value: n1 = " 
                + std::to_string(n1) + ", m2 = " + std::to_string(m2) + "\n";
            throw std::invalid_argument(dimensionError);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return matrix3;
}

std::vector<std::vector<double>> mtranspose(std::vector<std::vector<double>>& matrix1) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    std::vector<std::vector<double>> matrix3;
    int m3 = n1;
    int n3 = m1;
    matrix3.resize(m3);
    for (int i=0; i<m3; i++) { // @ 20231126 bug caught: was m1
        matrix3[i].resize(n3);
    }
    for (int i=0; i<m3; i++) {
        for (int j=0; j<n3; j++) {
            matrix3[i][j] = matrix1[j][i];
        }
    }
    return matrix3;
}

std::vector<std::vector<double>> asmechelonize(std::vector<std::vector<double>>& matrix1) { //asm stands for augmented square matrix
    // for non-square matrix, the try-exception logic will be a little different.
    // since this mmult is only set for solving determinant and inverse matrix,
    // we temparorily will not write code for non-square matrix.
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    std::vector<std::vector<double>> matrix3;
    try {
        matrix3.resize(m1); //@Keep matrix3 untouched before going into try
        for (int i=0; i<m1; i++) { //augmented matrix
            matrix3[i].resize(n1*2);
            for (int j=0; j<n1; j++) {
                matrix3[i][j] = matrix1[i][j];
            }
            matrix3[i][n1+i] = 1.0;
        }
        if (m1 != n1) {
            std::string nonSquareMatrixError = "NonSquareMatrixError: This is a non-square matrix.\n";
            throw std::invalid_argument(nonSquareMatrixError);
        }
        for (int j=0; j<n1; j++) {
            if ( fabs(matrix3[j][j] - 0.0) < eps ) { //if principal diagonal element is zero, we need to exchange two rows. 
                bool principal_diagonal_element_zero_flag = true;
                for (int i=j+1; i<m1; i++) {
                    if ( fabs(matrix3[i][j] - 0.0) >= eps ) { //@20240221 Note that fabs (for double) replaced abs.
                        std::vector<double> tmpRow = matrix3[i];
                        for (int k=0; k<n1; k++) {
                            matrix3[i][k] = - matrix3[j][k]; //exchanging two rows will change the sign of determinant
                            matrix3[j][k] = tmpRow[k]; 
                        }
                        principal_diagonal_element_zero_flag = false;
                        break;
                    }
                }
                if (principal_diagonal_element_zero_flag == true) {
                    std::string singularMatrixError = "SingularMatrixError: This is a singlular matrix.\n";
                    throw std::invalid_argument(singularMatrixError);
                }
            }
            for (int i=j+1; i<m1; i++) { // notice the upper and lower limits of i and j
                double scale = matrix3[i][j] *1.0 / matrix3[j][j]; 
                for (int k=j; k<n1*2; k++) {
                    matrix3[i][k] -= scale *1.0 * matrix3[j][k];
                }
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return matrix3;
}

double det(std::vector<std::vector<double>>& matrix1) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    std::vector<std::vector<double>> matrix3 = asmechelonize(matrix1);
    double determinantValue;
    if (matrix3.size()) {
        determinantValue = 1.0;
        for (int j=0; j<n1; j++) {
            determinantValue *= 1.0*matrix3[j][j];
        }
    } else {
        std::cerr << "Error: During calculation of the determinant, error(s) occur in augmented square matrix echelonization.\n";
    }
    return determinantValue;
}

std::vector<std::vector<double>> minverse(std::vector<std::vector<double>>& matrix1) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    std::vector<std::vector<double>> matrix3 = asmechelonize(matrix1);
    std::vector<std::vector<double>> matrix4;
    if (matrix3.size()) {
        for (int i=0; i<m1; i++) {
            double principal_diagonal_element_before_division = matrix3[i][i]; 
            //@ we need to extract this value before division, otherwise matrix3[i][i] changes after divided by itself
            for (int j=i; j<2*n1; j++) {
                matrix3[i][j] /= (1.0 * principal_diagonal_element_before_division);
            }
        }
        for (int j=n1-1; j>=1; j--) {
            for (int i=j-1; i>=0; i--) { // notice the upper and lower limits of i and j
                double scale = matrix3[i][j] *1.0 / matrix3[j][j]; 
                for (int k=j; k<n1*2; k++) {
                    matrix3[i][k] -= scale *1.0 * matrix3[j][k];
                }
            }
        }
        int m4 = m1;
        int n4 = n1;
        matrix4.resize(m4);
        for (int i=0; i<m4; i++) {
            matrix4[i].resize(n4);
            for (int j=0; j<n4; j++) {
                matrix4[i][j] = matrix3[i][j+n4];
            }
        }
    } else {
        std::cerr << "Error: During calculation of inverse matrix, error(s) occur in augmented square matrix echelonization.\n";
    }
    return matrix4;
}

std::vector<std::vector<double>> mblock(std::vector<std::vector<double>>& matrix1, int start_row, int start_col, int block_rows, int block_cols) {
    // Follow the custom of linear algebra, that first element of a matrix is a_ij, where i = j = 1.
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    std::vector<std::vector<double>> matrix3;
    try {
        if (start_row - 1 + block_rows <= m1 && start_col - 1 + block_cols <= n1) {
            matrix3.resize(block_rows);
            for (int i=0; i<block_rows; i++) {
                matrix3[i].resize(block_cols);
            }
            for (int i=0; i<block_rows; i++) {
                for (int j=0; j<block_cols; j++) {
                    matrix3[i][j] = matrix1[start_row-1+i][start_col-1+j];
                }
            }
        } else {
            std::string dimensionError = "DimensionError: Current Block Will Exceed the Index of Matrix! Current Value: m1 = " 
                + std::to_string(m1) + ", n1 = " + std::to_string(n1) + "; start_row = " + std::to_string(start_row)
                + ", start_col = " + std::to_string(start_col) + ", block_rows = " + std::to_string(block_rows)
                + ", block_cols = " + std::to_string(block_cols) + "\n"; 
            throw std::invalid_argument(dimensionError);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }  
    return matrix3; 
}

double dot(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    double value3;
    try {
        value3 = 0.0;
        if (m1 == 1 && m2 == 1 && n1 == n2) {
            for (int k=0; k<n1; k++) {
                value3 += matrix1[0][k] * matrix2[0][k];
            }
        } else if (m1 == 1 && n2 == 1 && n1 == m2) {
            for (int k=0; k<n1; k++) {
                value3 += matrix1[0][k] * matrix2[k][0];
            }
        } else if (n1 == 1 && m2 == 1 && m1 == n2) {
            for (int k=0; k<m1; k++) {
                value3 += matrix1[k][0] * matrix2[0][k];
            }
        } else if (n1 == 1 && n2 == 1 && m1 == m2) {
            for (int k=0; k<m1; k++) {
                value3 += matrix1[k][0] * matrix2[k][0];
            }
        } else { //Either of the two matrices is not a vector
            std::string dimensionError = "DimensionError: Either of the two matrices is not a vector, Or two vectors have different length! Current Value: m1 = " 
                + std::to_string(m1) + ", n1 = " + std::to_string(n1) + "; m2 = " + std::to_string(m2) + ", n2 = " + std::to_string(n2) + "\n";
            throw std::invalid_argument(dimensionError);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return value3;
}

double norm(std::vector<std::vector<double>>& matrix1) {
    // Here we still presume that this matrix is a vector.
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    double value3;
    try {
        if (m1 == 1 || n1 == 1) {
            value3 = std::sqrt(dot(matrix1, matrix1));
        } else { 
            std::string dimensionError = "DimensionError: The matrix is not a vector! Current Value: m1 = " 
                + std::to_string(m1) + ", n1 = " + std::to_string(n1) + "\n";
            throw std::invalid_argument(dimensionError);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return value3;
}

std::vector<std::vector<double>> cross(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) {
    // Here we also allow row vectors.
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    std::vector<std::vector<double>> matrix5;
    try {
        if ((m1 == 1 || n1 == 1) && (m2 == 1 || n2 == 1)) {
            std::vector<std::vector<double>> matrix3;
            std::vector<std::vector<double>> matrix4;
            if (n1 == 1) {
                matrix3 = matrix1;
            } else {
                matrix3 = mtranspose(matrix1);
            }
            if (n2 == 1) {
                matrix4 = matrix2;
            } else {
                matrix4 = mtranspose(matrix2);
            }
            if (matrix3.size() == 3 && matrix4.size() == 3) {
                matrix5 = {{0.0}, {0.0}, {0.0}};
                matrix5[0][0] = matrix3[1][0] * matrix4[2][0] - matrix3[2][0] * matrix4[1][0]; // a2b3-a3b2
                matrix5[1][0] = matrix3[2][0] * matrix4[0][0] - matrix3[0][0] * matrix4[2][0]; // a2b3-a3b2
                matrix5[2][0] = matrix3[0][0] * matrix4[1][0] - matrix3[1][0] * matrix4[0][0]; // a2b3-a3b2
            } else {
                std::string dimensionError = "DimensionError: At least one of the vectors is not at size 3! Current Value: vector1.size() = " 
                    + std::to_string(matrix3.size()) + ", vector2.size() = " + std::to_string(matrix4.size()) + "\n";
                throw std::invalid_argument(dimensionError);
            }
        } else { 
            std::string dimensionError = "DimensionError: At least one of the matrices is not a vector! Current Value: m1 = " 
                + std::to_string(m1) + ", n1 = " + std::to_string(n1) + "; m2 = " + std::to_string(m2) + ", n2 = " + std::to_string(n2) + "\n";
            throw std::invalid_argument(dimensionError);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
    return matrix5;
}

void mprint(std::vector<std::vector<double>>& matrix1) {
    if (matrix1.size() == 0) {
        std::cerr << "Error: Cannot print the matrix called, which is currently empty. Please trace back errors.\n";
    } else {
        for (int i=0; i<matrix1.size(); i++) {
            for (int j=0; j<matrix1[0].size(); j++) {
                std::cout << " " << matrix1[i][j];
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}
