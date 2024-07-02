// matrix_operations_using_vector.h
#ifndef MATRIX_OPERATIONS_USING_VECTOR_H
#define MATRIX_OPERATIONS_USING_VECTOR_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

extern const double eps;

std::vector<std::vector<double>> madd(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2);
std::vector<std::vector<double>> msubtract(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2);
std::vector<std::vector<double>> mscalarmult(double scalar, std::vector<std::vector<double>>& matrix1);
std::vector<std::vector<double>> mmult(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2);
std::vector<std::vector<double>> mtranspose(std::vector<std::vector<double>>& matrix1);
std::vector<std::vector<double>> asmechelonize(std::vector<std::vector<double>>& matrix1); //asm stands for augmented square matrix
double det(std::vector<std::vector<double>>& matrix1); 
std::vector<std::vector<double>> minverse(std::vector<std::vector<double>>& matrix1);
std::vector<std::vector<double>> mblock(std::vector<std::vector<double>>& matrix1, int start_row, int start_col, int block_rows, int block_cols);
double dot(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2);
double norm(std::vector<std::vector<double>>& matrix1);
std::vector<std::vector<double>> cross(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2);

void mprint(std::vector<std::vector<double>>& matrix1);

#endif // MATRIX_OPERATIONS_USING_VECTOR_H