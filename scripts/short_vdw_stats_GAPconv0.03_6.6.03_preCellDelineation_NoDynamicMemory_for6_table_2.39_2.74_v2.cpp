//20240221 This version of "short_vdw_stats_preCellDelineation_for6_table_2.39_2.74.cpp" becomes a TongRenWen of the ones (Matrices and StructuralInformation) included in descriptor2Model
//20240227 v2 solved problems of uninitialized pair_graph for last atom - 256; removed completedly resize in inner loop
// matrix_operations_using_vector.cpp
// 注意这里还有不少函数没有作isNaN或者isEmpty的判断，先放掉。20231122
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>

const double eps = 1e-10;

std::vector<std::vector<double>> madd(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) {
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    std::vector<std::vector<double> > matrix3;
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
                    if ( fabs(matrix3[i][j] - 0.0) >= eps ) { //@20240221 fabs (for double) replaced abs (for int)
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


#define pi 3.1415926535897932

struct Pair {
	int atom1;
	int atom2;
	std::vector<int> move_direction_for_second_atom;
	double nearest_distance;
	std::vector<double> average_coordinate;
};

class StructuralInformation {
public:
	std::vector<std::vector<double>> lattice_parameter;
	std::vector<std::vector<double>> coordinates;
	std::vector<std::string> elements;
	std::vector<int> n_elements;
	bool isDirect; //不完善 
	std::vector<Pair> pair; 
    int pair_size;
    std::vector<std::vector<int>> pair_graph;
    int total_atoms_n;

    void initialize_pair_graph();
    std::vector<int> determine_move_direction_for_second_atom(int atom1, int atom2);
	int find_pair_index(int atom1, int atom2, int defDir1, int defDir2, int defDir3);
    double dist(int atom1, int atom2, int defDir1, int defDir2, int defDir3);
    std::vector<double> nearest_average(int atom1, int atom2);
	double angle(int atom1, int atom2, int atom3); //只解决最近原子对间的夹角
	double dihedral_angle(int atom1, int atom2, int atom3, int atom4); //只解决最近原子对间的二面角
};


// StructuralInformation.cpp

const double INF = 2147483647;

void StructuralInformation::initialize_pair_graph() {
    total_atoms_n = 0;
    for (int i=0; i<n_elements.size(); i++) {
        total_atoms_n += n_elements[i];
    }
    pair_graph.resize(total_atoms_n+1);
    for (int i=0; i<=total_atoms_n; i++) { //@20240227 3hours debug: this should be <=
        pair_graph[i].resize(total_atoms_n+1);
        for (int j=0; j<=total_atoms_n; j++) { //@20240227 3hours debug: this should be <=
            pair_graph[i][j] = -1;
        }
    }
    //reserve the pair seem no effects
    //pair.reserve(total_atoms_n*(total_atoms_n-1)/2);
    pair.resize(total_atoms_n*(total_atoms_n-1)/2);
    pair_size = 0;
    return;
}

std::vector<int> StructuralInformation::determine_move_direction_for_second_atom(int atom1, int atom2) {
	double cartesian_difference[27][3] = {0.0};
	double min_distance = INF;
	double cartesian_distance[27] = {0.0};
	std::vector<int> tmp_move_direction_for_second_atom(3);
	int cnt=0; 
	for (int di=-1; di<=1; di++) {
		for (int dj=-1; dj<=1; dj++) {
			for (int dk=-1; dk<=1; dk++) {
				cartesian_difference[cnt][0] = coordinates[atom1-1][0] - coordinates[atom2-1][0] - ( di*lattice_parameter[0][0] + dj*lattice_parameter[1][0] + dk*lattice_parameter[2][0] ); //20231231 revise sign
				cartesian_difference[cnt][1] = coordinates[atom1-1][1] - coordinates[atom2-1][1] - ( di*lattice_parameter[0][1] + dj*lattice_parameter[1][1] + dk*lattice_parameter[2][1] );
				cartesian_difference[cnt][2] = coordinates[atom1-1][2] - coordinates[atom2-1][2] - ( di*lattice_parameter[0][2] + dj*lattice_parameter[1][2] + dk*lattice_parameter[2][2] );
				for (int j=0; j<3; j++) {
					cartesian_distance[cnt] += cartesian_difference[cnt][j] * cartesian_difference[cnt][j];
				}
				cartesian_distance[cnt] = sqrt(cartesian_distance[cnt]);
				if (cartesian_distance[cnt] < min_distance) {
					min_distance = cartesian_distance[cnt];
					tmp_move_direction_for_second_atom[0] = di;
					tmp_move_direction_for_second_atom[1] = dj;
					tmp_move_direction_for_second_atom[2] = dk;
				} 
				cnt++;	
			}
		}
	}
	Pair tmp_pair;
	tmp_pair.atom1 = atom1;
	tmp_pair.atom2 = atom2;
	tmp_pair.move_direction_for_second_atom = tmp_move_direction_for_second_atom;
	
	pair[pair_size] = tmp_pair;
    pair_size++;
	return tmp_move_direction_for_second_atom;
}

int StructuralInformation::find_pair_index(int atom1, int atom2, int defDir1 = -2, int defDir2 = -2, int defDir3 = -2 ) {
	if (defDir1 == -2 ) { // No pre-defined direction was inputed
		bool pair_written_flag = false;
		int pair_index = -1;
		for (int i=0; i<pair_size; i++) {
			if (pair[i].atom1 == atom1 && pair[i].atom2 == atom2) {
				pair_written_flag = true;
				pair_index = i;
			}
		}
		if (pair_written_flag == false) {
			determine_move_direction_for_second_atom(atom1, atom2);
			for (int i=0; i<pair_size; i++) {
				if (pair[i].atom1 == atom1 && pair[i].atom2 == atom2) {
					pair_written_flag = true;
					pair_index = i;
				}
			}
		}
		return pair_index;
	} else { // Pre-defined direction was inputed: We first decided to calculate twice, but this will cause more problems. So we optimize this time. 
		bool pair_written_flag = false;
		int pair_index = pair_graph[atom1][atom2];
        if (pair_index != -1) {
            pair_written_flag = true;
        }
		if (pair_written_flag == false) {
            pair_index = pair_size;
            if (defDir2 == -2 || defDir3 == -2) std::cout << "Error: Pre-defined Move Directions Are Not Consistent!\n";
            else {
                std::vector<int> tmp_move_direction_for_second_atom = {defDir1, defDir2, defDir3};
                Pair tmp_pair;
                tmp_pair.atom1 = atom1;
                tmp_pair.atom2 = atom2;
                tmp_pair.move_direction_for_second_atom = tmp_move_direction_for_second_atom;
                
                pair[pair_size] = tmp_pair;
                pair_size++;
                pair_graph[atom1][atom2] = pair_index;
                pair_graph[atom2][atom1] = pair_index;
            }
		}
        
		return pair_index;
	}
}


double StructuralInformation::dist(int atom1, int atom2, int defDir1 = -2, int defDir2 = -2, int defDir3 = -2 ) {// 
	int pair_index;
	if (defDir1 == -2 ) {
		pair_index = find_pair_index(atom1, atom2);
	} else {
		pair_index = find_pair_index(atom1, atom2, defDir1, defDir2, defDir3);
	}
	double cartesian_difference[3] = {0.0};
	cartesian_difference[0] = coordinates[atom1-1][0] - coordinates[atom2-1][0] - ( pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][0] ); //20231231 revise sign
	cartesian_difference[1] = coordinates[atom1-1][1] - coordinates[atom2-1][1] - ( pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][1] );
	cartesian_difference[2] = coordinates[atom1-1][2] - coordinates[atom2-1][2] - ( pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][2] );
	double cartesian_distance = 0.0;
	for (int j=0; j<3; j++) {
		cartesian_distance += cartesian_difference[j] * cartesian_difference[j];
	}
	cartesian_distance = std::sqrt(cartesian_distance);
	pair[pair_index].nearest_distance = cartesian_distance;
	return cartesian_distance;
}

std::vector<double> StructuralInformation::nearest_average(int atom1, int atom2) { //不完善，顺序交换时未判断
	int pair_index = find_pair_index(atom1, atom2);
	double atom2_nearest_cartesian_coordinates[3] = {0.0};
	atom2_nearest_cartesian_coordinates[0] = coordinates[atom2-1][0] + pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][0]; //20231231 revise sign
	atom2_nearest_cartesian_coordinates[1] = coordinates[atom2-1][1] + pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	atom2_nearest_cartesian_coordinates[2] = coordinates[atom2-1][2] + pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	std::vector<double> tmp_average_coordinates(3);
	for (int j=0; j<3; j++) {
		tmp_average_coordinates[j] = 0.5 * (coordinates[atom1-1][j] + atom2_nearest_cartesian_coordinates[j]);
	}
	return tmp_average_coordinates;
}

double StructuralInformation::angle(int atom1, int atom2, int atom3) {// 
	int pair_index[2] = {-1, -1};
	pair_index[0] = find_pair_index(atom2, atom1);
	pair_index[1] = find_pair_index(atom2, atom3);
	
	std::vector<std::vector<double>> two_vectors(2, std::vector<double>(3, 0.0));
	int i = 0;
	two_vectors[i][0] = coordinates[atom2-1][0] - coordinates[atom1-1][0] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0] ); //20231231 revise sign
	two_vectors[i][1] = coordinates[atom2-1][1] - coordinates[atom1-1][1] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1] );
	two_vectors[i][2] = coordinates[atom2-1][2] - coordinates[atom1-1][2] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2] );
	i = 1;
	two_vectors[i][0] = coordinates[atom2-1][0] - coordinates[atom3-1][0] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0] );
	two_vectors[i][1] = coordinates[atom2-1][1] - coordinates[atom3-1][1] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1] );
	two_vectors[i][2] = coordinates[atom2-1][2] - coordinates[atom3-1][2] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2] );
	
	std::vector<std::vector<double>> u1 = mblock(two_vectors, 1, 1, 1, 3);
	std::vector<std::vector<double>> u2 = mblock(two_vectors, 2, 1, 1, 3);

	double norm1 = norm(u1);
	double norm2 = norm(u2);

	double this_angle = 0.0;
	try {
		if (norm1 > eps && norm2 > eps) {
			this_angle = 180.0/pi*acos(dot(u1, u2)*1.0/(norm1 * norm2));
		} else {
			std::string zeroVectorError = "ZeroVectorError: At least one of the vectors is zero vector. Possible overlapping of the atoms chosen.\n";
			throw std::invalid_argument(zeroVectorError);
		}
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }  
	return this_angle;
}

double StructuralInformation::dihedral_angle(int atom1, int atom2, int atom3, int atom4) {// 
	int pair_index[3] = {-1, -1, -1};
	pair_index[0] = find_pair_index(atom1, atom2);
	pair_index[1] = find_pair_index(atom2, atom3);
	pair_index[2] = find_pair_index(atom3, atom4);
	
	std::vector<std::vector<double>> three_vectors(3, std::vector<double>(3, 0.0));
	int i = 0;
	three_vectors[i][0] = coordinates[atom1-1][0] - coordinates[atom2-1][0] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0] ); //20231231 revise sign
	three_vectors[i][1] = coordinates[atom1-1][1] - coordinates[atom2-1][1] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1] );
	three_vectors[i][2] = coordinates[atom1-1][2] - coordinates[atom2-1][2] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2] );
	i = 1;
	three_vectors[i][0] = coordinates[atom2-1][0] - coordinates[atom3-1][0] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0] );
	three_vectors[i][1] = coordinates[atom2-1][1] - coordinates[atom3-1][1] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1] );
	three_vectors[i][2] = coordinates[atom2-1][2] - coordinates[atom3-1][2] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2] );
	i = 2;
	three_vectors[i][0] = coordinates[atom3-1][0] - coordinates[atom4-1][0] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0] );
	three_vectors[i][1] = coordinates[atom3-1][1] - coordinates[atom4-1][1] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1] );
	three_vectors[i][2] = coordinates[atom3-1][2] - coordinates[atom4-1][2] - ( pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2] );
	
	std::vector<std::vector<double>> u1 = mblock(three_vectors, 1, 1, 1, 3);
	std::vector<std::vector<double>> u2 = mblock(three_vectors, 2, 1, 1, 3);
	std::vector<std::vector<double>> u3 = mblock(three_vectors, 3, 1, 1, 3);
	
	std::vector<std::vector<double>> cross_product1 = cross(u1, u2);
	std::vector<std::vector<double>> cross_product2 = cross(u2, u3);

	double norm1 = norm(cross_product1);
	double norm2 = norm(cross_product2);

	double this_dihedral_angle = 0.0;
	try {
		if (norm1 > eps && norm2 > eps) {
			this_dihedral_angle = 180.0/pi*acos(dot(cross_product1, cross_product2)*1.0/(norm1 * norm2));
		} else {
			std::string zeroVectorError = "ZeroVectorError: At least one of the normal vectors of the dihedral plane is zero vector. Possible collinearity of the atoms chosen.\n";
			throw std::invalid_argument(zeroVectorError);
		}
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }  
	return this_dihedral_angle;
}



//#define pi 3.1415926535897932
#define left -1
#define right 1
#define DESCRIPTOR_CNT 2615
#define PCLU_GROUP 17 // 17 coordinates for consideration of a P_Cluster

#define SHORT_VDW_LOWER_LIMIT 2.39
#define SHORT_VDW_UPPER_LIMIT 2.74
#define SHORT_VDW_SAFETY_CONST 0.1

struct node {
	int atom_i=-1;
	int atom_j=-1;
	double atom_pair_dist=-1.0;
};

std::vector<std::vector<double>> lattice_parameters = {
	{18.7999992371, 0.0, 0.0},
    {-4.2539077811, 12.4957698264, 0.0},
    {-3.0286706819, -0.6510739326, 22.6894962298}
};

std::vector<std::vector<double>> lattice_parameters_inverse; 

int main () {
	
	lattice_parameters_inverse = minverse(lattice_parameters);
	double cell_volume = det(lattice_parameters);
	std::vector<std::vector<double>> tmp_vector_a = mblock(lattice_parameters, 1,1,1,3);
	std::vector<std::vector<double>> tmp_vector_b = mblock(lattice_parameters, 2,1,1,3);
	std::vector<std::vector<double>> tmp_vector_c = mblock(lattice_parameters, 3,1,1,3);
	std::vector<std::vector<double>> tmp_vector_1 = cross(tmp_vector_b, tmp_vector_c);
	std::vector<std::vector<double>> tmp_vector_2 = cross(tmp_vector_c, tmp_vector_a);
	std::vector<std::vector<double>> tmp_vector_3 = cross(tmp_vector_a, tmp_vector_b);
	double lattice_spacing_a = cell_volume * 1.0 / fabs(norm(tmp_vector_1)); //@ remember to use fabs instead of abs for double values
	double lattice_spacing_b = cell_volume * 1.0 / fabs(norm(tmp_vector_2));
	double lattice_spacing_c = cell_volume * 1.0 / fabs(norm(tmp_vector_3));
	std::cout << "cell_volume and lattice spacings in a, b, c axes are: \n" << cell_volume << "  " << lattice_spacing_a 
		<< "  " << lattice_spacing_b << "  " << lattice_spacing_c << "\n";
	
	std::vector<std::string> filepath_vector(2615);
    std::string input_filepath = "E:\\Cal\\G16W\\MyProjects\\WORK\\fii\\test_output6.6.03\\Analysis\\POSCAR_filepath_all_GAPconv0.03.txt";
	std::ifstream myfileIN0;
	myfileIN0.open(input_filepath); //, ios_base::app); 这是缀加
    for (int loopi=0; loopi<2615; loopi++) {
        myfileIN0 >> filepath_vector[loopi];
    }
    myfileIN0.close();
    
    std::string output_file = "E:\\Cal\\G16W\\MyProjects\\WORK\\fii\\test_output6.6.03\\Analysis\\coordinate_stats_GAPconv0.03_preCellDelineation_v2.txt";
	std::ofstream myfileOUT0;
	myfileOUT0.open(output_file); //, std::ios_base::app); //这是缀加

    StructuralInformation structural_information;
    structural_information.lattice_parameter = lattice_parameters;
    structural_information.elements = {"P"};
    structural_information.n_elements = {256};
    structural_information.coordinates = std::vector<std::vector<double>>(256, std::vector<double>(3)); 
    structural_information.isDirect = false;
    structural_information.initialize_pair_graph();

    //@20240227 brought the preparation containers for Step 2 Specify which atom belong to which subcell and vice versa. out
    int subcell_num_a = lattice_spacing_a / (SHORT_VDW_UPPER_LIMIT + SHORT_VDW_SAFETY_CONST);
    int subcell_num_b = lattice_spacing_b / (SHORT_VDW_UPPER_LIMIT + SHORT_VDW_SAFETY_CONST);
    int subcell_num_c = lattice_spacing_c / (SHORT_VDW_UPPER_LIMIT + SHORT_VDW_SAFETY_CONST);
    if (subcell_num_a < 3 || subcell_num_b < 3 || subcell_num_c < 3 ) {
        std::cout << "Sorry, The Unit Cell Is Too Small To Get Divided Into Subcells. Current Algorithm Doesn't Support It.\n";
        return 1;
    }
    std::vector<std::vector<int>> atom2cell_list = std::vector<std::vector<int>>(256+1, std::vector<int>(3));
    std::vector<std::vector<std::vector<std::vector<int>>>> cell2atom_list = 
        std::vector<std::vector<std::vector<std::vector<int>>>>(subcell_num_a+2, 
            std::vector<std::vector<std::vector<int>>>(subcell_num_b+2, 
                std::vector<std::vector<int>>(subcell_num_c+2))); //, std::vector<double>(256)))); // 256 can be decreased in the future
    std::vector<std::vector<std::vector<int>>> cell2atom_list_each_size = 
        std::vector<std::vector<std::vector<int>>>(subcell_num_a+2, 
            std::vector<std::vector<int>>(subcell_num_b+2, std::vector<int>(subcell_num_c+2)));
    //reserve the vector seem no effects
    for (int ii=0; ii<=subcell_num_a+1; ii++) {
        for (int jj=0; jj<=subcell_num_b+1; jj++) {
            for (int kk=0; kk<=subcell_num_c+1; kk++) {
                //cell2atom_list[ii][jj][kk].reserve((totalAtomsNum/(subcell_num_a*subcell_num_b*subcell_num_c)+1) *3);
                cell2atom_list[ii][jj][kk].resize((structural_information.total_atoms_n/(subcell_num_a*subcell_num_b*subcell_num_c)+1) *3);
                cell2atom_list_each_size[ii][jj][kk] = 0;
            }
        }
    }



    std::vector<node> shortVdwList(300); //30 is enough according to stats; actually it should be below 22 for GAP optimized results.
    //Here we have to set it 300 for GMX-level.
	
	for (int loopi=0; loopi<2615; loopi++) {

        //auto start = std::chrono::high_resolution_clock::now();

		std::string input_filename = filepath_vector[loopi];
		//std::string input_filename = "D:\\Buzz\\FP_EM_fromScratch\\EM_round1\\Analysis\\test_ase\\2615_POSCAR.vasp";
		std::ifstream myfileIN;
		//myfileIN.open(input_filename); 

        try {
            myfileIN.open(input_filename);
            if (!myfileIN.is_open()) {  // Check if file opening failed
                throw std::runtime_error("File cannot be opened: " + input_filename);
            }
            //myfileIN.close();  // Close the file after processing
        } catch (const std::runtime_error& e) {
            std::cerr << "Exception caught: " << e.what() << std::endl;
            continue;
        }

		std::ofstream myfileOUT;
		myfileOUT.open(input_filename+"_GAPconv0.03_preCellDelineation_NoDynamicMemory_bondLength2.39_2.74_v2.txt"); 

        //auto stop = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by fileopen: " << duration.count() << " microseconds\n";
        // start = std::chrono::high_resolution_clock::now();
		
		// StructuralInformation structural_information;
	    // structural_information.lattice_parameter = lattice_parameters;
	    // structural_information.elements = {"P"};
	    // structural_information.n_elements = {256};
	    // structural_information.coordinates = std::vector<std::vector<double>>(256, std::vector<double>(3)); 
		// structural_information.isDirect = false;

        // stop = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        // std::cout << "Time taken by structural information initialization1: " << duration.count() << " microseconds\n";
        // start = std::chrono::high_resolution_clock::now();

        // structural_information.initialize_pair_graph();

        // stop = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        // std::cout << "Time taken by structural information initialization2: " << duration.count() << " microseconds\n";
        

        //start = std::chrono::high_resolution_clock::now();

        for (int i=0; i<=structural_information.total_atoms_n; i++) { //@20240227 3hours debug: this should be <=
            for (int j=0; j<=structural_information.total_atoms_n; j++) { //@20240227 3hours debug: this should be <=
                structural_information.pair_graph[i][j] = -1;
            }
        }
        structural_information.pair_size = 0; //It seems that there's no need to clear pair. 

        //stop = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by structural information ONLY pair graph initialization with -1: " << duration.count() << " microseconds\n";

        //myfileIN.close();// only for debug
        //start = std::chrono::high_resolution_clock::now();
	    
	    std::string tmp;
	    for (int i=1; i<=6; i++) {
	    	getline(myfileIN, tmp);
		}
		int atomsNum1; //, atomsNum2; 
		myfileIN >> atomsNum1; // >> atomsNum2;
		getline(myfileIN, tmp);
		getline(myfileIN, tmp);
		
		int totalAtomsNum = atomsNum1; //+atomsNum2;
	    //std::cout << totalAtomsNum << "\n" << tmp;
		
		for (int i=0; i<totalAtomsNum; i++) {
			for (int j=0; j<3; j++) {
				myfileIN >> structural_information.coordinates[i][j];
			}
		}
		myfileIN.close();

        //stop = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by file input: " << duration.count() << " microseconds\n";
        //start = std::chrono::high_resolution_clock::now();

		//New Features  
		// 1. Introduce Direct coordinates and update Cartesian coordinates
		std::vector<std::vector<double>> direct_coordinates = mmult(structural_information.coordinates, lattice_parameters_inverse);
		for (int i=0; i<totalAtomsNum; i++) {
			for (int j=0; j<3; j++) {
				while (direct_coordinates[i][j] > 1.0 + eps) {
					direct_coordinates[i][j] -= 1.0;
				}
				while (direct_coordinates[i][j] < 0.0 - eps) {
					direct_coordinates[i][j] += 1.0;
				}
			}
		}
		structural_information.coordinates = mmult(direct_coordinates, lattice_parameters);

        //stop = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by direct and cartesian coordinates conversion: " << duration.count() << " microseconds\n";
        //start = std::chrono::high_resolution_clock::now();
		
		// 2. Specify which atom belong to which subcell and vice versa.
		// int subcell_num_a = lattice_spacing_a / (SHORT_VDW_UPPER_LIMIT + SHORT_VDW_SAFETY_CONST);
		// int subcell_num_b = lattice_spacing_b / (SHORT_VDW_UPPER_LIMIT + SHORT_VDW_SAFETY_CONST);
		// int subcell_num_c = lattice_spacing_c / (SHORT_VDW_UPPER_LIMIT + SHORT_VDW_SAFETY_CONST);
		// if (subcell_num_a < 3 || subcell_num_b < 3 || subcell_num_c < 3 ) {
		// 	std::cout << "Sorry, The Unit Cell Is Too Small To Get Divided Into Subcells. Current Algorithm Doesn't Support It.\n";
		// 	return 1;
		// }
		// std::vector<std::vector<int>> atom2cell_list = std::vector<std::vector<int>>(256+1, std::vector<int>(3));
		// std::vector<std::vector<std::vector<std::vector<int>>>> cell2atom_list = 
		// 	std::vector<std::vector<std::vector<std::vector<int>>>>(subcell_num_a+2, 
		// 		std::vector<std::vector<std::vector<int>>>(subcell_num_b+2, 
		// 			std::vector<std::vector<int>>(subcell_num_c+2))); //, std::vector<double>(256)))); // 256 can be decreased in the future
		// std::vector<std::vector<std::vector<int>>> cell2atom_list_each_size = 
        //     std::vector<std::vector<std::vector<int>>>(subcell_num_a+2, 
		// 		std::vector<std::vector<int>>(subcell_num_b+2, std::vector<int>(subcell_num_c+2)));
        // //reserve the vector seem no effects
        for (int ii=0; ii<=subcell_num_a+1; ii++) {
            for (int jj=0; jj<=subcell_num_b+1; jj++) {
                for (int kk=0; kk<=subcell_num_c+1; kk++) {
                    //cell2atom_list[ii][jj][kk].reserve((totalAtomsNum/(subcell_num_a*subcell_num_b*subcell_num_c)+1) *3);
                    // cell2atom_list[ii][jj][kk].resize((totalAtomsNum/(subcell_num_a*subcell_num_b*subcell_num_c)+1) *3);
                    // 20240227 We don't resize every time. 
                    cell2atom_list_each_size[ii][jj][kk] = 0;
                }
            }
        }
        for (int i=1; i<=totalAtomsNum; i++) {
			atom2cell_list[i][0] = direct_coordinates[i-1][0] * subcell_num_a +1; 
			atom2cell_list[i][1] = direct_coordinates[i-1][1] * subcell_num_b +1; 
			atom2cell_list[i][2] = direct_coordinates[i-1][2] * subcell_num_c +1; 
			cell2atom_list[atom2cell_list[i][0]][atom2cell_list[i][1]][atom2cell_list[i][2]][cell2atom_list_each_size[atom2cell_list[i][0]][atom2cell_list[i][1]][atom2cell_list[i][2]]]=i;
            cell2atom_list_each_size[atom2cell_list[i][0]][atom2cell_list[i][1]][atom2cell_list[i][2]]++;
        }
		for (int ii=0; ii<=subcell_num_a+1; ii++) {
            //std::cout << "debug: ii = " << ii << "\n";
			int iii = ii;
			if (ii == 0) iii = subcell_num_a;
			if (ii == subcell_num_a+1) iii = 1;
			for (int jj=0; jj<=subcell_num_b+1; jj++) {
				int jjj = jj;
				if (jj == 0) jjj = subcell_num_b;
				if (jj == subcell_num_b+1) jjj = 1;
				for (int kk=0; kk<=subcell_num_c+1; kk++) {
					int kkk = kk;
					if (kk == 0) kkk = subcell_num_c;
					if (kk == subcell_num_c+1) kkk = 1;
					if (ii >=1 && ii <= subcell_num_a && jj >=1 && jj <= subcell_num_b && kk >=1 && kk <= subcell_num_c) continue;
					cell2atom_list[ii][jj][kk] = cell2atom_list[iii][jjj][kkk];
                    cell2atom_list_each_size[ii][jj][kk] = cell2atom_list_each_size[iii][jjj][kkk];
				}
			}
		}

        //stop = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by subcell and atom division: " << duration.count() << " microseconds\n";
        //start = std::chrono::high_resolution_clock::now();
	    
	    std::vector<int> coordinateStats(totalAtomsNum+1);
	    for (int i=1; i<=totalAtomsNum; i++) {
	    	coordinateStats[i] = 0;
		}
        //We moved this outside. This bug cost us an hour. Array boundary exceeded.
		//std::vector<node> shortVdwList(300); //30 is enough according to stats; actually it should be below 22 for GAP optimized results.
		
        int bondCnt = 0;
	    int shortVdwCnt = 0;
	    int coordinationNotThreeCnt = 0;
	    double tmpDist;
	    for (int i=1; i<=totalAtomsNum; i++) {
            //if (i == 56 || i==201 || i==255 || i==256) {
            //    std::cout << "debug : i = " << i << "\n";
            //}
	    	int thisSubcellI = atom2cell_list[i][0], thisSubcellJ = atom2cell_list[i][1], thisSubcellK = atom2cell_list[i][2]; 
	    	for (int ii=-1; ii<=1; ii++) {
	    		for (int jj=-1; jj<=1; jj++) {
	    			for (int kk=-1; kk<=1; kk++) {
	    				int thatSubcellI = thisSubcellI+ii;
	    				int thatSubcellJ = thisSubcellJ+jj;
	    				int thatSubcellK = thisSubcellK+kk;
	    				int iii=0, jjj=0, kkk=0;
	    				if (thatSubcellI == 0) iii = -1;
	    				if (thatSubcellJ == 0) jjj = -1;
	    				if (thatSubcellK == 0) kkk = -1;
	    				if (thatSubcellI == subcell_num_a+1) iii = +1;
	    				if (thatSubcellJ == subcell_num_b+1) jjj = +1;
	    				if (thatSubcellK == subcell_num_c+1) kkk = +1;
	    				int thatSubcellAtomNum = cell2atom_list_each_size[thatSubcellI][thatSubcellJ][thatSubcellK];
	    				
						for (int tmp_j=0; tmp_j<thatSubcellAtomNum; tmp_j++) {
	    					int j = cell2atom_list[thisSubcellI+ii][thisSubcellJ+jj][thisSubcellK+kk][tmp_j];
                            if (j<=i) continue;
                            //if (i == 56 && ii==-1 && jj==0 && kk==0) {
                            //    std::cout << "debug : i = " << i << "\n";
                            //}
	    					tmpDist = structural_information.dist(i, j, iii, jjj, kkk);
				    		if ( tmpDist <= SHORT_VDW_LOWER_LIMIT ) {
				    			coordinateStats[i]++;
				    			coordinateStats[j]++;
                                bondCnt++;
				    			myfileOUT << "  " << std::setw(3) << i << "  " << std::setw(3) << j << "  " << tmpDist << "\n";
							} else if (tmpDist < SHORT_VDW_UPPER_LIMIT) {
								shortVdwCnt++;
								shortVdwList[shortVdwCnt].atom_i = i;
								shortVdwList[shortVdwCnt].atom_j = j;
								shortVdwList[shortVdwCnt].atom_pair_dist = tmpDist;
							}
						}
	    			}
	    		}
	    	}
			if (coordinateStats[i] != 3) {
				coordinationNotThreeCnt++;
				myfileOUT0 << "   " << std::setw(3) << i << " Coordination Number: " << coordinateStats[i] << "\n";
			}
		}

        //stop = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by bond and shortVdW stats: " << duration.count() << " microseconds\n";
        //start = std::chrono::high_resolution_clock::now();
		
		myfileOUT << "\n---- shortVdwList ----\n";
		for (int i=1; i<=shortVdwCnt; i++) {
			myfileOUT << "  " << std::setw(3) << shortVdwList[i].atom_i << "  " << std::setw(3) << shortVdwList[i].atom_j << "  " << shortVdwList[i].atom_pair_dist << "\n";
		}

        myfileOUT << "\n";
        myfileOUT << "The total count of bond is: " << bondCnt << "\n";
        myfileOUT << "The total count of short vdW is: " << shortVdwCnt << "\n";
		myfileOUT0 << filepath_vector[loopi] << ": coordinationNotThreeCnt = " << coordinationNotThreeCnt << "\n" << "\n\n";
		myfileOUT.close();

        //stop = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        //std::cout << "Time taken by file output: " << duration.count() << " microseconds\n";
        if ((loopi+1) % 100 == 0)
            std::cout << "------Loop " << loopi+1 << " Finishes.------\n\n";
        
	}
	myfileOUT0.close();
	
	return 0;
} 
