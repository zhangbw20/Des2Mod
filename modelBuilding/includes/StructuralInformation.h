// StructuralInformation.h
#ifndef STRUCTURALINFORMATION_H
#define STRUCTURALINFORMATION_H

#include <vector>
#include <string>
#include <cmath>
#include ".\matrix_operations_using_vector.h"

#define pi 3.1415926535897932

extern const double INF; //declare but not define

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
	std::vector<int> atomCategory;
	std::vector<std::vector<int>> bondList;
	std::vector<double> bondLengthList;
	std::vector<std::vector<int>> angleList;
	std::vector<double> angleDegreeList;
	std::vector<std::vector<int>> vdwList;

	void initialize_coordinates_size(int numberOfAtoms);
    std::vector<int> determine_move_direction_for_second_atom(int atom1, int atom2);
	int find_pair_index(int atom1, int atom2);
    double dist(int atom1, int atom2);
    std::vector<double> nearest_average(int atom1, int atom2);
	double angle(int atom1, int atom2, int atom3); //只解决最近原子对间的夹角
	double dihedral_angle(int atom1, int atom2, int atom3, int atom4); //只解决最近原子对间的二面角
};

#endif // STRUCTURALINFORMATION_H
