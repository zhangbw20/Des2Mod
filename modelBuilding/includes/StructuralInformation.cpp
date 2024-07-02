// StructuralInformation.cpp

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include "StructuralInformation.h"

const double INF = 2147483647;

void StructuralInformation::initialize_coordinates_size(int numberOfAtoms) {
	coordinates = std::vector<std::vector<double>>(numberOfAtoms, std::vector<double>(3));
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
				cartesian_difference[cnt][0] = coordinates[atom1-1][0] - coordinates[atom2-1][0] + di*lattice_parameter[0][0] + dj*lattice_parameter[1][0] + dk*lattice_parameter[2][0];
				cartesian_difference[cnt][1] = coordinates[atom1-1][1] - coordinates[atom2-1][1] + di*lattice_parameter[0][1] + dj*lattice_parameter[1][1] + dk*lattice_parameter[2][1];
				cartesian_difference[cnt][2] = coordinates[atom1-1][2] - coordinates[atom2-1][2] + di*lattice_parameter[0][2] + dj*lattice_parameter[1][2] + dk*lattice_parameter[2][2];
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
	
	pair.push_back(tmp_pair);
	return tmp_move_direction_for_second_atom;
}

int StructuralInformation::find_pair_index(int atom1, int atom2) {
	bool pair_written_flag = false;
	int pair_index = -1;
	for (int i=0; i<pair.size(); i++) {
		if (pair[i].atom1 == atom1 && pair[i].atom2 == atom2) {
			pair_written_flag = true;
			pair_index = i;
		}
	}
	if (pair_written_flag == false) {
		determine_move_direction_for_second_atom(atom1, atom2);
		for (int i=0; i<pair.size(); i++) {
			if (pair[i].atom1 == atom1 && pair[i].atom2 == atom2) {
				pair_written_flag = true;
				pair_index = i;
			}
		}
	}
	return pair_index;
}


double StructuralInformation::dist(int atom1, int atom2) {// 
	int pair_index = find_pair_index(atom1, atom2);
	double cartesian_difference[3] = {0.0};
	cartesian_difference[0] = coordinates[atom1-1][0] - coordinates[atom2-1][0] + pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	cartesian_difference[1] = coordinates[atom1-1][1] - coordinates[atom2-1][1] + pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	cartesian_difference[2] = coordinates[atom1-1][2] - coordinates[atom2-1][2] + pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	double cartesian_distance = 0.0;
	for (int j=0; j<3; j++) {
		cartesian_distance += cartesian_difference[j] * cartesian_difference[j];
	}
	cartesian_distance = std::sqrt(cartesian_distance);
	pair[pair_index].nearest_distance = cartesian_distance;
	return cartesian_distance;
}

std::vector<double> StructuralInformation::nearest_average(int atom1, int atom2) { //≤ªÕÍ…∆£¨À≥–ÚΩªªª ±Œ¥≈–∂œ
	int pair_index = find_pair_index(atom1, atom2);
	double atom2_nearest_cartesian_coordinates[3] = {0.0};
	atom2_nearest_cartesian_coordinates[0] = coordinates[atom2-1][0] - pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][0] - pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][0] - pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	atom2_nearest_cartesian_coordinates[1] = coordinates[atom2-1][1] - pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][1] - pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][1] - pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	atom2_nearest_cartesian_coordinates[2] = coordinates[atom2-1][2] - pair[pair_index].move_direction_for_second_atom[0]*lattice_parameter[0][2] - pair[pair_index].move_direction_for_second_atom[1]*lattice_parameter[1][2] - pair[pair_index].move_direction_for_second_atom[2]*lattice_parameter[2][2];
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
	two_vectors[i][0] = coordinates[atom2-1][0] - coordinates[atom1-1][0] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	two_vectors[i][1] = coordinates[atom2-1][1] - coordinates[atom1-1][1] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	two_vectors[i][2] = coordinates[atom2-1][2] - coordinates[atom1-1][2] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	i = 1;
	two_vectors[i][0] = coordinates[atom2-1][0] - coordinates[atom3-1][0] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	two_vectors[i][1] = coordinates[atom2-1][1] - coordinates[atom3-1][1] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	two_vectors[i][2] = coordinates[atom2-1][2] - coordinates[atom3-1][2] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	
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
	three_vectors[i][0] = coordinates[atom1-1][0] - coordinates[atom2-1][0] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	three_vectors[i][1] = coordinates[atom1-1][1] - coordinates[atom2-1][1] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	three_vectors[i][2] = coordinates[atom1-1][2] - coordinates[atom2-1][2] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	i = 1;
	three_vectors[i][0] = coordinates[atom2-1][0] - coordinates[atom3-1][0] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	three_vectors[i][1] = coordinates[atom2-1][1] - coordinates[atom3-1][1] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	three_vectors[i][2] = coordinates[atom2-1][2] - coordinates[atom3-1][2] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	i = 2;
	three_vectors[i][0] = coordinates[atom3-1][0] - coordinates[atom4-1][0] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][0] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][0] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][0];
	three_vectors[i][1] = coordinates[atom3-1][1] - coordinates[atom4-1][1] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][1] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][1] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][1];
	three_vectors[i][2] = coordinates[atom3-1][2] - coordinates[atom4-1][2] + pair[pair_index[i]].move_direction_for_second_atom[0]*lattice_parameter[0][2] + pair[pair_index[i]].move_direction_for_second_atom[1]*lattice_parameter[1][2] + pair[pair_index[i]].move_direction_for_second_atom[2]*lattice_parameter[2][2];
	
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