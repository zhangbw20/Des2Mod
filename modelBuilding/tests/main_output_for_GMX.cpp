//when compiling use command: g++ -std=c++17 -o main_output_etot_poscar.exe main_output_etot_poscar.cpp ../includes/matrix_operations_using_vector.cpp ../includes/StructuralInformation.cpp
//formula_corr means correction of the "x### - x### + ###" part. Two Faults Corrected: Move direction when covering the boundary; vdW index mischosen when output.

// --------------------------------------------------------------------- //
//  N     N    OOO    TTTTTTT  EEEEEEE    Differences:
//  NN    N   O   O      T     E            LEFT, RIGHT was capitalized.
//  N N   N  O     O     T     E            Bond Length Unit: nm. 
//  N  N  N  O     O     T     EEEEEEE                  zbw 20240124
//  N   N N  O     O     T     E      
//  N    NN   O   O      T     E      
//  N     N    OOO       T     EEEEEEE
// --------------------------------------------------------------------- //



#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "./../includes/matrix_operations_using_vector.h"
#include "./../includes/StructuralInformation.h"
#include <fstream>
#include <iomanip>
//#include <filesystem>
#include <sstream>

using namespace std;

#define pi 3.1415926535897932
#define LEFT -1 //20240124 capitalized
#define RIGHT 1
#define DESCRIPTOR_CNT 2615
#define PCLU_GROUP 17 // 17 coordinates for consideration of a P_Cluster
//Noth that PCLU_GROUP == 17 cannot be casually changed here since vectors 
//later used may be initialized with RIGHT 17 elements.

string k_bond = "1.506240E+05";
string k_bond_1_2 = "7.531200E+04";
string k_angle = "5.020800E+02";
string k_angle_1_2 = "2.510400E+02";

struct PCLU { // abbreviation for P_Cluster
	int N_atoms = -1; // init as -1, then set to 8 or 9
	int twist_angle = -1; // init as int, for table lookup 
	double delicate_twist_angle = -1.0; // To locate P9 better, 
	// and solve overlapping problems so as to closer to FP
	double wavy_shift = 0.0; // c_direction
	double wavy_rotation = 0.0; // in_wOc_plane
	vector<vector<double>> coordinates = vector<vector<double>>(PCLU_GROUP, vector<double>(3, 0.0) ); //the coordinates relative to P0: 
	//0, 1~8, 9-10 reserve for P9, 11-14 in use for pre-determined P2, 15-16 for averaged P2
	// Detailed: 0 for P0, 10 for P_ortho in P9. For P8 these to overlap.
	// 11 -- 1, 12 -- 2, 13 -- 7, 14 -- 8;
    // Importantly, at molecular mechanics preparation, we overwrite 11 with linked atom 9 in the neighboring cluster; 13 and 14 for averaged P2 on the other side.
	double P2_deviations[3] = {-1.0, -1.0, -1.0};
    int overall_index[PCLU_GROUP] = {0};
};

//Parameters_Table
const double slice13_wavy_shift = -2.29; //c_direction, 10% of c axis
const double slice2_wavy_shift = 0.0; // @20231123 added slice2
const double slice1_wavy_rotation = -10.0; //in_wOc_plane set to 10 deg maybe? positive for anticlockwise 
const double slice2_wavy_rotation = 0.0; // @20231123 added slice2
const double slice3_wavy_rotation = +10.0; //Correct order: First twist, and then wavy rotation.

map<double, vector<vector<double>>> wavy_rotation_matrix_T;
map<double, vector<vector<double>>> wavy_translation_matrix_T;
//下两行可能会被抛弃
vector<vector<double>> wavy_rotation_matrix1_T(3, vector<double>(3) ) ;
vector<vector<double>> wavy_rotation_matrix3_T(3, vector<double>(3) ) ;

vector<vector<double>> center_position(25, vector<double>(3)); 
vector<vector<vector<double>>> center_position_3D_matrix_T(25, vector<vector<double>>(17, vector<double>(3))); 
vector<vector<double>> P8_initial_position(PCLU_GROUP, vector<double>(3)); //Paste to coordinates array after initialization.
vector<vector<double>> P9_initial_position(PCLU_GROUP, vector<double>(3)); //Paste to coordinates array after initialization.
vector<vector<double>> predict_P2_for_P8(3);
vector<vector<double>> predict_P2_for_P9(3);

vector<vector<vector<double>>> twist_rotation_matrix_T(361, vector<vector<double>>(3, vector<double>(3, 0.0) ) );
map<double, vector<vector<double>>> augmented_twist_rotation_matrix_T; //augmented twist matrix with angle being double

double cubic_difference_of_bond_length (double bond_length);
int int_abs(int n);
int neighbor_slice_cluster(int i0, int direction);
void initialize_center_position_matrix_T();
void initialize_P8_and_P9_position_matrix_T();
void initialize_twist_rotation_matrix_T();
void initialize_wavy_rotation_matrix_T(double wavy_rotation_theta);
void initialize_augmented_twist_rotation_matrix_T(double twist_rotation_theta);
void initialize_wavy_translation_matrix_T(double slice_wavy_shift);
int determine_P8_twist_angle();
int determine_P8_delicate_twist_angle();
void fill_bond_list(StructuralInformation& structural_information, vector<vector<bool>>& bondGraph, int& bondIndexCnt, int overall_index_k, int overall_index_l);

PCLU thisPCLU[25];
//w方向上三个磷团簇的相对位置关系可以用(x ± 8 + 24)%24处理 

//relatively less important information
vector<vector<double>> lattice_parameters = {
	{18.7999992371, 0.0, 0.0},
    {-4.2539077811, 12.4957698264, 0.0},
    {-3.0286706819, -0.6510739326, 22.6894962298}
};

int findLinked9Array[8+1][6+1] = {
	{0, 0, 0, 0, 0, 0, 0},
	{0, 6, 7, 1, 2, 3, 1},
	{0, 6, 5, -1, 2, 1, -1},
	{0, -2, -1, 1, 2, 3, 1},
	{0, -2, -3, -1, 2, 1, -1},
	{0, -2, -1, 1, 2, 3, 1},
	{0, -2, -3, -1, 2, 1, -1},
	{0, -2, -1, 1, -6, -5, 1},
	{0, -2, -3, -1, -6, -7, -1}
};
int findLinked9 (int thisClusterID, char thisDescriptor) { //find the neighboring cluster id of linked atom 9
	int thisClusterIDMod8 = (thisClusterID-1) % 8 +1 ;
	return thisClusterID + findLinked9Array[thisClusterIDMod8][thisDescriptor-'0']; 
}

std::vector<std::vector<int>> bondAndVdwSecondAtomMoveDirectionList;
std::vector<std::vector<double>> bondAndVdwSecondAtomMoveDirectionList_Double;
std::vector<std::vector<double>> bondAndVdwSecondAtomMoveCartesianDirectionList;
std::vector<std::vector<std::vector<double>>> moveDirectDirectionGraph;
std::vector<std::vector<std::vector<double>>> moveCartesianDirectionGraph;

int main () {

	initialize_center_position_matrix_T();

	initialize_P8_and_P9_position_matrix_T();

	initialize_twist_rotation_matrix_T();

	initialize_wavy_rotation_matrix_T(slice1_wavy_rotation);
	initialize_wavy_rotation_matrix_T(slice2_wavy_rotation);
	initialize_wavy_rotation_matrix_T(slice3_wavy_rotation);

	initialize_wavy_translation_matrix_T(slice13_wavy_shift);
	initialize_wavy_translation_matrix_T(slice2_wavy_shift);
	
	std::cout << "debug: what's happening?\n";

    vector<string> descriptor_vector(2615);
    string input_filename = "E:\\Cal\\G16W\\MyProjects\\WORK\\fii\\descriptor_vector_2615_input.txt";
	ifstream myfileIN;
	myfileIN.open(input_filename); //, ios_base::app); 这是缀加
    for (int loopi=0; loopi<2615; loopi++) {
        myfileIN >> descriptor_vector[loopi];
    }
    myfileIN.close();

    std::ifstream itp_header("E:\\Cal\\home\\yanqf\\zbw\\WORK\\test\\vscode\\gettingstarted\\descriptor2Model\\modelBuilding\\includes\\itp_header.txt");
    std::stringstream buffer;
    if (itp_header) {
        // Read the contents of template0.txt into a string buffer
        buffer << itp_header.rdbuf();
        itp_header.close();
    } else {
        std::cerr << "Failed to open E:\\Cal\\home\\yanqf\\zbw\\WORK\\test\\vscode\\gettingstarted\\descriptor2Model\\modelBuilding\\includes\\itp_header.txt\n";
        return 1;
    }

    int descriptor_cnt = 0;
	int model_cnt = 0;
	vector<int> model_number_of_each_descriptor;
	model_number_of_each_descriptor.resize(DESCRIPTOR_CNT);
	//尚未进循环
    //vector<double> GoF(2615);
    
    for (int loopi=3; loopi<2615; loopi++) { //debug
    //for (int loopi=0; loopi<2615; loopi++) {
        string descriptor0; // cin >> descriptor0; 
        descriptor0 = descriptor_vector[loopi];
        string descriptor = "_" + descriptor0.substr(0, 8) + descriptor0.substr(9, 8) + descriptor0.substr(18, 8);
        
        // Record the Twist Angle of P9
        for (int i=1; i<=8; i++) {
            if (descriptor[i] == '8') {
                thisPCLU[i].N_atoms = 8;
            } else { // descriptor[i] - '0' == 1~6
                thisPCLU[i].N_atoms = 9;
                thisPCLU[i].twist_angle = (descriptor[i] - '0') * 60;
                switch (descriptor[i] - '0') {
                    case 1: ;
                    case 4: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 18.40; //Close to FP P9 structure as much as possible
                        break;
                    case 2: ;
                    case 5: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 18.70;
                        break;
                    case 3: ;
                    case 6: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 19.25;
                        break;
                }
            }
        }

        for (int i=9; i<=16; i++) {
            if (descriptor[i] == '8') {
                thisPCLU[i].N_atoms = 8;
            } else { // descriptor[i] - '0' == 1~6
                thisPCLU[i].N_atoms = 9;
                thisPCLU[i].twist_angle = (descriptor[i] - '0') * 60;
                switch (descriptor[i] - '0') {
                    case 1: ;
                    case 4: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 17.65; //Close to FP P9 structure as much as possible
                        break;
                    case 2: ;
                    case 5: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 17.55;
                        break;
                    case 3: ;
                    case 6: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 19.55;
                        break;
                }
            }
        }

        for (int i=17; i<=24; i++) {
            if (descriptor[i] == '8') {
                thisPCLU[i].N_atoms = 8;
            } else { // descriptor[i] - '0' == 1~6
                thisPCLU[i].N_atoms = 9;
                thisPCLU[i].twist_angle = (descriptor[i] - '0') * 60;
                switch (descriptor[i] - '0') {
                    case 1: ;
                    case 4: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 11.45; //Close to FP P9 structure as much as possible
                        break;
                    case 2: ;
                    case 5: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 11.65;
                        break;
                    case 3: ;
                    case 6: ;
                        thisPCLU[i].delicate_twist_angle = thisPCLU[i].twist_angle + 19.85;
                        break;
                }
            }
        }

        // Only a quick definition of totalNumberOfAtoms
        int totalNumberOfAtoms = 48; //Initialize with number of P atoms on P2
        for (int i=1; i<=24; i++) {
            totalNumberOfAtoms += thisPCLU[i].N_atoms;
        }

        //Determine the Twist Angle of P8
        int this_model_number = determine_P8_delicate_twist_angle();
        model_number_of_each_descriptor[descriptor_cnt] = this_model_number; //From zero count to (DESCRIPTOR_CNT - 1)
        descriptor_cnt += 1;
        model_cnt += this_model_number;

        //Record the Wavy Shift and Rotation of Slice 1 and Slice 3
        for (int i=1; i<=8; i++) {
            thisPCLU[i].wavy_shift = slice13_wavy_shift;
            thisPCLU[i].wavy_rotation = slice1_wavy_rotation;
        }
        for (int i=17; i<=24; i++) {
            thisPCLU[i].wavy_shift = slice13_wavy_shift;
            thisPCLU[i].wavy_rotation = slice3_wavy_rotation;
        }
        
        //Rotation (and Translation) of All P8/P9 Units 
        for (int i=1; i<=24; i++) {
            vector<vector<double>> tmp_coordinates; //@ 20231123 P8 and P9 initial positions was not filled into thisPCLU
            if (!augmented_twist_rotation_matrix_T[thisPCLU[i].delicate_twist_angle].size()) {
                initialize_augmented_twist_rotation_matrix_T(thisPCLU[i].delicate_twist_angle);
            }
            if (thisPCLU[i].N_atoms == 8) {
                tmp_coordinates = mmult(P8_initial_position, augmented_twist_rotation_matrix_T[thisPCLU[i].delicate_twist_angle]);
            } else { // thisPCLU[i].N_atoms == 9
                tmp_coordinates = mmult(P9_initial_position, augmented_twist_rotation_matrix_T[thisPCLU[i].delicate_twist_angle]);
            }
            tmp_coordinates = mmult(tmp_coordinates, wavy_rotation_matrix_T[thisPCLU[i].wavy_rotation]);
            tmp_coordinates = madd(tmp_coordinates, wavy_translation_matrix_T[thisPCLU[i].wavy_shift]);
            thisPCLU[i].coordinates = madd(tmp_coordinates, center_position_3D_matrix_T[i]);
        }

        //Determine the P2 Positions according to its neighboring P8/P9 on both sides. (Not finished. 1/2)
        for (int i=1; i<=24; i++) {
            vector<vector<double>> tmp_P1 = mblock(thisPCLU[i].coordinates, 1+1, 0+1, 1, 3);
            vector<vector<double>> tmp_P2 = mblock(thisPCLU[i].coordinates, 2+1, 0+1, 1, 3);
            vector<vector<double>> tmp_P4 = mblock(thisPCLU[i].coordinates, 4+1, 0+1, 1, 3);
            vector<vector<double>> tmp_P5 = mblock(thisPCLU[i].coordinates, 5+1, 0+1, 1, 3);
            vector<vector<double>> tmp_P7 = mblock(thisPCLU[i].coordinates, 7+1, 0+1, 1, 3);
            vector<vector<double>> tmp_P8 = mblock(thisPCLU[i].coordinates, 8+1, 0+1, 1, 3);

            vector<vector<double>> tmp_P0 = mblock(thisPCLU[i].coordinates, 0+1, 0+1, 1, 3);//@ 20231123 忘记减P0了...
            tmp_P1 = msubtract(tmp_P1, tmp_P0);
            tmp_P2 = msubtract(tmp_P2, tmp_P0);
            tmp_P4 = msubtract(tmp_P4, tmp_P0);
            tmp_P5 = msubtract(tmp_P5, tmp_P0);
            tmp_P7 = msubtract(tmp_P7, tmp_P0);
            tmp_P8 = msubtract(tmp_P8, tmp_P0);
            
            vector<vector<double>> tmp_add_P4_P5 = madd(tmp_P4, tmp_P5);
            vector<vector<double>> tmp_pre_P2_11_coordinates = {tmp_P1[0], tmp_P2[0], tmp_add_P4_P5[0]};
            vector<vector<double>> tmp_pre_P2_12_coordinates = {tmp_P2[0], tmp_P1[0], tmp_add_P4_P5[0]};
            vector<vector<double>> tmp_pre_P2_13_coordinates = {tmp_P7[0], tmp_P8[0], tmp_add_P4_P5[0]};
            vector<vector<double>> tmp_pre_P2_14_coordinates = {tmp_P8[0], tmp_P7[0], tmp_add_P4_P5[0]};
            if (thisPCLU[i].N_atoms == 8) {
                tmp_pre_P2_11_coordinates = mmult(predict_P2_for_P8, tmp_pre_P2_11_coordinates);
                tmp_pre_P2_12_coordinates = mmult(predict_P2_for_P8, tmp_pre_P2_12_coordinates);
                tmp_pre_P2_13_coordinates = mmult(predict_P2_for_P8, tmp_pre_P2_13_coordinates);
                tmp_pre_P2_14_coordinates = mmult(predict_P2_for_P8, tmp_pre_P2_14_coordinates);
            } else { //thisPCLU[i].N_atoms == 9
                tmp_pre_P2_11_coordinates = mmult(predict_P2_for_P9, tmp_pre_P2_11_coordinates);
                tmp_pre_P2_12_coordinates = mmult(predict_P2_for_P9, tmp_pre_P2_12_coordinates);
                tmp_pre_P2_13_coordinates = mmult(predict_P2_for_P9, tmp_pre_P2_13_coordinates);
                tmp_pre_P2_14_coordinates = mmult(predict_P2_for_P9, tmp_pre_P2_14_coordinates);
            }
            tmp_pre_P2_11_coordinates = madd(tmp_pre_P2_11_coordinates, tmp_P0); //@20231123还得变换回去
            tmp_pre_P2_12_coordinates = madd(tmp_pre_P2_12_coordinates, tmp_P0);
            tmp_pre_P2_13_coordinates = madd(tmp_pre_P2_13_coordinates, tmp_P0);
            tmp_pre_P2_14_coordinates = madd(tmp_pre_P2_14_coordinates, tmp_P0);

            thisPCLU[i].coordinates[11] = tmp_pre_P2_11_coordinates[0]; 
            thisPCLU[i].coordinates[12] = tmp_pre_P2_12_coordinates[0];
            thisPCLU[i].coordinates[13] = tmp_pre_P2_13_coordinates[0];
            thisPCLU[i].coordinates[14] = tmp_pre_P2_14_coordinates[0];
        }
        
        
        //Check the distance of the two P2s. Average them if in rational range. (Done. 2/2)
            //Here we may need to adjust the wavy_shift and wavy_rotation values. Let P atoms being not too close nor too far away. 
        /*double SumOfSuqaredDeviations = 0.0;
        double GoF = 0.0;
        vector<double> bond_length_distribution(13);
        for (int bld = 0; bld < 13; bld++) {
            bond_length_distribution[bld] = 0;
        }*/
        for (int i=1; i<=24; i++) {
            // Note that P1 and P2 are at the RIGHT side of this PCLU, so we put their neighboring 
            // two atoms into this struct.
            vector<vector<double>> P1_neighbor = mblock(thisPCLU[i].coordinates, 11+1, 0+1, 1, 3); 
            vector<vector<double>> P2_neighbor = mblock(thisPCLU[i].coordinates, 12+1, 0+1, 1, 3); 
            vector<vector<double>> RIGHT_P7_neighbor = mblock(thisPCLU[neighbor_slice_cluster(i, RIGHT)].coordinates, 13+1, 0+1, 1, 3); 
            vector<vector<double>> RIGHT_P8_neighbor = mblock(thisPCLU[neighbor_slice_cluster(i, RIGHT)].coordinates, 14+1, 0+1, 1, 3); 
            
            StructuralInformation sub_structural_information;
            sub_structural_information.lattice_parameter = lattice_parameters;
            sub_structural_information.elements = {"P"};
            sub_structural_information.n_elements = {4};
            sub_structural_information.coordinates = {
                P1_neighbor[0],
                P2_neighbor[0],
                RIGHT_P7_neighbor[0],
                RIGHT_P8_neighbor[0]
            };
            sub_structural_information.isDirect = false;
            
            
            thisPCLU[i].P2_deviations[1] = sub_structural_information.dist(1, 3);
            thisPCLU[i].P2_deviations[2] = sub_structural_information.dist(2, 4);
            /*SumOfSuqaredDeviations += thisPCLU[i].P2_deviations[1] * thisPCLU[i].P2_deviations[1];
            SumOfSuqaredDeviations += thisPCLU[i].P2_deviations[2] * thisPCLU[i].P2_deviations[2];*/
            if (1) { //记得写偏差判断语句和调整偏差语句。
                vector<double> tmpP15 = sub_structural_information.nearest_average(1, 3);
                vector<double> tmpP16 = sub_structural_information.nearest_average(2, 4);
                thisPCLU[i].coordinates[15] = tmpP15;
                thisPCLU[i].coordinates[16] = tmpP16;
            }
        }

        

        /*
        for (int i=1; i<=24; i++) {
            // Evaluate GoF
            vector<vector<double>> thisP1 = mblock(thisPCLU[i].coordinates, 1+1, 0+1, 1, 3); 
            vector<vector<double>> thisP2 = mblock(thisPCLU[i].coordinates, 2+1, 0+1, 1, 3); 
            vector<vector<double>> RIGHTP7 = mblock(thisPCLU[neighbor_slice_cluster(i, RIGHT)].coordinates, 7+1, 0+1, 1, 3); 
            vector<vector<double>> RIGHTP8 = mblock(thisPCLU[neighbor_slice_cluster(i, RIGHT)].coordinates, 8+1, 0+1, 1, 3);
            vector<vector<double>> thisP15 = mblock(thisPCLU[i].coordinates, 15+1, 0+1, 1, 3); 
            vector<vector<double>> thisP16 = mblock(thisPCLU[i].coordinates, 16+1, 0+1, 1, 3); 
            
            StructuralInformation sub_structural_information;
            sub_structural_information.lattice_parameter = lattice_parameters;
            sub_structural_information.elements = {"P"};
            sub_structural_information.n_elements = {6};
            sub_structural_information.coordinates = {
                thisP1[0],
                thisP2[0],
                RIGHTP7[0],
                RIGHTP8[0],
                thisP15[0],
                thisP16[0]
            };
            sub_structural_information.isDirect = false;   

            
            pair<int, int> couples[] = {{1,5}, {2,6}, {3,5}, {4,6}, {5,6}};
            for (const auto& couple : couples) {
                double tmp_bld = sub_structural_information.dist(couple.first, couple.second);
                GoF += cubic_difference_of_bond_length(tmp_bld);
                if (tmp_bld > 3.1) {
                    bond_length_distribution[12]++;
                } else if (tmp_bld > 2.9) {
                    bond_length_distribution[11]++;
                } else if (tmp_bld > 2.7) {
                    bond_length_distribution[10]++;
                } else if (tmp_bld > 2.5) {
                    bond_length_distribution[9]++;
                } else if (tmp_bld > 2.3) {
                    bond_length_distribution[8]++;
                } else if (tmp_bld >= 2.15) {
                    bond_length_distribution[7]++; //Perfect Spot
                } else if (tmp_bld >= 1.95) {
                    bond_length_distribution[6]++;
                } else if (tmp_bld >= 1.75) {
                    bond_length_distribution[5]++;
                } else if (tmp_bld >= 1.55) {
                    bond_length_distribution[4]++;
                } else if (tmp_bld >= 1.35) {
                    bond_length_distribution[3]++;
                } else if (tmp_bld >= 1.15) {
                    bond_length_distribution[2]++;
                } else if (tmp_bld >= 0.95) {
                    bond_length_distribution[1]++;
                } else { // tmp_bld < 0.95
                    bond_length_distribution[0]++;
                }
            }
        }


        double rootSSD = sqrt(SumOfSuqaredDeviations); //There's no need to get divided by Total number of atoms, since no matter the ratio of P8 to P9, the pairs of P2 is the same. 
        GoF = cbrt(GoF); //cubic root
        std::cout << descriptor0 << ":   " << rootSSD << "   " << GoF;
        for (int bld = 0; bld < 13; bld++) {
            std::cout << "  " << bond_length_distribution[bld];
        } 
        std::cout << "\n";
        */


        // ## Prepare general structural information, including bond, angle, vdw

        StructuralInformation structural_information;
        structural_information.lattice_parameter = lattice_parameters;
        structural_information.elements = {"P"};
        structural_information.n_elements = {totalNumberOfAtoms};
        structural_information.initialize_coordinates_size(totalNumberOfAtoms); //Here since only P exists, then we directly use totalNumberOfAtoms.
        structural_information.isDirect = false;

        structural_information.atomCategory.resize(totalNumberOfAtoms);
        structural_information.bondList = vector<vector<int>>(totalNumberOfAtoms*3/2, vector<int>(2)); //for phosphorus; the same below.
        structural_information.bondLengthList.resize(totalNumberOfAtoms*3/2);
        structural_information.angleList = vector<vector<int>>(totalNumberOfAtoms*3, vector<int>(3)); //(Format: center atom at middle)
        structural_information.angleDegreeList.resize(totalNumberOfAtoms*3);
        //structural_information.vdwList's size cannot be determined now. 
        structural_information.vdwList = vector<vector<int>> (totalNumberOfAtoms*(totalNumberOfAtoms-1)/2, vector<int>(2));
        
        
        //Fill coordinates, atomCategory, with detemining the overall index
        int overallIndexCnt = 0;
        for (int i=1; i<=24; i++) {
            vector<int> index_list;
            int index_list_length = -1;
            vector<int> atom_category_list;
            if (descriptor[i] == '8') {
                index_list = {1, 2, 3, 4, 5, 6, 7, 8, 15, 16};
                index_list_length = 8;
                atom_category_list = {8, 8, 6, 7, 7, 6, 8, 8, 2, 2};
            } else {
                index_list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16};
                index_list_length = 9;
                atom_category_list = {9, 9, 0, 1, 1, 0, 9, 9, 4, 2, 2};
            }
            for (int j=0; j<index_list_length; j++) {
                overallIndexCnt++ ;
                thisPCLU[i].overall_index[index_list[j]] = overallIndexCnt;
                structural_information.coordinates[overallIndexCnt-1] = thisPCLU[i].coordinates[index_list[j]];
                structural_information.atomCategory[overallIndexCnt-1] = atom_category_list[j];
            }
        }
        //Don't forget P2
        for (int i=1; i<=24; i++) {
            for (int k=15; k<=16; k++) {
                overallIndexCnt++ ;
                thisPCLU[i].overall_index[k] = overallIndexCnt;
                structural_information.coordinates[overallIndexCnt-1] = thisPCLU[i].coordinates[k];
                structural_information.atomCategory[overallIndexCnt-1] = 2;
            }
        }

        //Checkpoint
        if (overallIndexCnt == totalNumberOfAtoms) {
            std::cout << "Checked that overallIndexCnt is the same as totalNumberOfAtoms: " << loopi+1 << "/" << 2615 << "\n";
        } else {
            std::cout << "WARNING: Checked that overallIndexCnt is different from totalNumberOfAtoms: " << loopi+1 << "/" << 2615 << "\n";
            std::cout << "overallIndexCnt = " << overallIndexCnt << ", totalNumberOfAtoms = " << totalNumberOfAtoms << "\n";
        }
        

        //Fill bond, angle and vdw parameters
            // Fill bond pair 
        vector<vector<bool>> bondGraph = vector<vector<bool>>(totalNumberOfAtoms+1, vector<bool>(totalNumberOfAtoms+1));
        int bondIndexCnt = 0;
        for (int i=1; i<=24; i++) {
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[3]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[4]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[15]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[3]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[5]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[16]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[5]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[7]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[8]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[7]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[8]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[15], thisPCLU[i].overall_index[16]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[7], thisPCLU[neighbor_slice_cluster(i, LEFT)].overall_index[15]);
            fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[8], thisPCLU[neighbor_slice_cluster(i, LEFT)].overall_index[16]);
            if (descriptor[i] == '8') {
                fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[6]);
            } else {
                fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[9]);
                fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[9]);
                if ( i < findLinked9(i, descriptor[i]) ) { //neigboring cluster id linked by Atom 9 larger
                    fill_bond_list(structural_information, bondGraph, bondIndexCnt, thisPCLU[i].overall_index[9], thisPCLU[findLinked9(i, descriptor[i])].overall_index[9]);
                }
            }
        }

        //Fill angle pair
        int angleIndexCnt = 0;
        for (int i=1; i<=24; i++) {
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[4]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[15]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[15]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[5]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[16]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[16]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[2]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[5]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[7]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[7]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[4]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[8]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[8]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[7], thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[8]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[7], thisPCLU[i].overall_index[6]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[8], thisPCLU[i].overall_index[6]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[15], thisPCLU[i].overall_index[16]}; //Avoid Copying Fault: Not 2, should be 16
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[16], thisPCLU[i].overall_index[15]}; //Avoid Copying Fault: Not 1, should be 15
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[4], thisPCLU[i].overall_index[7], thisPCLU[neighbor_slice_cluster(i, LEFT)].overall_index[15]}; 
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[7], thisPCLU[neighbor_slice_cluster(i, LEFT)].overall_index[15]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[5], thisPCLU[i].overall_index[8], thisPCLU[neighbor_slice_cluster(i, LEFT)].overall_index[16]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[8], thisPCLU[neighbor_slice_cluster(i, LEFT)].overall_index[16]};
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[15], thisPCLU[neighbor_slice_cluster(i, RIGHT)].overall_index[7]}; //R7 should be RIGHT 7 instead of RIGHT 15
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[16], thisPCLU[i].overall_index[15], thisPCLU[neighbor_slice_cluster(i, RIGHT)].overall_index[7]}; //R7 should be RIGHT 7 instead of RIGHT 15
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[16], thisPCLU[neighbor_slice_cluster(i, RIGHT)].overall_index[8]}; //R8 should be RIGHT 8 instead of RIGHT 16
            structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[15], thisPCLU[i].overall_index[16], thisPCLU[neighbor_slice_cluster(i, RIGHT)].overall_index[8]}; //R8 should be RIGHT 8 instead of RIGHT 16

            if (descriptor[i] == '8') {
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[6]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[6]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[7]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[8]};    
            } else {
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[1], thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[9]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[2], thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[9]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[7], thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[9]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[8], thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[9]};  
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[9], thisPCLU[i].overall_index[6]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[3], thisPCLU[i].overall_index[9], thisPCLU[findLinked9(i, descriptor[i])].overall_index[9]};
                structural_information.angleList[angleIndexCnt++] = {thisPCLU[i].overall_index[6], thisPCLU[i].overall_index[9], thisPCLU[findLinked9(i, descriptor[i])].overall_index[9]};
            }
        }

        //Fill vdw pair: We don't fill vdW here
        // int vdwIndexCnt = 0;
        // for (int i=1; i<=totalNumberOfAtoms; i++) {
        //     for (int j=i+1; j<=totalNumberOfAtoms; j++) {
        //         if (bondGraph[i][j] == true) continue;
        //         if (structural_information.dist(i, j) < 5.0) {
        //             structural_information.vdwList[vdwIndexCnt++] = {i, j};
        //         }
        //     }
        // }


        //Fill bond length
        int k=-1, l=-1, m=-1;
        int ck=-1, cl=-1, cm=-1; //category of k and l
        int cmin=-1, cmax=-1; 
        for (int i=0; i<bondIndexCnt; i++) {
            k = structural_information.bondList[i][0];
            l = structural_information.bondList[i][1];       
            ck = structural_information.atomCategory[k-1];
            cl = structural_information.atomCategory[l-1];
            if (ck <= cl) {
                cmin = ck;
                cmax = cl;
            } else {
                cmin = cl;
                cmax = ck;
            }
            if (cmin == 1 && cmax == 1) {
                structural_information.bondLengthList[i] = 0.2269; // Unit was turned into nm 20240124
            } else if (cmin == 2 && cmax == 2) {
                structural_information.bondLengthList[i] = 0.2267;
            } else if (cmin == 2 && (cmax == 8 || cmax == 9)) {
                structural_information.bondLengthList[i] = 0.2198;
            } else if ( (cmin == 6 && cmax == 6) || (cmin == 7 && cmax == 7) ) {
                structural_information.bondLengthList[i] = 0.2296;
            } else { // all the rest 8(67)s && 9(01)s
                structural_information.bondLengthList[i] = 0.2211;
            }
        }


        //Fill angle degree
        for (int i=0; i<angleIndexCnt; i++) {
            k = structural_information.angleList[i][0];
            l = structural_information.angleList[i][1];
            m = structural_information.angleList[i][2];
            ck = structural_information.atomCategory[k-1];
            cl = structural_information.atomCategory[l-1];
            cm = structural_information.atomCategory[m-1];
            if (ck <= cm) {
                cmin = ck;
                cmax = cm;
            } else {
                cmin = cm;
                cmax = ck;
            }
            if (cl == 2) {
                if (cmin == 2 && (cmax==8||cmax==9) ) { //228 or 229 CenterWasPutAtFirstInComments
                    structural_information.angleDegreeList[i] = 104.55;
                } else if ( (cmin == 8 && cmax == 8) || (cmin == 8 && cmax == 9) || (cmin == 9 && cmax == 9) ) { //289, [288 and 299] 20240129 added
                    structural_information.angleDegreeList[i] = 98.92;
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 8) {
                if (cmin == 2) {
                    if (cmax == 7) { //827
                        structural_information.angleDegreeList[i] = 95.92;
                    } else if (cmax == 6) { //826
                        structural_information.angleDegreeList[i] = 107.47;
                    } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
                } else if (cmin == 6 && cmax==7) { //867
                    structural_information.angleDegreeList[i] = 97.50;
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 7) {
                if (cmin == 7 && cmax==8) { //778 // (cmax == 7) { //787 <- @20240124 Bug Detected
                    structural_information.angleDegreeList[i] = 104.17;
                } else if (cmin == 8 && cmax==8) { //788
                    structural_information.angleDegreeList[i] = 87.63;
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 6) {
                if (cmin == 6 && cmax==8) { //668
                    structural_information.angleDegreeList[i] = 99.62;
                } else if (cmin == 8 && cmax==8) { //688
                    structural_information.angleDegreeList[i] = 97.89;
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 9) {
                if (cmin == 1 && cmax==2) { //912
                    structural_information.angleDegreeList[i] = 94.34;
                } else if (cmin == 0) { 
                    if (cmax == 2) { //902
                        structural_information.angleDegreeList[i] = 98.07;
                    } else if (cmax == 1) { //901
                        structural_information.angleDegreeList[i] = 106.19;
                    } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 1) {
                if (cmin == 1 && cmax==9) { //119
                    structural_information.angleDegreeList[i] = 104.53;
                } else if (cmin == 9 && cmax==9) { //199
                    structural_information.angleDegreeList[i] = 94.93;
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 0) {
                if (cmin == 4 && cmax==9) { //049  20240124 error detected -> (cmin == 0) { //009
                    structural_information.angleDegreeList[i] = 106.68;
                } else if (cmin == 9 && cmax==9) { //099
                    structural_information.angleDegreeList[i] = 98.98;
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else if (cl == 4) { // cl == 4
                if (cmax == 0 && cmin==0) { //400  20240124 error detected -> (cmin == 9) { //499
                    structural_information.angleDegreeList[i] = 111.83;
                } else if (cmax == 4 && cmin==0) { //404  20240124 error detected -> //449
                    structural_information.angleDegreeList[i] = 94.15;
                    //cout << "This branch was entered.\n";
                } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
            } else cout << "Error: Wrong angle type detected. " << "center-atom1-atom2: " << cl << cmin << cmax << "\n";
        }
        
        vector<vector<int>> bondAndVdwList = structural_information.bondList;
        //bondAndVdwList.insert(bondAndVdwList.end(), structural_information.vdwList.begin(), structural_information.vdwList.begin() + vdwIndexCnt);
        bondAndVdwSecondAtomMoveDirectionList.resize(bondAndVdwList.size());
        bondAndVdwSecondAtomMoveDirectionList_Double.resize(bondAndVdwList.size());
        bondAndVdwSecondAtomMoveCartesianDirectionList.resize(bondAndVdwList.size());
        moveDirectDirectionGraph.resize(totalNumberOfAtoms);
        moveCartesianDirectionGraph.resize(totalNumberOfAtoms);
        for (int i=0; i<totalNumberOfAtoms; i++) {
            moveDirectDirectionGraph[i].resize(totalNumberOfAtoms);
            moveCartesianDirectionGraph[i].resize(totalNumberOfAtoms);
            for (int j=0; j<totalNumberOfAtoms; j++) {
            	moveDirectDirectionGraph[i][j].resize(3);
                moveCartesianDirectionGraph[i][j].resize(3);
            }
        }
        
        for (int i=0; i<bondAndVdwList.size(); i++) {
            // if (i==477) {
            //     bool thisBreakPoint = true;
            // }
            // if (i==478) {
            //     bool thisBreakPoint = true;
            // }
            bondAndVdwSecondAtomMoveDirectionList[i] = structural_information.determine_move_direction_for_second_atom(
                bondAndVdwList[i][0], bondAndVdwList[i][1]);
            bondAndVdwSecondAtomMoveDirectionList_Double[i].resize(3);
            for (int j=0; j<3; j++) {
                bondAndVdwSecondAtomMoveDirectionList_Double[i][j] = bondAndVdwSecondAtomMoveDirectionList[i][j];
            }
            
            //deal with graph 
            moveDirectDirectionGraph[bondAndVdwList[i][0]-1][bondAndVdwList[i][1]-1] = bondAndVdwSecondAtomMoveDirectionList_Double[i];
            moveDirectDirectionGraph[bondAndVdwList[i][1]-1][bondAndVdwList[i][0]-1] = bondAndVdwSecondAtomMoveDirectionList_Double[i];
            for (int j=0; j<3; j++) {
                moveDirectDirectionGraph[bondAndVdwList[i][1]-1][bondAndVdwList[i][0]-1][j] = - moveDirectDirectionGraph[bondAndVdwList[i][1]-1][bondAndVdwList[i][0]-1][j];
            }
            
            std::vector<std::vector<double>> tmpMatrix1 = mblock(bondAndVdwSecondAtomMoveDirectionList_Double, i+1, 1, 1, 3);
            std::vector<std::vector<double>> tmpMatrix3 = mmult(tmpMatrix1, structural_information.lattice_parameter);
            
            //std::cout << bondAndVdwList[i][0] << " " << bondAndVdwList[i][1];
            bondAndVdwSecondAtomMoveCartesianDirectionList[i].resize(3); //@20240116 This sentence should certainly be outside the loop of j!
            for (int j=0; j<3; j++) {
                
                //std::cout << "debug: " << tmpMatrix3.size() << " " << tmpMatrix3[0].size();
                bondAndVdwSecondAtomMoveCartesianDirectionList[i][j] = tmpMatrix3[0][j];
                //std::cout << " " << bondAndVdwSecondAtomMoveCartesianDirectionList[i][j]; 
            }
            
            //moveCartesianDirectionGraph[bondAndVdwList[i][0]-1][bondAndVdwList[i][1]-1].resize(3);
            //moveCartesianDirectionGraph[bondAndVdwList[i][1]-1][bondAndVdwList[i][0]-1].resize(3);
            for (int j=0; j<3; j++) {
                moveCartesianDirectionGraph[bondAndVdwList[i][0]-1][bondAndVdwList[i][1]-1][j] = tmpMatrix3[0][j];
                moveCartesianDirectionGraph[bondAndVdwList[i][1]-1][bondAndVdwList[i][0]-1][j] = - tmpMatrix3[0][j];
            }
            
            
            //std::cout << "\n";
        }


        // Write to ase file
            // One file test and then apply to all.
                //First let's try to write to vasp.
        
        
        //Quick Output 
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << (loopi+1);
        string dirpath = "C:\\Users\\zhbw\\Downloads\\test_output6\\" + oss.str() + "_" + descriptor0; //Don't forget the "_"

	//Output Etot Expression
		//3k-2, 3k-1, 3k
	//int k=-1, l=-1, m=-1;


    // Write Characteristic Parts


    // Output pdb
    ofstream OUTpdb;
    OUTpdb.open(dirpath + "\\" + "FormII_GMX.pdb");
    string pdb_line1 = "CRYST1   18.800   13.200   22.900  89.10  97.60 108.80 P 1           1";
    OUTpdb << pdb_line1 << "\n";
    for (int i=0; i<totalNumberOfAtoms; i++) {
		//for (int j=0; j<3; j++) {
            OUTpdb << "ATOM  " << std::setw(5) << std::right << (i+1) << " " << "PK" << structural_information.atomCategory[i] << "               " << std::fixed << std::setprecision(3) << std::setw(8) << structural_information.coordinates[i][0] << std::setw(8) << structural_information.coordinates[i][1] << std::setw(8) << structural_information.coordinates[i][2] << "  1.00  0.00           P" << "\n";
		//}
	}
    //HETATM    1 P1                   6.638   3.284   1.230  1.00  0.00           P
    OUTpdb << "END\n";
    OUTpdb.close();


    //Output itp
    ofstream OUTitp;
    string filepath_and_itpname = dirpath + "\\" + "FormII_GMX.itp"; 
    OUTitp.open(filepath_and_itpname);

    // Write header
    // OUTitp << itp_header.rdbuf();
    OUTitp << buffer.str(); //rdbuf can't be used multiple times. 

    // Write atoms
    OUTitp << "[ atoms ]" << "\n";
    OUTitp << ";  Index   type   residue  resname   atom        cgnr     charge       mass" << "\n";
    //     1     PK9        1      MOL     P1            1    0.00000000   30.973762
    for (int i=0; i<totalNumberOfAtoms; i++) {
		//for (int j=0; j<3; j++) {
            OUTitp << " " << std::setw(5) << std::right << (i+1) << "     " << "PK" << structural_information.atomCategory[i] << "        1      MOL     " << "PK" << structural_information.atomCategory[i] << "       " << std::setw(5) << std::right << (i+1) << "    0.00000000   30.973762" << "\n";
		//}
	}
    OUTitp << " " << "\n";

    // Write bonds
    OUTitp << "[ bonds ]" << "\n";
    OUTitp << "; atom_i  atom_j  functype      r0 (nm)    k (kJ/mol/nm^2)" << "\n";
    //    1       5         1        0.221100     1.506240E+05     ; P1-P3, prebuilt PK9-PK1
    for (int i=0; i<bondIndexCnt; i++) {
		k = structural_information.bondList[i][0];
		l = structural_information.bondList[i][1];
        OUTitp << std::setw(5) << std::right << k << "   " << std::setw(5) << std::right << l << "         1        " << std::fixed << std::setprecision(6) << std::setw(8) << structural_information.bondLengthList[i] << "     ";
        if (structural_information.atomCategory[k-1] == 2 || structural_information.atomCategory[l-1] == 2 ) {
            OUTitp << k_bond_1_2;
        } else {
            OUTitp << k_bond;
        }
        OUTitp << "     ";
        OUTitp << "; " << "PK" << structural_information.atomCategory[k-1] << "-" << "PK" << structural_information.atomCategory[l-1] << "\n"; 
	}
    OUTitp << " " << "\n";

    // Write angles
    OUTitp << "[ angles ]" << "\n";
    OUTitp << "; atom_i  atom_j  atom_k  functype    a0 (Deg.)  k (kJ/mol/rad^2)" << "\n";
    //    5       1      29         1        94.340      5.020800E+02     ; P3-P1-P15, prebuilt PK1-PK9-PK2
	for (int i=0; i<angleIndexCnt; i++) {
		k = structural_information.angleList[i][0];
		l = structural_information.angleList[i][1];
		m = structural_information.angleList[i][2];
		OUTitp << std::setw(5) << std::right << k << "   " << std::setw(5) << std::right << l << "   " << std::setw(5) << std::right << m << "         1      " << std::fixed << std::setprecision(3) << std::setw(8) << structural_information.angleDegreeList[i] << "      ";;
        if (structural_information.atomCategory[k-1] == 2 || structural_information.atomCategory[l-1] == 2 || structural_information.atomCategory[m-1] == 2) {
            OUTitp << k_angle_1_2;
        } else {
            OUTitp << k_angle;
        }
        OUTitp << "     ";
        OUTitp << "; " << "PK" << structural_information.atomCategory[k-1] << "-" << "PK" << structural_information.atomCategory[l-1] << "-" << "PK" << structural_information.atomCategory[m-1] << "\n"; 
	}
    OUTitp << " " << "\n";

    OUTitp << "; No pair needs to generate" << "\n\n";

    OUTitp << "; No improper needs to generate" << "\n";

    OUTitp.close();

        
    } // end of the loopi LOOP
	



	string wait_any_input_and_enter_string;
	std::cout << "Please input anything and then enter to end this program.\n";
	std::cin >> wait_any_input_and_enter_string;


	return 0;
} 

//--------------------------//

double cubic_difference_of_bond_length (double bond_length) {
    double cubic = 0.0;
    double difference;
    if ( bond_length > 2.3 ) {
        difference = bond_length - 2.3;
    } else if ( bond_length < 2.15 ) {
        difference = 2.15 - bond_length;
    }
    cubic = difference * difference * difference;
    return cubic;
}

int int_abs(int n) {
	if ( n >= 0 ) 
		return n;
	else // n < 0
		return -n;
}

void initialize_center_position_matrix_T() {
	center_position = {
	{0.0, 0.0, 0.0}, //this line is only for uccupying the position
	{0.483246453, -0.569689691, 19.853309201},
	{-1.643707438, 5.678195222, 19.853309201},
	{-1.950016713, 8.964906162, 14.180935144},
	{0.176937178, 2.717021249, 14.180935144},
	{-0.129372097, 6.003732188, 8.508561086},
	{1.997581794, -0.244152725, 8.508561086},
	{1.691272519, 3.042558215, 2.836187029},
	{-0.435681372, 9.290443128, 2.836187029},
	{6.040928121, 1.512938675, 19.853309201},
	{3.913974167, 7.760823775, 19.853309201},
	{3.607664892, 11.047534715, 14.180935144},
	{5.734618909, 4.799649429, 14.180935144},
	{5.428309508, 8.086360741, 8.508561086},
	{7.555263461, 1.838475642, 8.508561086},
	{7.24895425, 5.125186395, 2.836187029},
	{5.122000233, 11.373071681, 2.836187029},
	{11.598609508, 3.595567042, 19.853309201},
	{9.471655744, 9.843451583, 19.853309201},
	{13.419254155, 0.634392975, 14.180935144},
	{11.29230036, 6.882277609, 14.180935144},
	{10.985991085, 10.168988549, 8.508561086},
	{13.112944849, 3.921104008, 8.508561086},
	{12.806635701, 7.207814576, 2.836187029},
	{14.933589496, 0.959929942, 2.836187029}};

	vector<vector<double>> ones_PCLU_GROUP(PCLU_GROUP, vector<double>(1, 1.0));

	for (int i=0; i<=24; i++) {
		vector<vector<double>> center_position_block = mblock(center_position, i+1, 0+1, 1, 3);
		center_position_3D_matrix_T[i] = mmult(ones_PCLU_GROUP, center_position_block);
	}

	return;
}

void initialize_P8_and_P9_position_matrix_T() {
	P8_initial_position = {
	{0.0, 0.0, 0.0},
	{1.606897936, 0.602152331, -1.71601527},
    {1.606897936, 0.602152331, 1.71601527},
    {0.535374141, 1.706088546, 0.0},
    {0.494679473, -1.320096897, -1.1},
    {0.494679473, -1.320096897, 1.1},
    {-1.524733086, 0.934105248, 0.0},
    {-1.606897936, -0.602152331, -1.71601527},
    {-1.606897936, -0.602152331, 1.71601527},
	{0.0, 0.0, 0.0}, //Line 9
	{0.0, 0.0, 0.0}, //Line 10
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0}};

	P9_initial_position = {
	{0.0, 0.0, 0.0},
    {1.67254838, 0.850269373, -1.78124913},
    {1.67254838, 0.850269373, 1.78124913},
    {1.204686507, 2.511646196, 0.0},
    {0.411926137, -1.099262139, -1.14799694},
    {0.411926137, -1.099262139, 1.14799694},
    {-2.558612066, 1.101426569, 0.0},
    {-1.819437874, -0.458281478, -1.78124913},
    {-1.819437874, -0.458281478, 1.78124913},
    {-1.123831501, 2.999045969, 0.0},
    {-0.0734447470116576, 0.195993947619753, 0.0}, //Line 10 -- P_ortho
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0}};

	predict_P2_for_P8 = {{1.404330253, 0.730723945, 0.414282426}}; //a, b, c for linear combination of P 1(7), 2(8), 4+5
	predict_P2_for_P9 = {{1.325930882, 0.655279857, 0.50452553}}; 

	return;
}

void initialize_twist_rotation_matrix_T() {
	twist_rotation_matrix_T[0] = { //@20231123  0 was added since % operation will cause 0
	{1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0}};
	
	twist_rotation_matrix_T[40] = {
	{0.971192616624968, 0.0768751070430124, -0.225555135898682},
	{0.0768751070430124, 0.79485182649401, 0.601914272829829},
	{0.225555135898682, -0.601914272829829, 0.766044443118978}};

	twist_rotation_matrix_T[60] = {
	{0.938434068933695, 0.16429425329296, -0.303889612523128},
	{0.16429425329296, 0.561565931066305, 0.810956906007059},
	{0.303889612523128, -0.81095690600706, 0.5}};

	twist_rotation_matrix_T[80] = {
	{0.898249761339454, 0.271529711214977, -0.345570516944263},
	{0.271529711214977, 0.275398416327476, 0.922186167870582},
	{0.345570516944263, -0.922186167870582, 0.17364817766693}};

	twist_rotation_matrix_T[120] = {
	{0.815302206801086, 0.492882759878881, -0.303889612523128},
	{0.492882759878881, -0.315302206801086, 0.81095690600706},
	{0.303889612523128, -0.81095690600706, -0.5}};

	twist_rotation_matrix_T[160] = {
	{0.76116203563775, 0.637360701499772, -0.120015381045581},
	{0.637360701499772, -0.700854656423658, 0.320271895040753},
	{0.120015381045581, -0.320271895040753, -0.939692620785908}};

	twist_rotation_matrix_T[180] = {
	{0.753736275734781, 0.657177013171841, 0.0},
	{0.657177013171841, -0.753736275734781, 0.0},
	{0.0, 0.0, -1.0}};

	twist_rotation_matrix_T[200] = {
	{0.76116203563775, 0.637360701499772, 0.120015381045581},
	{0.637360701499772, -0.700854656423658, -0.320271895040753},
	{-0.120015381045581, 0.320271895040753, -0.939692620785908}};

	twist_rotation_matrix_T[240] = {
	{0.815302206801086, 0.492882759878881, 0.303889612523128},
	{0.492882759878881, -0.315302206801086, -0.810956906007059},
	{-0.303889612523128, 0.810956906007059, -0.5}};

	twist_rotation_matrix_T[280] = {
	{0.898249761339454, 0.271529711214977, 0.345570516944263},
	{0.271529711214977, 0.275398416327476, -0.922186167870582},
	{-0.345570516944263, 0.922186167870582, 0.17364817766693}};

	twist_rotation_matrix_T[300] = {
	{0.938434068933695, 0.16429425329296, 0.303889612523128},
	{0.16429425329296, 0.561565931066305, -0.810956906007059},
	{-0.303889612523128, 0.81095690600706, 0.5}};

	twist_rotation_matrix_T[320] = {
	{0.971192616624968, 0.0768751070430125, 0.225555135898682},
	{0.0768751070430124, 0.79485182649401, -0.601914272829829},
	{-0.225555135898682, 0.601914272829829, 0.766044443118978}};

	twist_rotation_matrix_T[360] = {
	{1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0}};

	// 20 deg series also exist. Don't forget.
	twist_rotation_matrix_T[20] = {
	{0.992574240097032, 0.0198163116720692, -0.120015381045581},
	{0.0198163116720692, 0.947118380688877, 0.320271895040753},
	{0.120015381045581, -0.320271895040753, 0.939692620785908}};

	twist_rotation_matrix_T[100] = {
	{0.855486514395327, 0.385647301956864, -0.345570516944263},
	{0.385647301956864, -0.0291346920622575, 0.922186167870582},
	{0.345570516944263, -0.922186167870582, -0.17364817766693}};

	twist_rotation_matrix_T[140] = {
	{0.782543659109813, 0.580301906128829, -0.225555135898682},
	{0.580301906128828, -0.548588102228791, 0.601914272829829},
	{0.225555135898682, -0.601914272829829, -0.766044443118978}};

	twist_rotation_matrix_T[220] = {
	{0.782543659109813, 0.580301906128829, 0.225555135898682},
	{0.580301906128828, -0.548588102228791, -0.601914272829829},
	{-0.225555135898682, 0.601914272829829, -0.766044443118978}};

	twist_rotation_matrix_T[260] = {
	{0.855486514395327, 0.385647301956864, 0.345570516944263},
	{0.385647301956864, -0.0291346920622576, -0.922186167870582},
	{-0.345570516944263, 0.922186167870582, -0.17364817766693}};

	twist_rotation_matrix_T[340] = {
	{0.992574240097032, 0.0198163116720692, 0.120015381045581},
	{0.0198163116720692, 0.947118380688877, -0.320271895040753},
	{-0.120015381045581, 0.320271895040753, 0.939692620785908}};

	return;
}

void initialize_wavy_rotation_matrix_T(double wavy_rotation_theta) {
	vector<vector<double>> R2(3, vector<double>(3)); 
	vector<vector<double>> Ry_theta(3, vector<double>(3)); 
	vector<vector<double>> R2u(3, vector<double>(3)); 

	R2 = {
	{0.936412376, -0.350901499, 0.0},
	{0.350901499, 0.936412376, 0.0},
	{0.0, 0.0, 1.0}};

	Ry_theta = {
	{cos(wavy_rotation_theta*pi/180), 0.0, sin(wavy_rotation_theta*pi/180)},
	{0.0, 1.0, 0.0},
	{-sin(wavy_rotation_theta*pi/180), 0.0, cos(wavy_rotation_theta*pi/180)}};

	R2u = {
	{0.936412376, 0.350901499, 0.0},
	{-0.350901499, 0.936412376, 0.0},
	{0.0, 0.0, 1.0}};

	vector<vector<double>> Rtmp = mmult(R2, Ry_theta);
	vector<vector<double>> wavy_rotation_matrix = mmult(Rtmp, R2u);

	wavy_rotation_matrix_T[wavy_rotation_theta] = mtranspose(wavy_rotation_matrix);

	return;
}

void initialize_wavy_rotation_matrix_T_of_slice1_and_3(double wavy_rotation_theta) {
	vector<vector<double>> R2(3, vector<double>(3)); 
	vector<vector<double>> Ry_theta1(3, vector<double>(3)); 
	vector<vector<double>> R2u(3, vector<double>(3)); 

	R2 = {
	{0.936412376, -0.350901499, 0.0},
	{0.350901499, 0.936412376, 0.0},
	{0.0, 0.0, 1.0}};

	Ry_theta1 = {
	{cos(wavy_rotation_theta*pi/180), 0.0, sin(wavy_rotation_theta*pi/180)},
	{0.0, 1.0, 0.0},
	{-sin(wavy_rotation_theta*pi/180), 0.0, cos(wavy_rotation_theta*pi/180)}};

	R2u = {
	{0.936412376, 0.350901499, 0.0},
	{-0.350901499, 0.936412376, 0.0},
	{0.0, 0.0, 1.0}};

	vector<vector<double>> Rtmp1 = mmult(R2, Ry_theta1);
	vector<vector<double>> wavy_rotation_matrix1 = mmult(Rtmp1, R2u);

	wavy_rotation_matrix1_T = mtranspose(wavy_rotation_matrix1);


	vector<vector<double>> Ry_theta3(3, vector<double>(3)); 
	Ry_theta3 = {
	{cos(-wavy_rotation_theta*pi/180), 0.0, sin(-wavy_rotation_theta*pi/180)},
	{0.0, 1.0, 0.0},
	{-sin(-wavy_rotation_theta*pi/180), 0.0, cos(-wavy_rotation_theta*pi/180)}};

	vector<vector<double>> Rtmp3 = mmult(R2, Ry_theta3);
	vector<vector<double>> wavy_rotation_matrix3 = mmult(Rtmp3, R2u);

	wavy_rotation_matrix3_T = mtranspose(wavy_rotation_matrix3);

	return;
}

void initialize_augmented_twist_rotation_matrix_T(double twist_rotation_theta) {
	vector<vector<double>> R2(3, vector<double>(3)); 
	vector<vector<double>> Rx_theta(3, vector<double>(3)); 
	vector<vector<double>> R2u(3, vector<double>(3)); 

	R2 = {
	{0.936412376, -0.350901499, 0.0},
	{0.350901499, 0.936412376, 0.0},
	{0.0, 0.0, 1.0}};

	Rx_theta = {
    {1.0, 0.0, 0.0},
	{0.0, cos(twist_rotation_theta*pi/180), -sin(twist_rotation_theta*pi/180)},
	{0.0, sin(twist_rotation_theta*pi/180), cos(twist_rotation_theta*pi/180)}};

	R2u = {
	{0.936412376, 0.350901499, 0.0},
	{-0.350901499, 0.936412376, 0.0},
	{0.0, 0.0, 1.0}};

	vector<vector<double>> Rtmp = mmult(R2, Rx_theta);
	vector<vector<double>> twist_rotation_matrix = mmult(Rtmp, R2u);

	augmented_twist_rotation_matrix_T[twist_rotation_theta] = mtranspose(twist_rotation_matrix);

	return;
}

void initialize_wavy_translation_matrix_T(double slice_wavy_shift) {
	vector<vector<double>> tmp_matrix(PCLU_GROUP);
	for (int i=0; i<PCLU_GROUP; i++) {
		tmp_matrix[i] = {0.0, 0.0, slice_wavy_shift};
	}
	wavy_translation_matrix_T[slice_wavy_shift] = tmp_matrix;
	return;
}

int neighbor_slice_cluster(int i0, int direction) { // -1 for LEFT, 1 for RIGHT.
	int i = -1;
	if (direction == LEFT) {
		i = (i0 - 8 + 24 ) % 24;
		if (i == 0) i = 24; //@20231123  convert 0 to 24
		if (i0 <= 8) {
			if ( i % 2 == 0) {
				i = i - 1;
			} else {
				i = i + 1;
			}
		}
	} else {
		i = (i0 + 8 + 24 ) % 24;
		if (i == 0) i = 24; //@20231123  convert 0 to 24
		if (i0 > 16) {
			if ( i % 2 == 0) {
				i = i - 1;
			} else {
				i = i + 1;
			}
		}
	}
	return i;
}

int determine_P8_twist_angle() {
	int this_model_number = 1;
	bool visitedP8[25] = {false};
	for (int i=1; i<=24; i++) {
		if (thisPCLU[i].N_atoms == 8 && visitedP8[i] == false) {
			int continuous_P8_cnt = 1; 
			int cls = neighbor_slice_cluster(i, LEFT); //cluster of LEFT slice 
			while (thisPCLU[cls].N_atoms == 8) {
				continuous_P8_cnt++; // @20231123 Adjusted Sequence with next line
				cls = neighbor_slice_cluster(cls, LEFT); // @20231123 should not be (i, LEFT)
				if ( continuous_P8_cnt >= 3 ) {
					std::cout << "Model Building Error: continuous_P8_cnt reaches 3.\n";
					break;
				}
			}
			int crs = neighbor_slice_cluster(i, RIGHT); //cluster of RIGHT slice 
			while (thisPCLU[crs].N_atoms == 8) {
				continuous_P8_cnt++;
				crs = neighbor_slice_cluster(crs, RIGHT);
				if ( continuous_P8_cnt >= 3 ) {
					std::cout << "Model Building Error: continuous_P8_cnt reaches 3.\n";
					break;
				}
			}

			int twist_difference = thisPCLU[crs].twist_angle - thisPCLU[cls].twist_angle;
			if ( twist_difference > 180 ) {
				twist_difference -= 360;
			} else if ( twist_difference <= -180 ) {
				twist_difference += 360;
			}

			//隔两层转180°就意味着有两种可能。
			if ( twist_difference == 180 && continuous_P8_cnt == 2) {
				this_model_number *= 2;
			}

			//现在±180都会变成+180也就是沿逆时针的正方向转60°每Slice
			int twist_per_slice = twist_difference / (continuous_P8_cnt + 1);
			int j=neighbor_slice_cluster(cls, RIGHT); 
			while (j != crs) {
				visitedP8[j] = true;
				thisPCLU[j].twist_angle = (thisPCLU[neighbor_slice_cluster(j, LEFT)].twist_angle + twist_per_slice + 360) % 360;
				j = neighbor_slice_cluster(j, RIGHT);
			}
		}
	}
	return this_model_number;
}

int determine_P8_delicate_twist_angle() {
	int this_model_number = 1;
	bool visitedP8[25] = {false};
	for (int i=1; i<=24; i++) {
		if (thisPCLU[i].N_atoms == 8 && visitedP8[i] == false) {
			int continuous_P8_cnt = 1; 
			int cls = neighbor_slice_cluster(i, LEFT); //cluster of LEFT slice 
			while (thisPCLU[cls].N_atoms == 8) {
				continuous_P8_cnt++; // @20231123 Adjusted Sequence with next line
				cls = neighbor_slice_cluster(cls, LEFT); // @20231123 should not be (i, LEFT)
				if ( continuous_P8_cnt >= 3 ) {
					std::cout << "Model Building Error: continuous_P8_cnt reaches 3.\n";
					break;
				}
			}
			int crs = neighbor_slice_cluster(i, RIGHT); //cluster of RIGHT slice 
			while (thisPCLU[crs].N_atoms == 8) {
				continuous_P8_cnt++;
				crs = neighbor_slice_cluster(crs, RIGHT);
				if ( continuous_P8_cnt >= 3 ) {
					std::cout << "Model Building Error: continuous_P8_cnt reaches 3.\n";
					break;
				}
			}

			double delicate_twist_difference = thisPCLU[crs].delicate_twist_angle - thisPCLU[cls].delicate_twist_angle;
			if ( delicate_twist_difference > 180 + eps ) {
				delicate_twist_difference -= 360;
			} else if ( delicate_twist_difference <= -180 + eps) {
				delicate_twist_difference += 360;
			}

			//隔两层转180°就意味着有两种可能。
			if ( fabs(delicate_twist_difference - 180) < eps && continuous_P8_cnt == 2) {
				this_model_number *= 2;
			}

			//现在±180都会变成+180也就是沿逆时针的正方向转60°每Slice
			double delicate_twist_per_slice = delicate_twist_difference / (continuous_P8_cnt + 1);
			int j=neighbor_slice_cluster(cls, RIGHT); 
			while (j != crs) {
				visitedP8[j] = true;
				double substitute_mod360 = (thisPCLU[neighbor_slice_cluster(j, LEFT)].delicate_twist_angle + delicate_twist_per_slice + 360);
				while (substitute_mod360 > 360.0 + eps) {
					substitute_mod360 -= 360.0;
				} 
				thisPCLU[j].delicate_twist_angle = substitute_mod360;
				j = neighbor_slice_cluster(j, RIGHT);
			}
		}
	}
	return this_model_number;
}

void fill_bond_list(StructuralInformation& structural_information, vector<vector<bool>>& bondGraph, int& bondIndexCnt, int overall_index_k, int overall_index_l) {
    structural_information.bondList[bondIndexCnt++] = {overall_index_k, overall_index_l}; 
    bondGraph[overall_index_k][overall_index_l] = 1;
    bondGraph[overall_index_l][overall_index_k] = 1;
}