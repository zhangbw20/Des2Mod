// File: removeRepitition123.cpp
// Description: Via a combination of mirror, translation and leaning transformation, this file 
//				further removes repititive arrangements of P8, P9 blocks in Form II Red Phosphorus.
// Author: Bowen ZHANG
// Date: 2023-09-07
// License: ThU License
// Dependencies: <bits/stdc++.h> 
// Usage: A typical input file is named as "xxxx.txt" in
//		  the current directory. 
// Version: 
// Notes: 

#include <set>
#include <map>
#include <queue>
#include <stack>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#define MAXN 545
char G[MAXN][MAXN][MAXN]={0};
char inSlice8[MAXN][9];
int N, Ne; 
char tmpSlice[9][9];

//derivative global variables
int N_Check;
vector<vector<int>> checkList;
char checkSliceIJ[2][4][9];
map<char, char> charMap;
char checkSliceIJ_mtl[3][9][4][4][9]; //mtl: mirror + translation + leaning, slice, i
//@20240326 Original Order is Wrong... [9][4][3][4][9] should be [3][9][4][4][9]

vector<int> duplicate; //53376
map<string, int> des_2id;
vector<string> id2des_;

void fillcheckSliceIJ(int sliceId, int U, int V, int W) {
	for (int i=1; i<=8; i++) {
		checkSliceIJ[sliceId][1][i] = inSlice8[U][i];
		checkSliceIJ[sliceId][2][i] = inSlice8[V][i];
		checkSliceIJ[sliceId][3][i] = inSlice8[W][i];
	}
	return;
}

void writeAllMTL(int id) {
	//变换J与I相比较 
	//vector<vector<bool>> ifMTL(4, vector<bool>(9));
	int ifMTL_cnt = 0;
	
	char tmpChar[9] = {0}; 
	
	//第〇步先做E变换 
	for (int m=1; m<=2; m++) {
		for (int t=1; t<=8; t++) {
			for (int l=1; l<=3; l++) {
				for (int s=1; s<=3; s++) {
					for (int i=1; i<=8; i++) {
						checkSliceIJ_mtl[m][t][l][s][i] = checkSliceIJ[1][s][i]; 
					}
				}
			}
		}
	} 
	
	//第一步做lean变换 
	for (int m=1; m<=2; m++) {
		for (int t=1; t<=8; t++) {
			for (int s=1; s<=3; s++) {  
				//l=1: 无需再动 
				
				//l=2: 1-2互换, 5-6互换 需要tmpChar寄存！ 
				tmpChar[1] = charMap[checkSliceIJ_mtl[m][t][2][s][2]];
				tmpChar[2] = charMap[checkSliceIJ_mtl[m][t][2][s][1]];
				tmpChar[3] = charMap[checkSliceIJ_mtl[m][t][2][s][3]]; 
				tmpChar[4] = charMap[checkSliceIJ_mtl[m][t][2][s][4]];
				tmpChar[5] = charMap[checkSliceIJ_mtl[m][t][2][s][6]]; 
				tmpChar[6] = charMap[checkSliceIJ_mtl[m][t][2][s][5]]; 
				tmpChar[7] = charMap[checkSliceIJ_mtl[m][t][2][s][7]]; 
				tmpChar[8] = charMap[checkSliceIJ_mtl[m][t][2][s][8]];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][t][2][s][i] = tmpChar[i];
				}
				
				//l=3: 3-4互换, 7-8互换 
				tmpChar[1] = charMap[checkSliceIJ_mtl[m][t][3][s][1]];
				tmpChar[2] = charMap[checkSliceIJ_mtl[m][t][3][s][2]];
				tmpChar[3] = charMap[checkSliceIJ_mtl[m][t][3][s][4]]; 
				tmpChar[4] = charMap[checkSliceIJ_mtl[m][t][3][s][3]];
				tmpChar[5] = charMap[checkSliceIJ_mtl[m][t][3][s][5]]; 
				tmpChar[6] = charMap[checkSliceIJ_mtl[m][t][3][s][6]]; 
				tmpChar[7] = charMap[checkSliceIJ_mtl[m][t][3][s][8]]; 
				tmpChar[8] = charMap[checkSliceIJ_mtl[m][t][3][s][7]];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][t][3][s][i] = tmpChar[i];
				}
			}
		}
	}
	
	//第二步变换translation 
	for (int m=1; m<=2; m++) {
		for (int l=1; l<=3; l++) {
			for (int s=1; s<=3; s++) {
				//t=1: 无需再动

				//t=2
				tmpChar[1] = checkSliceIJ_mtl[m][2][l][s][2];
				tmpChar[2] = checkSliceIJ_mtl[m][2][l][s][1];
				tmpChar[3] = checkSliceIJ_mtl[m][2][l][s][4]; 
				tmpChar[4] = checkSliceIJ_mtl[m][2][l][s][3];
				tmpChar[5] = checkSliceIJ_mtl[m][2][l][s][6]; 
				tmpChar[6] = checkSliceIJ_mtl[m][2][l][s][5]; 
				tmpChar[7] = checkSliceIJ_mtl[m][2][l][s][8]; 
				tmpChar[8] = checkSliceIJ_mtl[m][2][l][s][7];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][2][l][s][i] = tmpChar[i];
				}
				
				//t=3
				tmpChar[1] = checkSliceIJ_mtl[m][3][l][s][3];
				tmpChar[2] = checkSliceIJ_mtl[m][3][l][s][4];
				tmpChar[3] = checkSliceIJ_mtl[m][3][l][s][5]; 
				tmpChar[4] = checkSliceIJ_mtl[m][3][l][s][6];
				tmpChar[5] = checkSliceIJ_mtl[m][3][l][s][7]; 
				tmpChar[6] = checkSliceIJ_mtl[m][3][l][s][8]; 
				tmpChar[7] = checkSliceIJ_mtl[m][3][l][s][1]; 
				tmpChar[8] = checkSliceIJ_mtl[m][3][l][s][2];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][3][l][s][i] = tmpChar[i];
				}
				
				//t=4
				tmpChar[1] = checkSliceIJ_mtl[m][4][l][s][4];
				tmpChar[2] = checkSliceIJ_mtl[m][4][l][s][3];
				tmpChar[3] = checkSliceIJ_mtl[m][4][l][s][6]; 
				tmpChar[4] = checkSliceIJ_mtl[m][4][l][s][5];
				tmpChar[5] = checkSliceIJ_mtl[m][4][l][s][8]; 
				tmpChar[6] = checkSliceIJ_mtl[m][4][l][s][7]; 
				tmpChar[7] = checkSliceIJ_mtl[m][4][l][s][2]; 
				tmpChar[8] = checkSliceIJ_mtl[m][4][l][s][1];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][4][l][s][i] = tmpChar[i];
				}
				
				//t=5
				tmpChar[1] = checkSliceIJ_mtl[m][5][l][s][5];
				tmpChar[2] = checkSliceIJ_mtl[m][5][l][s][6];
				tmpChar[3] = checkSliceIJ_mtl[m][5][l][s][7]; 
				tmpChar[4] = checkSliceIJ_mtl[m][5][l][s][8];
				tmpChar[5] = checkSliceIJ_mtl[m][5][l][s][1]; 
				tmpChar[6] = checkSliceIJ_mtl[m][5][l][s][2]; 
				tmpChar[7] = checkSliceIJ_mtl[m][5][l][s][3]; 
				tmpChar[8] = checkSliceIJ_mtl[m][5][l][s][4];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][5][l][s][i] = tmpChar[i];
				}
				
				//t=6
				tmpChar[1] = checkSliceIJ_mtl[m][6][l][s][6];
				tmpChar[2] = checkSliceIJ_mtl[m][6][l][s][5];
				tmpChar[3] = checkSliceIJ_mtl[m][6][l][s][8]; 
				tmpChar[4] = checkSliceIJ_mtl[m][6][l][s][7];
				tmpChar[5] = checkSliceIJ_mtl[m][6][l][s][2]; 
				tmpChar[6] = checkSliceIJ_mtl[m][6][l][s][1]; 
				tmpChar[7] = checkSliceIJ_mtl[m][6][l][s][4]; 
				tmpChar[8] = checkSliceIJ_mtl[m][6][l][s][3];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][6][l][s][i] = tmpChar[i];
				}
				
				//t=7
				tmpChar[1] = checkSliceIJ_mtl[m][7][l][s][7];
				tmpChar[2] = checkSliceIJ_mtl[m][7][l][s][8];
				tmpChar[3] = checkSliceIJ_mtl[m][7][l][s][1]; 
				tmpChar[4] = checkSliceIJ_mtl[m][7][l][s][2];
				tmpChar[5] = checkSliceIJ_mtl[m][7][l][s][3]; 
				tmpChar[6] = checkSliceIJ_mtl[m][7][l][s][4]; 
				tmpChar[7] = checkSliceIJ_mtl[m][7][l][s][5]; 
				tmpChar[8] = checkSliceIJ_mtl[m][7][l][s][6];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][7][l][s][i] = tmpChar[i];
				}
				
				//t=8
				tmpChar[1] = checkSliceIJ_mtl[m][8][l][s][8];
				tmpChar[2] = checkSliceIJ_mtl[m][8][l][s][7];
				tmpChar[3] = checkSliceIJ_mtl[m][8][l][s][2]; 
				tmpChar[4] = checkSliceIJ_mtl[m][8][l][s][1];
				tmpChar[5] = checkSliceIJ_mtl[m][8][l][s][4]; 
				tmpChar[6] = checkSliceIJ_mtl[m][8][l][s][3]; 
				tmpChar[7] = checkSliceIJ_mtl[m][8][l][s][6]; 
				tmpChar[8] = checkSliceIJ_mtl[m][8][l][s][5];
				for (int i=1; i<=8; i++) {
					checkSliceIJ_mtl[m][8][l][s][i] = tmpChar[i];
				}
			}
		}
	}
	
	//第三步变换mirror 
	for (int t=1; t<=8; t++) {
		for (int l=1; l<=3; l++) {
			//m=1无需再动
			
			//m=2 
			//Slice 2的值不动，Slice1,3交换 
			for (int i=1; i<=8; i++) { //记录Slice 1的值，即将被覆盖 
				tmpChar[i] = checkSliceIJ_mtl[2][t][l][1][i];
			}
			for (int i=1; i<=8; i++) { //将Slice 3的值赋给Slice 1 
				checkSliceIJ_mtl[2][t][l][1][i] = checkSliceIJ_mtl[2][t][l][3][i];
			}
			for (int i=1; i<=8; i++) { //将tmpChar的值赋给Slice 3 
				checkSliceIJ_mtl[2][t][l][3][i] = tmpChar[i];
			}
			
		}
	}
	//cout << "debug: m t l ifMTL_cnt\n";
	for (int m=1; m<=2; m++) {
		for (int t=1; t<=8; t++) {
			for (int l=1; l<=3; l++) {
				//ifMTL_cnt = 0; //置零 ifMTL 
				if (m==1 && t==1 && l==1) continue;
				string tmp_des_ = "";

				for (int s=1; s<=3; s++) {
					for (int i=1; i<=8; i++) {
						//if (checkSliceIJ_mtl[m][t][l][s][i] == checkSliceIJ[0][s][i]) {
						//	ifMTL_cnt++;
						//}
						tmp_des_ = tmp_des_ + checkSliceIJ_mtl[m][t][l][s][i];
						//cout << checkSliceIJ_mtl[m][t][l][s][i];
					}
					tmp_des_ = tmp_des_ + "_";
					//cout << "_";
				}
				//cout << "debug: " << tmp_des_ << "\n";
				//cout << "  ";
				//cout << "debug: " << m << " " << t << " " << l << " " << ifMTL_cnt << "\n";
				auto it = des_2id.find(tmp_des_);
				if (it != des_2id.end() && it->second > id) {
					duplicate[it->second] = id;
				}
				//if (ifMTL_cnt == 24) return true;
			}
		}
	} 
	return;
	//return false;
}

int main() {
	charMap['1'] = '2'; charMap['2'] = '1'; 
	charMap['4'] = '5'; charMap['5'] = '4'; 
	charMap['3'] = '6'; charMap['6'] = '3'; 
	charMap['8'] = '8'; 
	
	ifstream myfileIN;
	myfileIN.open("60deg_removeRepetition_Input.txt");
	myfileIN >> N;
	string s_tmp;
	for (int i=0; i<N; i++) {
		myfileIN >> s_tmp; 
		for (int j=0; j<8; j++) {
			inSlice8[i][j+1] = s_tmp[j];
		}
	}
	myfileIN >> N_Check;
	checkList.resize(N_Check, vector<int>(4));
	duplicate.resize(N_Check);
	for (int i=0; i<N_Check; i++) {
		duplicate[i] = -1;
	}
	id2des_.resize(N_Check);
	for (int i=0; i<N_Check; i++) {
		checkList[i][0] = -1; //默认为通过 
		string tmp_des_ = "";
		for (int s=1; s<=3; s++) {
			myfileIN >> checkList[i][s];
			for (char c=1; c<=8; c++) {
				tmp_des_ = tmp_des_ + inSlice8[checkList[i][s]][c];
			}
			tmp_des_ = tmp_des_ + "_";
		}
		id2des_[i] = tmp_des_;
		des_2id[tmp_des_] = i;
	}
	myfileIN.close();
	cout << "Progress:\n";
	for (int i=0; i<N_Check; i++) {
		if ((i+1)%100 == 0) {
			cout << "  " << i+1 << " / " << N_Check << " ...\n";
		}
		//for (int j=i+1; j<N_Check; j++) {
			//fillcheckSliceIJ(0, checkList[i][1], checkList[i][2], checkList[i][3]);
			//fillcheckSliceIJ(1, checkList[j][1], checkList[j][2], checkList[j][3]);
			fillcheckSliceIJ(1, checkList[i][1], checkList[i][2], checkList[i][3]);
			writeAllMTL(i);
			//if (checkIfMTL() == true) {
			//	checkList[j][0] = i;
			//}
		//}
	} 
	
	ofstream myfileOUT;
	myfileOUT.open("60deg_removeRepetition_123_Output.txt");
	for (int i=0; i<N_Check; i++) {
		myfileOUT << duplicate[i] << " " << checkList[i][1] << " " << checkList[i][2] << " " << checkList[i][3] <<  endl;
	}
	myfileOUT.close();
	
	return 0;
}
