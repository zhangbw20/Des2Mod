#include <set>
#include <map>
#include <queue>
#include <stack>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
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

void fillcheckSliceIJ(int sliceId, int U, int V, int W) {
	for (int i=1; i<=8; i++) {
		checkSliceIJ[sliceId][1][i] = inSlice8[U][i];
		checkSliceIJ[sliceId][2][i] = inSlice8[V][i];
		checkSliceIJ[sliceId][3][i] = inSlice8[W][i];
	}
	return;
}

bool checkIfTranslate() {
	vector<bool> ifTranslate(4);
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// i==2 以2为起点 
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][2]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][1]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][4]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][3]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][6]
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][5]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][7]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;
	
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// 以3为起点 
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][3]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][4]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][5]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][6]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][7]
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][1]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][2]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;
	
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// 以4为起点
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][4]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][3]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][6]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][5]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][7]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][2]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][1]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;
	
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// 以5为起点	
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][5]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][6]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][7]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][1]
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][2]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][3]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][4]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;
	
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// 以6为起点
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][6]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][5]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][7]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][2] 
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][1]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][4]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][3]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;
	
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// 以7为起点	
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][7]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][1]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][2]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][3]
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][4]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][5]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][6]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;
	
	ifTranslate = {false, false, false, false};
	for (int s=1; s<=3; s++) {
		// 以8为起点
		if (   checkSliceIJ[1][s][1] == checkSliceIJ[0][s][8]
			&& checkSliceIJ[1][s][2] == checkSliceIJ[0][s][7]
			&& checkSliceIJ[1][s][3] == checkSliceIJ[0][s][2]
			&& checkSliceIJ[1][s][4] == checkSliceIJ[0][s][1]
			&& checkSliceIJ[1][s][5] == checkSliceIJ[0][s][4] 
			&& checkSliceIJ[1][s][6] == checkSliceIJ[0][s][3]
			&& checkSliceIJ[1][s][7] == checkSliceIJ[0][s][6]
			&& checkSliceIJ[1][s][8] == checkSliceIJ[0][s][5]) {
			ifTranslate[s] = true;
		}
	}
	if ( ifTranslate[1] && ifTranslate[2] && ifTranslate[3] ) return true;

	return false;
}


int main() {
	cin >> N;
	string s_tmp;
	for (int i=0; i<N; i++) {
		cin >> s_tmp; 
		for (int j=0; j<8; j++) {
			inSlice8[i][j+1] = s_tmp[j];
		}
	}
	cin >> N_Check;
	checkList.resize(N_Check, vector<int>(4));
	for (int i=0; i<N_Check; i++) {
		checkList[i][0] = -1; //默认为通过 
		for (int s=1; s<=3; s++) {
			cin >> checkList[i][s];
		}
	}
	for (int i=0; i<N_Check; i++) {
		for (int j=i+1; j<N_Check; j++) {
			fillcheckSliceIJ(0, checkList[i][1], checkList[i][2], checkList[i][3]);
			fillcheckSliceIJ(1, checkList[j][1], checkList[j][2], checkList[j][3]);
			if (checkIfTranslate() == true) {
				checkList[j][0] = i;
			}
		}
	} 
	
	for (int i=0; i<N_Check; i++) {
		cout << checkList[i][0] << " " << checkList[i][1] << " " << checkList[i][2] << " " << checkList[i][3] <<  endl;
	}
	
	return 0;
}
