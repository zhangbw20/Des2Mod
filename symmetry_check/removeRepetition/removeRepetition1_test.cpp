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
			if ( //checkList[i][0] == -1 //@20240325 There should be no -1 check. Use Disjoint Set Union.
				   checkList[j][2] == checkList[i][2]
				&& checkList[j][1] == checkList[i][3]
				&& checkList[j][3] == checkList[i][1]) {
				checkList[j][0] = i;
			}
		}
	} 
	
	for (int i=0; i<N_Check; i++) {
		cout << checkList[i][0] << " " << checkList[i][1] << " " << checkList[i][2] << " " << checkList[i][3] <<  endl;
	}
	
	return 0;
}
