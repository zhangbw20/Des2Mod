// File: Check_Form_II_Arrangement_In_Each_Slice_Pairs1-4.cpp
// Description: This file decides all ways of arrangements of P8, P9 blocks in Each Slice of 
//				Form II Red Phosphorus, according to Kim and Ma's Angew 2023.
// Author: Bowen ZHANG
// Date: 2023-09-05
// License: ThU License
// Dependencies: <bits/stdc++.h> 
// Usage: A typical input file is named as "Check_Form_II_Arrangement_In_Each_Slice_Input.txt" in
//		  the current directory. 
// Version: 
// Notes: We find it to be 545 ways in total.

#include <set>
#include <map>
#include <queue>
#include <stack>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

#define MAXN 9
int G[MAXN][MAXN]={0};
int N, Ne; 

bool visited[MAXN] = {false};
int cnt=0;
int Defined_depth;
vector<vector<int>> pairs;
vector<int> output_list(8);
int inSlice[545][9];

void dfs(int start_point, int depth) {
	if (depth == Defined_depth) {
		for (int i=0; i<depth*2; i++) {
			cout << output_list[i];
		}
		cout << endl;
		for (int k=0; k<depth; k++) {
			//讨论2*k和2*k+1 
			if (output_list[2*k+1] - output_list[2*k] == 2 
				|| output_list[2*k+1] - output_list[2*k] == -6) {
				inSlice[cnt][output_list[2*k]] = 4; //P9朝向为240° 
				inSlice[cnt][output_list[2*k+1]] = 1; //P9朝向为60° 
			} else if (output_list[2*k]%2==1) {//奇数 
				if (output_list[2*k+1] - output_list[2*k] == 1 
					|| output_list[2*k+1] - output_list[2*k] == -1) {
					inSlice[cnt][output_list[2*k]] = 6; //P9朝向为360° 
					inSlice[cnt][output_list[2*k+1]] = 3; //P9朝向为180° 
				}  
				if (0//output_list[2*k+1] - output_list[2*k] == 1
					|| output_list[2*k+1] - output_list[2*k] == 3
					|| output_list[2*k+1] - output_list[2*k] == -5
					|| output_list[2*k+1] - output_list[2*k] == -7) {
					inSlice[cnt][output_list[2*k]] = 5; //P9朝向为300°
					inSlice[cnt][output_list[2*k+1]] = 2; //P9朝向为120°
				}  
			} else { //偶数 
				if (0 //output_list[2*k+1] - output_list[2*k] == 1 
					|| output_list[2*k+1] - output_list[2*k] == -1) {
					inSlice[cnt][output_list[2*k]] = 6; //P9朝向为360° 
					inSlice[cnt][output_list[2*k+1]] = 3; //P9朝向为180° 
				}  
				if (output_list[2*k+1] - output_list[2*k] == 1
					|| output_list[2*k+1] - output_list[2*k] == 3
					|| output_list[2*k+1] - output_list[2*k] == -5
					|| output_list[2*k+1] - output_list[2*k] == -7) {
					inSlice[cnt][output_list[2*k]] = 5; //P9朝向为300°
					inSlice[cnt][output_list[2*k+1]] = 2; //P9朝向为120°
				}
			}
		}
		cnt++;
		return;
	}
	for (int i=start_point; i<Ne; i++) {
		if (visited[pairs[i][0]] == false && visited[pairs[i][1]] == false) {
			visited[pairs[i][0]] = true; 
			visited[pairs[i][1]] = true;
			output_list[depth*2] = pairs[i][0];
			output_list[depth*2+1] = pairs[i][1];

			dfs(i+1, depth+1);
			output_list[depth*2] = 0;
			output_list[depth*2+1] = 0;
			visited[pairs[i][0]] = false; 
			visited[pairs[i][1]] = false;
		}
	}
}

int main () {
	for (int i=0; i<545; i++) {
		for (int j=0; j<9; j++) {
			inSlice[i][j] = 8;
		}
	}
	
	cin >> N; cin >> Ne;
	pairs.resize(Ne, vector<int>(2));
	int tmp1, tmp2;
	for (int i=0; i<Ne; i++) {
		cin >> pairs[i][0] >> pairs[i][1];
	}
	//cin >> Defined_depth;
	
	//层数依赖的循环，只能dfs 
	
	for (Defined_depth=0; Defined_depth<=4; Defined_depth++) {
		dfs(0, 0);
	} 
	
	
	cout << cnt << endl << "---------------" << endl;
	
	for (int i=0; i<cnt; i++) {
		for (int j=1; j<=8; j++)
			cout << inSlice[i][j];
		cout << endl;
	}
	
	return 0;
}



