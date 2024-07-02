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

void fillTmpSlice(int U, int V, int W) {
	for (int i=1; i<=8; i++) {
		tmpSlice[1][i] = inSlice8[U][i];
		tmpSlice[2][i] = inSlice8[V][i];
		tmpSlice[3][i] = inSlice8[W][i];
		if (i%2==1) {//奇数
			tmpSlice[4][i] = inSlice8[U][i+1];
			tmpSlice[5][i] = inSlice8[V][i+1];
			tmpSlice[6][i] = inSlice8[W][i+1];
		} else {//偶数 
			tmpSlice[4][i] = inSlice8[U][i-1];
			tmpSlice[5][i] = inSlice8[V][i-1];
			tmpSlice[6][i] = inSlice8[W][i-1];
		}
		tmpSlice[7][i] = inSlice8[U][i];
		tmpSlice[8][i] = inSlice8[V][i];
	}
	return;
}

int checkTwist() {
	int twistExtent = 0, twistExtent2 = 0, twistExtent3 = 0;
	int tmp;
	for (int s=1; s<=7; s++) { //Two adjacent slices
		for (int i=1; i<=8; i++) {
			if (tmpSlice[s][i] != '8' && tmpSlice[s+1][i] != '8') {
				tmp = tmpSlice[s][i] - tmpSlice[s+1][i];
				if (tmp < 0) tmp = -tmp;
				if (6-tmp < tmp) tmp = 6-tmp;
				if ( tmp >= 2) {
					return 8;
				} 
				if (tmp > twistExtent) {
					twistExtent = tmp;
				}
			}
		}		
	} 
	for (int s=1; s<=6; s++) { //1st and 3rd slice (次近邻的 second nearest slices)
		for (int i=1; i<=8; i++) {
			if ( tmpSlice[s+1][i] == '8' 
				&& tmpSlice[s][i] != '8' && tmpSlice[s+2][i] != '8') {
				tmp = tmpSlice[s][i] - tmpSlice[s+2][i];
				if (tmp < 0) tmp = -tmp;
				if (6-tmp < tmp) tmp = 6-tmp;
				if ( tmp >= 3) {
					return 9; 
				} 
				if (tmp > twistExtent2) {
					twistExtent2 = tmp;
				}
			}
		}		
	}
	for (int s=1; s<=5; s++) { //1st and 4th slice (次次近邻的 third nearest slices)
		for (int i=1; i<=8; i++) {
			if ( tmpSlice[s+1][i] == '8' && tmpSlice[s+2][i] == '8'
				&& tmpSlice[s][i] != '8' && tmpSlice[s+3][i] != '8') {
				tmp = tmpSlice[s][i] - tmpSlice[s+3][i];
				if (tmp < 0) tmp = -tmp;
				if (6-tmp < tmp) tmp = 6-tmp;
				if (tmp > twistExtent3) {
					twistExtent3 = tmp; //180°靠3层转过来 
				}
			}
		}		
	}
	//if (twistExtent3 == 3) {
	//	return 3;
	//} else if (2*twistExtent > twistExtent2) {
	//	return 2*twistExtent;
	//} else return twistExtent2;
	twistExtent *= 6;
	twistExtent2*= 3;
	twistExtent3*= 2;
	int result = twistExtent;
	if (twistExtent2 > result) result = twistExtent2;
	if (twistExtent3 > result) result = twistExtent3;
	return result;
}

int checkHawthorn() {
	for (int s=1; s<=6; s++) {
		for (int i=1; i<=8; i++) {
			if (tmpSlice[s][i] == '8' && tmpSlice[s+1][i] == '8' && tmpSlice[s+2][i] == '8') {
				return 49; //'a'-'0'
			} else if (tmpSlice[s][i] != '8' && tmpSlice[s+1][i] != '8' && tmpSlice[s+2][i] != '8' ){
				return 17; //'A'-'0'
			}
		}
	}
	return 0;
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
	ofstream myfileOUT4;
	ofstream myfileOUT6;
	myfileOUT4.open("checkTwist4_output_12.txt");
	myfileOUT6.open("checkTwist6_output_53376.txt");
	
	int cnt0=0, cnt2=0, cnt3=0, cnt4=0, cnt6=0;
	for (int i=0; i<N; i++) {
		//cout << "\n\n----  i = " << i << "  ----\n";
		for (int j=0; j<N; j++){
			for (int k=0; k<N; k++) {
				fillTmpSlice(i, j, k);
				G[i][j][k] = checkTwist()+'0';
				G[i][j][k] += checkHawthorn();
				//cout << " " << G[i][j][k];
				if (G[i][j][k] == '0' ) {
					cnt0++; 
					//cout << G[i][j][k] << " " << i << " " << j << " " << k; //向屏幕输出Twist'0'的结果 
					//for (int s=1; s<=3; s++) {
						//cout << " "; //向屏幕输出Twist'0'的结果 
						//for (int t=1; t<=8; t++) {
							//cout << tmpSlice[s][t]; //向屏幕输出Twist'0'的结果 
						//}
					//}
					//cout << endl; //向屏幕输出Twist'0'的结果 
				}
				if (G[i][j][k] == '2') cnt2++; 
				if (G[i][j][k] == '3') cnt3++; 
				if (G[i][j][k] == '4') {
					cnt4++; 
					myfileOUT4 << G[i][j][k] << " " << i << " " << j << " " << k; //向文件输出Twist'4'的结果 
					for (int s=1; s<=3; s++) {
						myfileOUT4 << " "; 
						for (int t=1; t<=8; t++) {
							myfileOUT4 << tmpSlice[s][t]; 
						}
					}
					myfileOUT4 << endl; 
				}
				if (G[i][j][k] == '6') {
					cnt6++; 
					myfileOUT6 << G[i][j][k] << " " << i << " " << j << " " << k; //向文件输出Twist'4'的结果 
					for (int s=1; s<=3; s++) {
						myfileOUT6 << " "; 
						for (int t=1; t<=8; t++) {
							myfileOUT6 << tmpSlice[s][t]; 
						}
					}
					myfileOUT6 << endl; 
				}
			}
			//cout << endl;
		}
		//cout << endl;
	}
	
	
	cout << "cnt0: Total between all slices: " << cnt0 << endl;
	cout << "cnt2: Total between all slices: " << cnt2 << endl;
	cout << "cnt3: Total between all slices: " << cnt3 << endl;
	cout << "cnt4: Total between all slices: " << cnt4 << endl;
	cout << "cnt6: Total between all slices: " << cnt6 << endl;
	
	myfileOUT4.close();
	myfileOUT6.close();
	
	return 0;
}
