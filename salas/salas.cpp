#include<stdlib.h>
#include<stdio.h>
#include<vector.h>

using namespace std;

typedef struct {
	int n,
	int* up,
	int* dn,
	int upn,
	int dnn,
	float fitness
} Solution;

float evaluate(Solution x, float** A, float* tam){
	float cost = 0;
	float* pos_x = new float[x.n];
	for(int i = 0; i < x.n; i++){
		int found = 0;
		float sum_x = 0;
		for(int j = 0; j < x.upn; j++){
			if(x.up[j] == i){
				pos_x[i] = sum_x+tam[i]*0.5;
				found = 1;
				break;
			}
			sum_x += tam[x.up[j]];
		}
		if(!found){
			sum_x = 0;
			for(int j = 0; j < x.dnn; j++){
				if(x.dn[j] == i){
					pos_x[i] = sum_x+tam[i]*0.5;
					found = 1;
					break;
				}
				sum_x += tam[x.dn[j]];
			}
		}
	}
	for(int i = 0; i < x.n; i++){
		for(int j = 0; j < x.n; j++){
			cost += fabs(pos_x[j]-pos_x[i])*A[i][j];
		}
	}
	return cost;
}

int main(){
	
}