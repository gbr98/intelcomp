#include<stdlib.h>
#include<stdio.h>
#include<vector.h>

using namespace std;

typedef struct {
	int n,
	int* up,
	int* dn,
	float fitness
} Solution;

float evaluate(Solution x, float** A, float* tam){
	float** dist = new float[x.n][x.n];
	for(int i = 0; i < x.n; i++){
		for(int j = 0; j < x.n; j++){
			dist[i][j] = 
		}
	}
}

int main(){

}