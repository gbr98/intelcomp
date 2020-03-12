#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<vector>
#include<string>
#include<iostream>
#include<random>

using namespace std;

typedef struct {
	int n;
	int* up;
	int* dn;
	int upn;
	int dnn;
	float fitness;
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
					break;
				}
				sum_x += tam[x.dn[j]];
			}
		}
	}
	for(int i = 0; i < x.n; i++){
		for(int j = i; j < x.n; j++){
			cost += fabs(pos_x[j]-pos_x[i])*A[j][i];
		}
	}
	delete pos_x;
	return cost;
}

void physical_random(double** B, double** A, double* tam, int n){
	double* x = new double[n];
	double* v = new double[n];
	for(int i = 0; i < n; i++){
		x[i] = unif(0.0,1.0);
	}
	for(int it = 0; it < 100; it++){
		
	}
}

void string_to_vec(string strc, float* vec){
	string str = strc;
	int curr = str.find(","), prev = 0;
	int i;
	for(i = 0; curr != string::npos; i++){
		vec[i] = atof(str.substr(prev, curr-prev).c_str());
		prev = curr+1;
		curr = str.find(",", prev);
	}
	vec[i] = atof(str.substr(prev, curr-prev).c_str());
}

int main(int argc, char** argv){
	FILE* ins = NULL;
	if(argc > 1){
		ins = fopen(argv[1], "r");
	} else {
		cout << "Argumentos insuficientes. Terminando." << endl;
		exit(1);
	}
	int n;
	fscanf(ins, "%d", &n);
	float* tam = new float[n];
	float** A = new float*[n];
	for(int i = 0; i < n; i++){
		A[i] = new float[n];
	}
	char* strc = new char[4*n];
	fscanf(ins, "%s", strc);
	string_to_vec(strc, tam);
	for(int i = 0; i < n; i++){
		fscanf(ins, "%s", strc);
		string_to_vec(strc, A[i]);
	}

	Solution s = {n, new int[5]{5,1,9,0,2}, new int[5]{7,3,4,6,8}, 5, 5};
	s.fitness = evaluate(s, A, tam);
	//cout << s.fitness << endl;

	// construct
	int* disp = new int[n];
	for(int i = 0; i < n; i++){
		disp[i] = i;
	}
	vector<int> to_place(disp, disp+n);
	vector<int> up, dn;
	int pick = 0, side;
	while(!to_place.empty()){
		pick = rand()%to_place.size();
		side = rand()%2;
		if(side == 1){
			up.push_back(to_place[pick]);
		} else {
			dn.push_back(to_place[pick]);
		}
		to_place.erase(to_place.begin()+pick);
	}
	Solution sol = {n, &up[0], &dn[0], (int)up.size(), (int)dn.size()};
	sol.fitness = evaluate(sol, A, tam);
	cout << sol.fitness << endl;
	
	// deleting...
	for(int i = 0; i < n; i++){
		delete A[i];
	}
	delete A;
	delete tam;
	delete strc;
	delete disp;
	return 0;
}