#include <iostream>
#include <omp.h>
#include <algorithm> 
#include "graphGenerator.cpp"
#include <numeric>

using namespace std;

void getAllChildren(vector<vector<int>> &g, vector<int> &frontier, vector<int> &children, int* distances, int k){
	for(int i=0; i<frontier.size(); i++){
		int node = frontier[i];
		for(int j=0; j<g[node].size(); j++){
			int n = g[node][j];
			if(distances[n] == -1){
				distances[n] = k;
				children.push_back(n);
			}
		}
	}
}

int* bfs(vector<vector<int>> g, int source){
	vector<int> frontier = {source};
	vector<int> children = vector<int>{};
	
	int *distances = new int[g.size()];
	for(int i=0; i<g.size(); i++){
		distances[i] = -1;
	}
	
	int k = 1;
	while(frontier.size() > 0){
		getAllChildren(g, frontier, children, distances, k++);
		frontier = children;
		children = {};
	}
	return distances;
}


void serialBfsApsp(vector<vector<int>> g){
	int source = 0;
	float* diameter = new float[g.size()];
	float* distance = new float[g.size()];
	
	for(int source=0; source<g.size(); source++){
		int* distances = bfs(g, source);	
		diameter[source] = *max_element(distances, distances + g.size());
		distance[source] = accumulate(distances, distances + g.size(), 0.0) / g.size();
	}
	cout << "Max diameter: " << *max_element(diameter, diameter + g.size()) << endl;
	cout << "Avg distance: " << (accumulate(distance, distance + g.size(), 0.0) / g.size()) << endl;
}

void parallelBfsApsp(vector<vector<int>> g){
	
	int source = 0;
	float* diameter = new float[g.size()];
	float* distance = new float[g.size()];
	
	#pragma omp parallel for num_threads(8)
	for(int source=0; source<g.size(); source++){
		int* distances = bfs(g, source);	
		diameter[source] = *max_element(distances, distances + g.size());
		distance[source] = accumulate(distances, distances + g.size(), 0.0) / g.size();
	}
	cout << "Max diameter: " << *max_element(diameter, diameter + g.size()) << endl;
	cout << "Avg distance: " << (accumulate(distance, distance + g.size(), 0.0) / g.size()) << endl;
}



int main(int argc, char *argv[]) {
	
	string fileName = "n1024d4.random.edges"; 
	int order = 1024;
	
	vector<vector<int>> g = getAdjacencyListVector(fileName,order); 

	double tt = omp_get_wtime();
	serialBfsApsp(g);
	printf("sequential-bfs-apsp = %fs\n", omp_get_wtime() - tt);

	tt = omp_get_wtime();
	parallelBfsApsp(g);
	printf("parallel-bfs-apsp = %fs\n", omp_get_wtime() - tt);
	
}