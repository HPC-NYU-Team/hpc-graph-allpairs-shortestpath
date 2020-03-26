#include <iostream>
#include "graphGenerator.cpp"

using namespace std;
int main(int argc, char *argv[]) {
	string fileName = "n5d2.random.edges"; 
	int order = 5;
	
	cout<<"Adjacency matrix using arrays:\n"; 
	int **graph = getAdjacencyMatrixArray( fileName,order); 
	for(int i=0; i<order; i++)
	{    
		for(int j=0; j<order; j++)
			cout<<graph[i][j]<<" "; 
		cout<<endl;
	}


	vector<vector<int> > g; 

	cout<<"Adjacency matrix using vectors:\n"; 
	g = getAdjacencyMatrixVector(fileName,order); 
	for(int i=0; i<order; i++)
	{    
		for(int j=0; j<order; j++)
			cout<<g[i][j]<<" "; 
		cout<<endl;
	}

	cout<<"Adjacency list using vectors:\n"; 
	g = getAdjacencyListVector(fileName,order); 
	
	for(int i=0; i<order; i++)
	{    cout<<i<<" => "; 
		for(int j=0; j<g[i].size(); j++)
			cout<<g[i][j]<<" "; 
		cout<<endl;
	}
	return 0; 

}