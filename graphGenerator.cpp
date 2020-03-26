#include <iostream>
#include<string> 
#include <fstream> 
#include <sstream> 
#include <vector> 

using namespace std; 



int ** getAdjacencyMatrixArray(string fileName, int order){

    int **graph = new int*[order]; 
    for(int i = 0; i < order; i++)
        graph[i] = new int[order]; 

    for(int i=0; i<order; i++)
        for(int j=0; j<order; j++)
        graph[i][j] = 0; 

    ifstream graphFile; 
    graphFile.open(fileName);
    string buffer,word;
    int pid = 0;
    int u , v; 
    while(getline(graphFile,buffer)){
        stringstream ss(buffer); 
        ss>>word; 
        u = stoi(word); 
        ss>>word; 
        v = stoi(word); 
        graph[u][v] = 1; 
        graph[v][u] = 1; 

    }
    graphFile.close(); 
    return graph; 
}

vector<vector<int> > getAdjacencyMatrixVector(string fileName, int order){
    vector<vector<int> > graph(order, vector<int>(order, 0));
    ifstream graphFile; 
    graphFile.open(fileName);
    string buffer,word;
    int pid = 0;
    int u , v; 
    while(getline(graphFile,buffer)){
        stringstream ss(buffer); 
        ss>>word; 
        u = stoi(word); 
        ss>>word; 
        v = stoi(word); 
        graph[u][v] = 1; 
        graph[v][u] = 1; 

    }
    graphFile.close(); 
    return graph; 
}
vector<vector<int> > getAdjacencyListVector(string fileName, int order){
    vector<vector<int> > graph(order); 
    ifstream graphFile; 
    graphFile.open(fileName);
    string buffer,word;
    int pid = 0;
    int u , v; 
    while(getline(graphFile,buffer)){
        stringstream ss(buffer); 
        ss>>word; 
        u = stoi(word); 
        ss>>word; 
        v = stoi(word); 
        graph[u].push_back(v); 
        graph[v].push_back(u); 

    }
    graphFile.close(); 
    return graph; 
}


    //TO FREE:
      // free
    // for(int i = 0; i < order; ++i)
    //     delete [] graph[i];
    // delete [] graph;