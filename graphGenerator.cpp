#include <iostream>
#include<string> 
#include <fstream> 
#include <sstream> 
#include <vector> 

using namespace std; 



int ** getAdjacencyMatrixArray(  string fileName, int order){

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

int main(){
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



    //TO FREE:
      // free
    // for(int i = 0; i < order; ++i)
    //     delete [] graph[i];
    // delete [] graph;