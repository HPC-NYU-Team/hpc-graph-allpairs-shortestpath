#include <iostream>
#include <climits> 
#include "graphGenerator.cpp"

template <class T>
void print(int order, T **A){
    for(int i =  0 ; i<order; i++){    
        for(int j = 0; j<order; j++)
            cout<<A[i][j]<<" "; 
        cout<<endl; 
    }
}

int main(){
    string fileName = "n5d2.random.edges"; 
    int order  = 5;  
    int **graph = getAdjacencyMatrixArray( fileName,5); 
    vector<vector<int> > neighbours = getAdjacencyListVector(fileName,5);

    for(int i = 0; i < neighbours.size(); i++){
        cout<<endl; 
        for(int j = 0; j < neighbours[i].size() ; j++)
            cout<<neighbours[i][j]<<" ";
    }
    
    cout<<endl; 
    bool **A = new bool*[order]; 
    for(int i = 0; i < order; i++)
        A[i] = new bool[order]; 

    bool **B = new bool*[order]; 
    for(int i = 0; i < order; i++)
        B[i] = new bool[order]; 

    int **ASPS = new int*[order]; 
    for(int i = 0; i < order; i++)
        ASPS[i] = new int[order]; 
    
    for(int i = 0; i < order; i++)
        for(int j = 0; j < order ; j++)
            if(j==order-1-i)
                ASPS[i][j] = 0; 
            else 
                ASPS[i][j] = INT_MAX; 


    
    int diameter = 1; 
    int distance = order * (order-1); 

    
   // cout<<"HI"; 
    //INITIALISING A AND B 
    for(int i = 0; i < order; i++)
        for(int j = 0; j < order; j++)
        {    
            if(j==order-1-i){
                A[i][j] = 1; 
                B[i][j] = 1;

            }
            else {
                A[i][j] = 0; 
                B[i][j] = 0; 
            }
        }

    int num; 
    for(int k = 0 ; k < order-1 ; k++){
        for(int i = 0; i < order ; i++){
            for(auto n: neighbours[i]){
                for(int j = 0 ; j < order ; j ++){
                    B[i][j] = B[i][j] || A[n][j]; 
                }

            }
        }
        num = 0; 
        for(int i = 0; i < order; i++)
            for(int j = 0; j < order ; j++){
                if(B[i][j]==1){
                    if(ASPS[i][j] > k)
                        ASPS[i][j] = k+1; 
                    num++;
                }
            }
    cout<<endl; 
    cout<<"B"<<endl; 
    print(5,B); 
    cout<<"ASPS"<<endl; 
    print(5,ASPS); 

    if(num == order*order) break; 
    swap(A,B); 
    diameter++; 
    distance = distance+(order * order - num); 

    }

    double average_distance = (double)distance/((order-1)* order); 
    cout<<diameter; 
    cout<<average_distance; 

 

    return 0; 
}





