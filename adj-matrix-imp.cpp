#include <iostream>
#include <climits> 
#include <omp.h>
#include "graphGenerator.cpp"

template <class T>
void print(int order, T **A){
    for(int i =  0 ; i<order; i++){    
        for(int j = 0; j<order; j++)
            cout<<A[i][j]<<" "; 
        cout<<endl; 
    }
}

void initialiseAB(bool **A, bool **B, int order){
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

}

void initialiseAPSP(int **APSP, int order){
    for(int i = 0; i < order; i++)
        for(int j = 0; j < order ; j++)
            if(j==order-1-i)
                APSP[i][j] = 0; 
            else 
                APSP[i][j] = INT_MAX; 
}


pair<int, double> serailAdjAPSP(int **graph, int order, vector<vector<int> > &neighbours, int **APSP){
    bool **A = new bool*[order]; 
    for(int i = 0; i < order; i++)
        A[i] = new bool[order]; 

    bool **B = new bool*[order]; 
    for(int i = 0; i < order; i++)
        B[i] = new bool[order]; 

    initialiseAB(A,B,order);
    initialiseAPSP(APSP,order); 

    int diameter = 1; 
    int distance = order * (order-1); 
    double average_distance = 0.0; 

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
                    if(APSP[i][j] > k)
                        APSP[i][j] = k+1; 
                    num++;
                }
            }
        if(num == order*order) break; 
        swap(A,B); 
        diameter++; 
        distance = distance+(order * order - num); 

    }

    average_distance = (double)distance/((order-1)* order); 
    
    return make_pair(diameter,average_distance); 
}

int main(int argc, char *argv[]){
    string fileName = "n10d4.random.edges"; 
    int order  = 10;  
    int **graph = getAdjacencyMatrixArray( fileName,order); 
     
    int **APSP = new int*[order]; 
    for(int i = 0; i < order; i++)
        APSP[i] = new int[order];  


    vector<vector<int> > neighbours = getAdjacencyListVector(fileName,order);
    double tt = omp_get_wtime();
	pair<int, double> diamAvgDist = serailAdjAPSP(graph, order, neighbours, APSP);
	printf("sequential-adj-apsp = %fs\n", omp_get_wtime() - tt);
    

    return 0; 
}





