#include <iostream>
#include <climits> 
#include <omp.h>
#include <cmath>
#include "graphGenerator.cpp"
#define CHUNK_SIZE 100
using namespace std; 
template <class T>
void print(int order, T **A){
    for(int i =  0 ; i<order; i++){    
        for(int j = 0; j<order; j++)
            cout<<A[i][j]<<" "; 
        cout<<endl; 
    }
}

void initialiseAB(bool **A, bool **B, int col_beg, int col_end,int order){
    for(int i = 0; i < order; i++)
        for(int j = col_beg; j < col_end; j++)
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

void initialiseAPSP(int **APSP,int col_beg, int col_end,int order){
    for(int i = 0; i < order; i++)
        for(int j = col_beg; j < col_end ; j++)
            if(j==order-1-i)
                APSP[i][j] = 0; 
            else 
                APSP[i][j] = INT_MAX; 
}


pair<int, double> serailAdjAPSP(int **graph, int order, vector<vector<int> > &neighbours, int **APSP, bool **A, bool **B){

    initialiseAB(A,B,0,order,order);
    initialiseAPSP(APSP,0,order,order); 

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
    cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(diameter,average_distance); 
}


pair<int, double> serailDividedAdjAPSP(int **graph, int order, vector<vector<int> > &neighbours, int **APSP, bool **A, bool **B){
    int parSize = order/CHUNK_SIZE; 
    int diameter = 1; 
    int distance = order * (order-1); 
    double average_distance = 0.0; 
    int  num = 0; 
    //cout<<parSize; 
    //#pragma omp parallel for 
    for(int p = 0; p < parSize ;  p++){
        int k_outer;
        initialiseAB(A,B,p*CHUNK_SIZE,p*CHUNK_SIZE+CHUNK_SIZE,order);
        initialiseAPSP(APSP,p*CHUNK_SIZE,p*CHUNK_SIZE+CHUNK_SIZE,order); 
        for(int k = 0 ; k < order-1 ; k++){
            #pragma omp parallel for
            for(int i = 0; i < order ; i++){
                for(auto n: neighbours[i]){
                    for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                        B[i][j] = B[i][j] || A[n][j]; 
                    }

                }
            }
            num = 0; 
            #pragma omp parallel for reduction(+:num)
            for(int i = 0; i < order; i++)
                for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                    if(B[i][j]==1){
                        if(APSP[i][j] > k)
                            APSP[i][j] = k+1; 
                        num++;
                    }
                }
            if(num == order*CHUNK_SIZE) break;
            swap(A,B); 
            k_outer = k; 
            distance = distance+(order * CHUNK_SIZE - num); 

    }
    diameter = max(diameter,(k_outer+1)+1); 
    }

    average_distance = (double)distance/((order-1)* order); 
    cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(diameter,average_distance); 
}

int main(int argc, char *argv[]){
    string fileName = "n1000d5.random.edges"; 
    int order  = 1000;  
    int **graph = getAdjacencyMatrixArray( fileName,order); 
     
    int **APSP_serial = new int*[order]; 
    for(int i = 0; i < order; i++)
        APSP_serial[i] = new int[order];

    int **APSP_serial_div = new int*[order]; 
    for(int i = 0; i < order; i++)
        APSP_serial_div[i] = new int[order];  

    bool **A = new bool*[order]; 
    for(int i = 0; i < order; i++)
        A[i] = new bool[order]; 

    bool **B = new bool*[order]; 
    for(int i = 0; i < order; i++)
        B[i] = new bool[order]; 


    vector<vector<int> > neighbours = getAdjacencyListVector(fileName,order);
    double tt = omp_get_wtime();
	pair<int, double> diamAvgDist = serailAdjAPSP(graph, order, neighbours, APSP_serial,A,B);
	printf("sequential-adj-apsp = %fs\n", omp_get_wtime() - tt);

    tt = omp_get_wtime();
	pair<int, double> diamAvgDistDiv = serailDividedAdjAPSP(graph, order, neighbours, APSP_serial_div,A,B);
	printf("sequential-div-adj-apsp = %fs\n", omp_get_wtime() - tt);

    int error; 
    for(int i = 0; i < order; i++){
        for(int j = 0; j <  order; j++){
            error += fabs(APSP_serial[i][j]-APSP_serial_div[i][j]);
        }
    }

    printf("The calculated error is %d",error); 

    return 0; 
}





