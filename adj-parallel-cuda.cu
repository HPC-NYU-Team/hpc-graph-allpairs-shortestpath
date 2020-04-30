#include <iostream>
#include <climits> 
#include <omp.h>
#include <cstring>
#include <cmath>
#include "graphGenerator.cpp"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
using namespace std; 
int CHUNK_SIZE=5; 
int order = 10; 
int degree = 3; 
int nT = 4; 
double g_dist = 0.0; 
int g_diam = 0;
bool S_flag = false; 
bool O_flag = false; 
bool M_flag = false; 
string fileName = ""; 

//mpirun --mca btl_tcp_if_include eth0 --mca btl '^openib' --np 8 mpi 8192 5 1024

__global__ void logicalkernel(double *A, double *B, int *neighbours, int order ,int degree) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx <order)
        for(int i=0 ; i <  degree ; i++){
            int n = neighbours[idx*degree + i ]; 
            for(int j = 0; j < order; j++){
                B[i * order+ j] = B[i*order+j] || A[n*order+j]; 
            }
        }
    }
}
void initialiseABcuda(bool *A,bool *B, int order){
    for(int i = 0; i < order; i++)
    for(int j = 0; j < order; j++){
        if(j==order-1-i){
            A[i*order +j] = 1; 
            B[i*order +j] = 1;
        }
        else {
            A[i*order +j] = 0; 
            B[i*order +j] = 0; 
        }
    }
}
void initialiseAB(bool **A, bool **B, int r, int c, int offset, int arr_off){
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++){
            if(offset*c+j==r-1-i){
                A[i][c*arr_off +j] = 1; 
                B[i][c*arr_off +j] = 1;
            }
            else {
                A[i][c*arr_off +j] = 0; 
                B[i][c*arr_off +j] = 0; 
            }
        }
}


pair<int, double> serailAdjAPSP(int order, vector<vector<int> > &neighbours, bool **A, bool **B){
    initialiseAB(A,B,order,order,0,0);
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
                    num++;
                }
            }
        if(num == order*order) break; 
        swap(A,B); 
        diameter++; 
        distance = distance+(order * order - num); 
    }
    average_distance = (double)distance/((order-1)* order); 
    //cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(diameter,average_distance); 
}

pair<int, double> cudaAPSP(int order, bool *A_host, bool*B_host,bool* A_device, bool *B_device, int* nbr , int* nbr_device ){
    cudaMemcpy(nbr, nbr_device, order*degree*sizeof(int), cudaMemcpyHostToDevice);
    initialiseABcuda(A,B, order); 
    cudaMemcpy(A_host,A_device, order*order*sizeof(bool), cudaMemcpyHostToDevice); 
    cudaMemcpy(B_host,B_device, order*order*sizeof(bool), cudaMemcpyHostToDevice); 
    int diameter = 1; 
    int distance = order * (order-1); 
    double average_distance = 0.0; 

    int num; 
    for(int k = 0 ; k < order-1 ; k++){
        logicalkernel<<<order/1024 + 1, 1024>>>(A_device,B_device,nbr_device,order,degree); 
        // for(int i = 0; i < order ; i++){
        //     for(auto n: neighbours[i]){
        //         for(int j = 0 ; j < order ; j ++){
        //             B[i][j] = B[i][j] || A[n][j]; 
        //         }
        //     }
        // }
        num = 0; 
        #pragma omp parallel for reduction(+:num) num_threads(nT)
        for(int i = 0; i < order; i++)
            for(int j = 0; j < order ; j++){
                if(B[i][j]==1){
                    num++;
                }
            }
        if(num == order*order) break; 
        cudaMemcpy(A_device,B_device,order*order*sizeof(boolean),cudaMemcpyDeviceToDevice )
        //swap(A,B); 
        diameter++; 
        distance = distance+(order * order - num); 
    }
    average_distance = (double)distance/((order-1)* order); 
    //cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(diameter,average_distance); 
}

pair<int, double> serailDividedAdjAPSP(int order, vector<vector<int> > &neighbours,  bool **A, bool **B){
    int parSize = order/CHUNK_SIZE; 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int diameter = 1; 
    int distance = order * (order-1); 
    double average_distance = 0.0; 
    int  num = 0; 
    for(int p = 0; p < parSize ;  p++){
        int k_outer;
        initialiseAB(A,B,order,CHUNK_SIZE,p,p);
        for(int k = 0 ; k < order-1 ; k++){
            #pragma omp parallel for num_threads(nT)
            for(int i = 0; i < order ; i++){
                for(auto n: neighbours[i]){
                    for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                        B[i][j] = B[i][j] || A[n][j]; 
                    }
                }
            }
            num = 0; 
            #pragma omp parallel for reduction(+:num) num_threads(nT)
            for(int i = 0; i < order; i++)
                for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                    if(B[i][j]==1){
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
    //cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(diameter,average_distance); 
}


void get_file_name(string order, string degree); 

void getOptions(int argc, char** argv); 

int main(int argc, char *argv[]){
    double time_Serial = 0.0; 
    double time_OpenMP = 0.0; 
    double time_cuda = 0.0; 
    
    /*----Setting parameters-------*/ 
    getOptions(argc, argv); 
 
    /*----Getting neighbors-------*/ 
    vector<vector<int> > neighbours = getAdjacencyListVector(fileName,order);
    //int [] [] neighbours; 
    double tt; 
    {
        if(S_flag || O_flag){
            bool **A = new bool*[order]; 
            for(int i = 0; i < order; i++)
                A[i] = new bool[order]; 

            bool **B = new bool*[order]; 
            for(int i = 0; i < order; i++)
                B[i] = new bool[order]; 

            if(S_flag){
                tt = omp_get_wtime();
	            pair<int, double> diamAvgDist = serailAdjAPSP(order, neighbours,A,B);
                time_Serial = omp_get_wtime() - tt; 
	            //printf("sequential-adj-apsp = %fs\n", time_Serial);
            }

            if(O_flag){
                tt = omp_get_wtime();
	            pair<int, double> diamAvgDistDiv = serailDividedAdjAPSP(order, neighbours,A,B);
                time_OpenMP = omp_get_wtime() - tt; 
	            //printf("sequential-div-adj-apsp = %fs\n", time_OpenMP); 
            }

            for(int i = 0; i <order;i++){
                free(A[i]); 
                free(B[i]); 
            }
            free(A); 
            free(B); 
        }
    } 

    
    if(M_flag)
    { 
       
        int *nbr = new int[order*degree]; 
        for(int i=0; i<order; i++)
            for(int j=0; j<degree; j++){
                nbr[i*degree+j] = neighbours[i][j]; 
            }
        //change this to an array
        bool *A_host = new bool[order*order]; 
        bool *B_host = new bool[oder*order];
        int *nbr_device; 
        bool *A_device, *B_device; 
        cudaMalloc(&nbr_device, order*degree*sizeof(int));  
        cudaMalloc(&A_device, order*order*sizeof(bool)); 
        cudaMalloc(&B_device, order*order*sizeof(bool)); 
        
        tt = omp_get_wtime();
        cudaAPSP(order,A_host,B_host,A_device,B_device,nbr,nbr_device);
        time_cuda =  omp_get_wtime() - tt; 
            //printf("Parallel-div-adj-apsp = %fs\n", time_MPI);
        free(A_host);
        free(B_host); 
    }

    printf(" Order  Degree  Chunk_size  Threads Time_Serial Time_OpenMP   Time_MPI  Dist Diam\n");
    printf("%5d %5d %10d %10d %12f %12f %12f %12f %12d\n",order,degree,CHUNK_SIZE,nT,time_Serial,time_OpenMP,time_cuda,g_dist,g_diam); 
    

    return 0; 
}





void get_file_name( string order, string degree){
    //fileName = "../Data/n";
    fileName = "n"; 
    fileName += order; 
    fileName += "d"; 
    fileName += degree;
    fileName += ".random.edges"; 
}

void getOptions(int argc, char** argv){
    int opt;  
    int arg_number=1; 
    string flags; 
    while((opt = getopt(argc, argv, "o:d:c:t:f:")) != -1)  
    {  
        switch(opt)  
        {  
            case 'o':{
                arg_number++; 
                sscanf(optarg, "%d",&order); 
                break; 

            } 
            case 'd':{
                arg_number++; 
                sscanf(optarg, "%d",&degree); 
                break;
            }
            case 'c':{
                arg_number++; 
                sscanf(optarg, "%d",&CHUNK_SIZE); 
                break;
            }
            case 't':{
                arg_number++;
                sscanf(optarg, "%d",&nT); 
                break; 
            }
            case 'f':{
                string options; 
                arg_number++; 
                flags = optarg;  
                for(int i = 0; i < flags.size(); i++){
                    switch(flags[i]){
                        case 'S':{
                            S_flag = true; 
                            break; 
                        }
                        case 'O':{
                            O_flag = true; 
                            break; 
                        }
                        case 'M':{
                            M_flag = true; 
                            break; 
                        }
                    }
                }
                break;

            }
            case ':':  { cout<<"option needs a value\n";   
                break;  }
            case '?':{ cout<<"unknown option:"<<(char)opt; 
                break;  }
        }  
    }  

    get_file_name(to_string(order),to_string(degree)); 
}
