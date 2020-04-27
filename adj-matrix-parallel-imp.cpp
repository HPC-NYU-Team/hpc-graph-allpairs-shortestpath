#include <iostream>
#include <climits> 
#include <omp.h>
#include <mpi.h> 
#include <cmath>
#include "graphGenerator.cpp"
#define CHUNK_SIZE 5
using namespace std; 
template <class T>
void print(int order, T **A){
    for(int i =  0 ; i<order; i++){    
        for(int j = 0; j<order; j++)
            cout<<A[i][j]<<" "; 
        cout<<endl; 
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
// void initialiseAB(bool **A, bool **B, int col_beg, int col_end,int order){
//     for(int i = 0; i < order; i++)
//         for(int j = col_beg; j < col_end; j++)
//         {    
//             if(j==order-1-i){
//                 A[i][j] = 1; 
//                 B[i][j] = 1;
//             }
//             else {
//                 A[i][j] = 0; 
//                 B[i][j] = 0; 
//             }
//         }

// }

// void initialiseABparallel(bool **A, bool **B,int chunk, int offset,int order){
//     for(int i = 0; i < order; i++)
//         for(int j = 0; j < chunk; j++)
//         {    
//             if(offset*chunk+j==order-1-i){
//                 A[i][j] = 1; 
//                 B[i][j] = 1;
//             }
//             else {
//                 A[i][j] = 0; 
//                 B[i][j] = 0; 
//             }
//         }
// }

void initialiseAPSP(int **APSP,int col_beg, int col_end,int order){
    for(int i = 0; i < order; i++)
        for(int j = col_beg; j < col_end ; j++)
            if(j==order-1-i)
                APSP[i][j] = 0; 
            else 
                APSP[i][j] = INT_MAX; 
}


void initialiseAPSPparallel(int *APSP,int order,int chunk,int offset){
   for(int i = 0; i < order; i++)
        for(int j = 0; j < chunk; j++)
        {    
            if(offset*chunk+j==order-1-i){
                APSP[i*chunk + j] = 0; 
            }
            else {
                APSP[i*chunk + j] = INT_MAX;  
            }
        } 
}

pair<int, double> serailAdjAPSP(int order, vector<vector<int> > &neighbours, int **APSP, bool **A, bool **B){

    initialiseAB(A,B,order,order,0,0);
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
pair<int, double> serailDividedAdjAPSP(int order, vector<vector<int> > &neighbours, int **APSP, bool **A, bool **B){
    int parSize = order/CHUNK_SIZE; 
    int diameter = 1; 
    int distance = order * (order-1); 
    double average_distance = 0.0; 
    int  num = 0; 
    for(int p = 0; p < parSize ;  p++){
        int k_outer;
        initialiseAB(A,B,order,CHUNK_SIZE,p,p);
        //initialiseAB(A,B,p*CHUNK_SIZE,p*CHUNK_SIZE+CHUNK_SIZE,order);
        initialiseAPSP(APSP,p*CHUNK_SIZE,p*CHUNK_SIZE+CHUNK_SIZE,order); 
        for(int k = 0 ; k < order-1 ; k++){
            #pragma omp parallel for num_threads(4)
            for(int i = 0; i < order ; i++){
                for(auto n: neighbours[i]){
                    for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                        B[i][j] = B[i][j] || A[n][j]; 
                    }
                }
            }
            num = 0; 
            #pragma omp parallel for reduction(+:num) num_threads(4)
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


pair<int, double> parallelDividedAdjAPSP( int order,int chunk, vector<vector<int> > &neighbours, int *APSP, bool **A, bool **B){

    int diameter = 1; 
    double distance = 0; 
    double average_distance = 0.0; 
    int  num = 0; 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    initialiseAB(A,B,order,chunk,rank,0);
    //initialiseABparallel(A,B,chunk, rank, order); 
    initialiseAPSPparallel(APSP, order , chunk , rank); 
    for(int k = 0 ; k < order-1 ; k++){
        #pragma omp parallel for num_threads(4)
        for(int i = 0; i < order ; i++){
            for(auto n: neighbours[i]){
                for(int j = 0 ; j < chunk ; j ++){
                    B[i][j] = B[i][j] || A[n][j]; 
                }
            }
        }
        num = 0; 
        #pragma omp parallel for reduction(+:num) num_threads(4)
        for(int i = 0; i < order; i++)
            for(int j = 0 ; j < chunk ; j ++){
                if(B[i][j]==1){
                    if(APSP[i*chunk + j] > k)
                        APSP[i*chunk + j] = k+1; 
                    num++;
                }
            }
        if(num == order*chunk) break;
        swap(A,B); 
        diameter++; 
        distance = distance+(order * chunk - num); 
    }
    int global_diam;
    MPI_Allreduce(&diameter, &global_diam, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&distance, &average_distance, 1, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD);
    average_distance /= ((order-1)* order);
    average_distance +=1;  
    cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(global_diam,average_distance); 
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv); 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD , &rank);
    string fileName = "n10d3.random.edges"; 
    int order  = 10;  
     
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
    if(rank == 0){
    double tt = omp_get_wtime();
	pair<int, double> diamAvgDist = serailAdjAPSP(order, neighbours, APSP_serial,A,B);
	printf("sequential-adj-apsp = %fs\n", omp_get_wtime() - tt);

    tt = omp_get_wtime();
	pair<int, double> diamAvgDistDiv = serailDividedAdjAPSP(order, neighbours, APSP_serial_div,A,B);
	printf("sequential-div-adj-apsp = %fs\n", omp_get_wtime() - tt);

    int error; 
    for(int i = 0; i < order; i++){
        for(int j = 0; j <  order; j++){
            error += fabs(APSP_serial[i][j]-APSP_serial_div[i][j]);
        }
    }

    printf("The calculated error is %d\n",error); 
    } 
    MPI_Barrier(MPI_COMM_WORLD); 
    //MPI part 
    {
        cout<<"\n\n Final computation\n\n"; 
        int npes; 
        MPI_Comm_size(MPI_COMM_WORLD , &npes);
        int chunk = order/npes; 

        bool **A_sub = new bool*[order]; 
        for(int i = 0; i < order; i++)
            A_sub[i] = new bool[chunk]; 

        bool **B_sub = new bool*[order]; 
        for(int i = 0; i < order; i++)
            B_sub[i] = new bool[chunk]; 

        int *APSP_parallel_div_sub = new int[order*chunk]; 

        parallelDividedAdjAPSP( order,chunk, neighbours, APSP_parallel_div_sub, A_sub, B_sub);  
        MPI_Barrier(MPI_COMM_WORLD); 

        MPI_Send(APSP_parallel_div_sub, chunk*order , MPI_INT,0,0,MPI_COMM_WORLD); 
        if(rank == 0){
            
            int **APSP_parallel_div = new int*[order]; 
            for(int i = 0; i < order; i++)
                APSP_parallel_div[i] = new int[order];

            for(int i = 0 ; i < order ; i++){
                    //cout<<endl; 
                    for(int j = 0 ; j < chunk ; j++ ){
                       //cout<<i<<j<<" ";
                       APSP_parallel_div[i][j] = APSP_parallel_div_sub[i*chunk + j]; 
                        
                    }
                }
            
            int *temp = new int[order*chunk]; 
            for(int n = 1 ; n < npes ; n++){
                MPI_Recv(temp, chunk*order, MPI_INT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int i = 0 ; i < order ; i++){
                    for(int j = 0 ; j < chunk ; j++ ){
                       APSP_parallel_div[i][n*chunk + j] = temp[i*chunk + j]; 
                        
                    }
                }
            }

            int error=0; 
            for(int i = 0; i < order; i++){
                //cout<<endl;
                for(int j = 0; j <  order; j++){
                    //cout<<APSP_serial[i][j]<<" "; 
                    error += fabs(APSP_serial[i][j]-APSP_parallel_div[i][j]);
                }
            }
            printf("The calculated error is %d\n",error);
        }

        
    }
    MPI_Finalize();
    return 0; 
}





