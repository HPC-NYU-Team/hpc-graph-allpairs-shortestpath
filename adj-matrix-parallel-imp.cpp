#include <iostream>
#include <climits> 
#include <omp.h>
#include <mpi.h> 
#include <cstring>
#include <cmath>
#include "graphGenerator.cpp"
using namespace std; 
int CHUNK_SIZE=5; 

template <class T>
void print(int order, T **A); 

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


void initialiseAPSP(int *APSP,int r,int col_max,int c,int offset, int arr_off){
   for(int i = 0; i < r; i++){
       //cout<<endl;
        for(int j = 0; j < c; j++)
        {     
            if(offset*c +j==r-1-i){
                //cout<<i<<j<<" "; 
                APSP[i*col_max + c*arr_off +j] = 0; 
            }
            else {
                //cout<<"x"<<"  "; 
                APSP[i*col_max + c*arr_off +j] = INT_MAX;  
            }
        } 
   }
}

pair<int, double> serailAdjAPSP(int order, vector<vector<int> > &neighbours, int *APSP, bool **A, bool **B){

    initialiseAB(A,B,order,order,0,0);
    initialiseAPSP(APSP,order,order,order,0,0); 

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
                    if(APSP[i*order + j] > k)
                        APSP[i*order + j] = k+1; 
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
pair<int, double> serailDividedAdjAPSP(int order, vector<vector<int> > &neighbours, int *APSP, bool **A, bool **B){
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
        initialiseAPSP(APSP,order,order,CHUNK_SIZE,p,p);
        for(int k = 0 ; k < order-1 ; k++){
            #pragma omp parallel for num_threads(8)
            for(int i = 0; i < order ; i++){
                for(auto n: neighbours[i]){
                    for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                        B[i][j] = B[i][j] || A[n][j]; 
                    }
                }
            }
            num = 0; 
            #pragma omp parallel for reduction(+:num) num_threads(8)
            for(int i = 0; i < order; i++)
                for(int j = p*CHUNK_SIZE ; j < p*CHUNK_SIZE+CHUNK_SIZE ; j ++){
                    if(B[i][j]==1){
                        if(APSP[i*order + j] > k)
                            APSP[i*order + j] = k+1; 
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
    initialiseAPSP(APSP,order , chunk , chunk, rank,0); 
    for(int k = 0 ; k < order-1 ; k++){
        #pragma omp parallel for num_threads(8)
        for(int i = 0; i < order ; i++){
            for(auto n: neighbours[i]){
                for(int j = 0 ; j < chunk ; j ++){
                    B[i][j] = B[i][j] || A[n][j]; 
                }
            }
        }
        num = 0; 
        #pragma omp parallel for reduction(+:num) num_threads(8)
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
    if(rank == 0)
        cout<<diameter<<" "<<average_distance<<endl; 
    return make_pair(global_diam,average_distance); 
}
void get_file_name(string &file , string order, string degree); 

int main(int argc, char *argv[]){
    /*----Initialising openMPI-------*/ 
    MPI_Init(&argc, &argv); 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD , &rank);

    /*----Setting parameters-------*/ 
    string fileName; 
    int order = 10; 
    int degree = 3; 
    if(argc < 4 ){
        if(rank==0)
            cout<<"Too few arguments, taking default values"<<endl; 
        fileName = "n10d3.random.edges"; 
    }
    else{
        get_file_name(fileName , argv[1], argv[2]); 
        order = atoi(argv[1]);
        degree = atoi(argv[2]);  
        CHUNK_SIZE = atoi(argv[3]); 
    }
 
    /*----Getting neighbors-------*/ 
    vector<vector<int> > neighbours = getAdjacencyListVector(fileName,order);
    int *APSP_serial; 
    double tt; 
    if(rank == 0){
        APSP_serial = new int[order*order]; 

        int *APSP_serial_div = new int[order*order]; 
        bool **A = new bool*[order]; 
        for(int i = 0; i < order; i++)
            A[i] = new bool[order]; 

        bool **B = new bool*[order]; 
        for(int i = 0; i < order; i++)
            B[i] = new bool[order]; 
    
        tt = omp_get_wtime();
	    pair<int, double> diamAvgDist = serailAdjAPSP(order, neighbours, APSP_serial,A,B);
	    printf("sequential-adj-apsp = %fs\n", omp_get_wtime() - tt);

        tt = omp_get_wtime();
	    pair<int, double> diamAvgDistDiv = serailDividedAdjAPSP(order, neighbours, APSP_serial_div,A,B);
	    printf("sequential-div-adj-apsp = %fs\n", omp_get_wtime() - tt);

        int error; 
        for(int i = 0; i < order; i++){
            for(int j = 0; j <  order; j++){
                error += fabs(APSP_serial[i*order + j]-APSP_serial_div[i*order + j]);
            }
        }

        printf("The calculated error is %d\n",error); 
        for(int i = 0; i <order;i++){
            free(A[i]); 
            free(B[i]); 
        }
        free(A); 
        free(B); 
        free(APSP_serial_div);

    } 
    MPI_Barrier(MPI_COMM_WORLD); 
    
    //MPI part 
    { 
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
        tt = omp_get_wtime();
        parallelDividedAdjAPSP( order,chunk, neighbours, APSP_parallel_div_sub, A_sub, B_sub); 
        MPI_Barrier(MPI_COMM_WORLD); 
        //printf("Parallel-div-adj-apsp = %fs\n", omp_get_wtime() - tt);
        if(rank!=0)
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
            //end
            printf("Parallel-div-adj-apsp = %fs\n", omp_get_wtime() - tt);
            int error=0; 
            for(int i = 0; i < order; i++){
                //cout<<endl;
                for(int j = 0; j <  order; j++){
                    //cout<<APSP_serial[i][j]<<" "; 
                    error += fabs(APSP_serial[i*order + j]-APSP_parallel_div[i][j]);
                }
            }
            printf("The calculated error is %d\n",error);
            
            free(APSP_serial); 
        }

        
    }
    MPI_Finalize();
    
    return 0; 
}




template <class T>
void print(int order, T **A){
    for(int i =  0 ; i<order; i++){    
        for(int j = 0; j<order; j++)
            cout<<A[i][j]<<" "; 
        cout<<endl; 
    }
} 

void get_file_name(string &file , string order, string degree){
    file = "n";
    file += order; 
    file += "d"; 
    file += degree;
    file += ".random.edges"; 

}