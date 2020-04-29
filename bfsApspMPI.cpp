#include <iostream>
#include <omp.h>
#include <algorithm> 
#include "graphGenerator.cpp"
#include <numeric>
#include <mpi.h>
#include <math.h>
#include <cstring>
#include <getopt.h>

using namespace std;


int RANK;
int PROCS;
int THREADS;


void get_file_name(string &file , int order, int degree){
    file = "n";
    file += to_string(order); 
    file += "d"; 
    file += to_string(degree);
    file += ".random.edges"; 
}


void getAllChildren(vector<vector<int>> &g, vector<int> &frontier, vector<int> &children, int* distances, int k){
	//printf("%d\n", frontier.size());
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

void getAllChildrenParallel(vector<vector<int>> &g, vector<int> &frontier, vector<int> &children, int* distances, int k){
	int count = 0;
	vector<int> local_children = vector<int>{};
	#pragma omp parallel private(local_children) num_threads(THREADS)
	{
		#pragma omp for nowait
		for(int i=0; i<frontier.size(); i++){
				int node = frontier[i];
				for(int j=0; j<g[node].size(); j++){
						int n = g[node][j];
						if(distances[n] == -1){
								distances[n] = k;
								local_children.push_back(n);
						}
				}
		}
		#pragma omp critical
		move(local_children.begin(), local_children.end(), back_inserter(children));
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

int* bfsParallel(vector<vector<int>> g, int source){
        vector<int> frontier = {source};
        vector<int> children = vector<int>{};

        int *distances = new int[g.size()];
        for(int i=0; i<g.size(); i++){
                distances[i] = -1;
        }

        int k = 1;
        while(frontier.size() > 0){
                getAllChildrenParallel(g, frontier, children, distances, k++);
                frontier = children;
                children = {};
        }
        return distances;
}

pair<float, float> serialBfsApsp(vector<vector<int>> g){
	int chunkSize = ceil(1.0*g.size()/PROCS);
	int startSrc = chunkSize * RANK;
        int endSrc   = chunkSize * (RANK+1);
	if(startSrc > g.size())
		startSrc = g.size();
        if(endSrc > g.size())
                endSrc = g.size();
	int currSize = endSrc-startSrc;
	//printf("%d %d %d", startSrc, endSrc, currSize);
	
	float* diameter = new float[currSize];
	float* distance = new float[currSize];
	
	for(int i=0; i<currSize; i++){
		int* distances = bfs(g, startSrc+i);	
		diameter[i] = *max_element(distances, distances + g.size());
		distance[i] = accumulate(distances, distances + g.size(), 0.0) / g.size();
	}
	double maxDia = *max_element(diameter, diameter + currSize);
	double avgDist = accumulate(distance, distance + currSize, 0.0) / currSize;
	return make_pair(maxDia, avgDist);
}

pair<float, float> parallelBfsApsp(vector<vector<int>> g){
	if(PROCS == 1){
		int source = 0;
		float* diameter = new float[g.size()];
		float* distance = new float[g.size()];
		
		#pragma omp parallel for num_threads(THREADS)
		for(int source=0; source<g.size(); source++){
			int* distances = bfsParallel(g, source);	
			diameter[source] = *max_element(distances, distances + g.size());
			distance[source] = accumulate(distances, distances + g.size(), 0.0) / g.size();
		}
		float global_max_dia = *max_element(diameter, diameter + g.size());
		float global_avg_dist = (accumulate(distance, distance + g.size(), 0.0) / g.size());
		return make_pair(global_max_dia, global_avg_dist);
	}
	int chunkSize = ceil(1.0*g.size()/PROCS);
	int startSrc = chunkSize * RANK;
        int endSrc   = chunkSize * (RANK+1);
	if(startSrc > g.size())
                startSrc = g.size();
        if(endSrc > g.size())
                endSrc = g.size();
	int currSize = endSrc-startSrc;
	float* diameter = new float[currSize];
	float* distance = new float[currSize];
	
	#pragma omp parallel for num_threads(THREADS)
	for(int i=0; i<currSize; i++){
		int* distances = bfsParallel(g, startSrc+i);
		diameter[i] = *max_element(distances, distances + g.size());
		distance[i] = accumulate(distances, distances + g.size(), 0.0) / g.size();
	}
	float local_max_dia = *max_element(diameter, diameter + currSize);
	float local_avg_dist = (accumulate(distance, distance + currSize, 0.0) / currSize);
	//cout << "Proc:" << RANK << "  " << "Local Max diameter: " << local_max_dia << endl;
	//cout << "Proc:" << RANK << "  " << "Local Avg distance: " << local_avg_dist << endl;
	float global_max_dia, global_avg_dist;
	MPI_Reduce(&local_max_dia, &global_max_dia, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local_avg_dist, &global_avg_dist, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(RANK == 0){
		// The ceil might give us an extra empty 0, take that into account, dont avg it. Could be  PROCS/PROCS-1. TODO
		global_avg_dist /= (1.0*PROCS);
	}
	return make_pair(global_max_dia, global_avg_dist);
}



int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &PROCS);
	MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
	THREADS = 1;
	
	string fileName;
	int order = -1;
	int degree = -1;
	int serial = false;
	int parallelBfs = false;
	int threads = 1;
	int procs = 1;

	while (1) {
		static struct option long_options[] =
			{
				{"serial",   	no_argument,       &serial, 1},
				{"parallelbfs",no_argument,       &parallelBfs, 1},
				{"order",  	required_argument, 0, 'o'},
				{"degree",  	required_argument, 0, 'd'},
				{"threads",  	required_argument, 0, 't'},
				{"procs",  	required_argument, 0, 'p'},
				{0, 0, 0, 0}
			};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		int c = getopt_long (argc, argv, "o:d:t:p:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			break;

			case 'o':
			order = atoi(optarg);
			break;

			case 'd':
			degree = atoi(optarg);
			break;

			case 't':
			threads = atoi(optarg);
			break;

			case 'p':
			procs = atoi(optarg);
			break;

			case '?':
			/* getopt_long already printed an error message. */
			break;

			default:
			abort ();
		}
	}

	THREADS = threads;

	
	if(RANK == 0){
		if(order == -1 || degree == -1){
			printf("Order and degree is not set\n");
			MPI_Finalize();
			return 0;
		}
	}
	
    
	vector<vector<int>> g(order, vector<int>(degree));
	int gsize = order*degree;

	if(RANK == 0){
		get_file_name(fileName, order, degree); 
	    g = getAdjacencyListVector(fileName,order);
	}
		
	// If Serial Run	
	if(serial && RANK == 0){
		double tt = omp_get_wtime();
		pair<float, float> res = serialBfsApsp(g);
		double timeTaken = omp_get_wtime() - tt;
		printf("Sequential: Order = %d\t Degree = %d\tMaxDia = %f\tAvgDist = %f\tTime = %f\n", order, degree, res.first, res.second, timeTaken);
	}

	
	vector<int> newg;
	if(RANK == 0){
		for(auto && v : g){
			newg.insert(newg.end(), v.begin(), v.end());
		}
	}
	newg.resize(gsize);

	MPI_Bcast(&newg[0], gsize, MPI_INT, 0, MPI_COMM_WORLD);
	if(RANK != 0){
		for(int i=0; i<newg.size(); i++){
			g[i/degree][i%degree] = newg[i];
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 
	double tt2 = MPI_Wtime();
	pair<float, float> res = parallelBfsApsp(g);
	if(RANK == 0){
		double timeTaken = MPI_Wtime() - tt2;
		printf("Parallel: Procs= %d\tThreads = %d\tOrder = %d\t Degree = %d\tMaxDia = %f\tAvgDist = %f\tTime = %f\n", procs, threads, order, degree, res.first, res.second, timeTaken);
	}

	
	MPI_Finalize();
	return 0;

}


