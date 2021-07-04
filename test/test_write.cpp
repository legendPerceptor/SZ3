//
// Created by apple on 2021/7/2.
//

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <utils/FileUtil.h>

int main(int argc, char** argv) {
    std::cout << "argc: " << argc << std::endl;
    if(argc != 3) {
        std::cout << "provide input and output files" << std::endl;
        exit(-1);
    }
    srand(time(0));
    // MPI Initialization
    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    size_t num=0;
    std::unique_ptr<float[]> data;
    data = SZ::readfile<float>(argv[1], num);
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime, endTime;
    double cost;
    if(world_rank==0){
        std::cout<< "Write time record" << std::endl;
        startTime = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    char filename[100];
    sprintf(filename, "%s/file_%d.out", argv[2], world_rank);
    SZ::writefile(filename, data.get(), num);
    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank==0){
        endTime = MPI_Wtime();
        cost = endTime - startTime;
        std::cout << "The write time is " << cost << " seconds" << std::endl;
    }
    return 0;
}