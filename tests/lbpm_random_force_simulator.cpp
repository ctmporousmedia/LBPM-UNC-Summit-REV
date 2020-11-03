#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include "fnGenerateSphereTest.h"
//#include "fn_lbpm_segmented_pp.h"
#include "fn_lbpm_random_pp.h"


#include "models/ColorModel.h"

using namespace std;

//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
//*************************************************************************

int main(int argc, char **argv)
{
    // Initialize MPI
    int provided_thread_support = -1;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided_thread_support);
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD,&comm);
    int rank = comm_rank(comm);
    int nprocs = comm_size(comm);
    if ( rank==0 && provided_thread_support<MPI_THREAD_MULTIPLE )
        std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
    { // Limit scope so variables that contain communicators will free before MPI_Finialize
        
        if (rank == 0){
            printf("\n========================================\n");
            printf("LATTICE BOLTZMANN SIMULATOR\n");
            printf("========================================\n");
        }
        PROFILE_ENABLE(1);
        //PROFILE_ENABLE_TRACE();
        //PROFILE_ENABLE_MEMORY();
        PROFILE_SYNCHRONIZE();
        PROFILE_START("Main");
        Utilities::setErrorHandlers();
        
        
        auto filename = argv[1];
        
        
        
        auto db = std::make_shared<Database>( filename );
        auto color_db = db->getDatabase( "Color" );
        int FluidInitialization = color_db->getScalar<int>( "FluidInitialization" );
        
        // If restart: Don't do pre-processing
        bool Restart;
        Restart = color_db->getScalar<bool>( "Restart" );
        
        
        if (Restart == false) {
            if (FluidInitialization == 1) lbpm_random_pp(argc, argv, rank,nprocs,comm);
            else if (FluidInitialization == 2) lbpm_partial_saturation_pp(argc, argv, rank,nprocs,comm);
            else if (FluidInitialization == 3) lbpm_morphdrain_pp(argc,argv,rank,nprocs,comm);
            else if (FluidInitialization == 4) lbpm_morphopen_pp(argc,argv,rank,nprocs,comm);
        }
        
        
        // End
        ScaLBL_ColorModel ColorModel(rank,nprocs,comm);
        ColorModel.ReadParams(filename);
        ColorModel.SetDomain();
        ColorModel.ReadInput();
        ColorModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
        ColorModel.Initialize();   // initializing the model will set initial conditions for variables
        ColorModel.Run();
        //ColorModel.WriteDebug();
        
        PROFILE_STOP("Main");
        PROFILE_SAVE("lbpm_color_simulator",1);
        // ****************************************************
        MPI_Barrier(comm);
    } // Limit scope so variables that contain communicators will free before MPI_Finialize
    MPI_Comm_free(&comm);
    MPI_Finalize();
}


