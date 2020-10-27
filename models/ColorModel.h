/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
Implementation of color lattice boltzmann model
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"
#include "common/MPI_Helpers.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"

class ScaLBL_ColorModel{
public:
    ScaLBL_ColorModel(int RANK, int NP, MPI_Comm COMM);
    ~ScaLBL_ColorModel();
    
    // functions in they should be run
    size_t analysis_interval;
    
    bool legacySDs;
    bool Restart,pBC;
    int timestep,timestepMax;
    int BoundaryCondition;
    double tauA,tauB,rhoA,rhoB,alpha,beta;
    double Fx,Fy,Fz,flux;
    double din,dout,inletA,inletB,outletA,outletB;
    
    size_t Nx,Ny,Nz,N,Np,Ni;
    
    int rank,nprocx,nprocy,nprocz,nprocs;
    double Lx,Ly,Lz;
    int method;

    std::shared_ptr<Domain> Dm;   // this domain is for analysis
    std::shared_ptr<Domain> Mask; // this domain is for lbm
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular;
    std::shared_ptr<TwoPhase> Averages;
    
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> color_db;
    std::shared_ptr<Database> analysis_db;

    IntArray Map,InactiveMap;
    char *id;
    int *NeighborList;
    int *InterpolationList;
    double * LIBBqA;
    double * LIBBqBC;
    double * LIBBqD;
    
    bool Berea = false;
    bool updateQ = false;
    
    bool * Boolfield;
    char * ID;
    int *dvcMap;
    int *dvcInactiveMap;
    double *fq, *Aq, *Bq, *Aq2, *Bq2, *Qdist;
    double * D_coefs;
    double * DD_coefs;
    double * Phi_fictitious;
    double * fq2;
    double * savedfq;
    double * savedAq, * savedBq;
    double * vfmask;
    double * dvcSignDist;
    double * SignDist;
    double * DenA, *DenB;
    
    double * laggedDenA, *laggedDenB;
    double * analysisDenA, *analysisDenB;
    
    double * MassA, *MassB;
    double * Phi;
    double * Phi2;
    double * FVA; // Finite Volume Area
    double * ColorGrad;
    double * Velocity;
    double * Pressure;
    int * Xi;
    double * force;
    
    double * GradPhiX;
    double * GradPhiY;
    double * GradPhiZ;
    
    double * gradphix;
    double * gradphiy;
    double * gradphiz;
    
public:
    double porosity;
    double affinity;
    double starttime,stoptime,cputime;
    MPI_Comm comm;
    
    int dist_mem_size;
    int neighborSize;
    
    size_t Npad;

    /* Fill in the qDistances arrays */
    size_t Nq;
    size_t NqN;
    size_t NqNpad;
    
    
    /* local CPU registers */
    double* qDistances;
    int* neighborList;
    int* interpolationList;
    double* libbqA;
    double* libbqBC;
    double* libbqD;
    double *PhaseLabel;
    bool* boolfield;
    double* empty;
    int *TmpMap;
    
    double af; // affinity 
    size_t NqNp;
    size_t libbSize;
    size_t signSize;
    size_t qsize;
    // filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);
    void AssignComponentLabels(double *phase);
    double MorphInit(const double beta, const double morph_delta);
    
    void ReadRestartFiles();
    void UpdateMapsAndLists();
    void InitializePhaseAndDensityFields();
    void ReadSignedDistanceField();
    
    void InitializeSolidVolumeFraction();
    void InitializeInterpolationArrays();
    
    void SaveMomentumDistribution_1();
    void SaveMomentumDistribution_2();
    
    // NonAA with LIBB
    void LIBBNonAAPhaseField_1();
    void LIBBNonAAPhaseField_2();
    
    void LIBBNonAAStreamAndCollide_1(int timestep);
    void LIBBNonAAStreamAndCollide_2(int timestep);
    
    // AA with LIBB
    void LIBBwithAAPhaseField_Odd();
    void LIBBwithAAPhaseField_Even();
    
    void FillInactiveDensities();
    
    void LIBBwithAAStreamAndCollideAndStream_Odd();
    void LIBBwithAACollide_Even();
    
    void Analyze();
    
        
public:
    void ReadParams(string filename);
      void ReadParams(std::shared_ptr<Database> db0);
      void SetDomain();
      void ReadInput();
      bool ReadInputAndValidate();
      void Create();
      bool GridIndependence();
      void Initialize();
      void Run();
      void WriteDebug();
    void RestartDistances_pp();
};



