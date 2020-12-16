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
#include "models/ColorModel.h"
#include "analysis/distance.h"
#include <vector>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

int CM_DOUBLE = 8;
int CM_INT = 1;




bool ScaLBL_ColorModel::GridIndependence() { return true; }
bool ScaLBL_ColorModel::ReadInputAndValidate(){ return 1; } // Not used
double ScaLBL_ColorModel::MorphInit(const double beta, const double morph_delta) { return 0; }
void ScaLBL_ColorModel::WriteDebug(){}
inline void EffectiveAreas(char* id, int N, double * VFmask, double * VolumeMask, int Nx, int Ny, int Nz) {}




inline void NormalizeQdistances(char * id, int N, double * qDistances, double * libbqA, double * libbqBC, double * libbqD) {
    for (int q = 0; q < 18; q++) {
        for (int n = 0; n < N; n++) {
            if (id[n] == 0) {
                libbqA[n+q*N] = 0;
                libbqBC[n+q*N] = 0;
                libbqD[n+q*N] = 0;
                qDistances[n + q*N] = 0.0;
            }
        }
    }
    /* Normalizing the diagonal directions so that the A,BC,D values can be computed according to < or >= 0.5 rule */
    for (int q = 0; q < 6; q++) { // on-axis directions
        for (int n = 0; n < N; n++) {
            if (qDistances[n + q*N] == -1.0) {
                qDistances[n + q*N] = 1.0; // elements that weren't modified, set to 1
            }
        }
    }
    for (int q = 6; q < 18; q++) { // diagonal directions
        for (int n = 0; n < N; n++) {
            if (qDistances[n + q*N] == -1.0) {
                qDistances[n + q*N] = sqrt(2); // elements that weren't modified, set to sqrt(2)
            }
        }
    }
    
    /* Should I normalize the diagonal directions (note: directions 0,5 are already between 0 and 1 */
    for (int q = 6; q < 18; q++) {
        for (size_t n = 0; n < N; n++) {
            qDistances[n + q*N] /= sqrt(2); // Normalize diagonal directions for the A,BC,D arrays
        }
    }
}

inline void ComputeFiniteDifferenceStencilCoefficients(double *qDistances, double *d_coefs, double * dd_coefs, int N) {}


inline void save_3Dchar_parallel(int rank, char* field, int Nx, int Ny, int Nz, int depth, std::string s)
{
    assert(depth >= 0);
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    FILE *f2 = fopen(LocalRankFilename, "w");
    
    char var = 0;
    size_t c = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z;
                
                var = field[n];
                fprintf(f2, "%d %d %d %d\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}



inline double calculate_porosity_double(double * id, int Nx, int Ny, int Nz, Domain &Dm, int rank, double af) {
    int nprocx,nprocy,nprocz;
    
    nprocx = Dm.nprocx();
    nprocy = Dm.nprocy();
    nprocz = Dm.nprocz();
    size_t count = 0;
    size_t totalGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                size_t n = k*Nx*Ny+j*Nx+i;
                if (id[n] == af)  count++;
            }
        }
    }
    
    // total Global is the number of nodes in the pore-space
    MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,Dm.Comm);
    double volume=double(nprocx*nprocy*nprocz)*double(Nx-2)*double(Ny-2)*double(Nz-2);
    double porosity=double(totalGlobal)/volume;
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    double dcount = 0;
    double dcountGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                for (int p=0;p<8;p++){
                    int n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                    if (id[n] == af)  {  dcount += 0.125; }
                }
            }
        }
    }
    MPI_Allreduce(&dcount,&dcountGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    porosity=double(dcountGlobal)/volume;
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
    
    return 1.0-porosity;
}

inline double calculate_porosity(char * id, int Nx, int Ny, int Nz, Domain &Dm, int rank) {
    int nprocx,nprocy,nprocz;
    
    nprocx = Dm.nprocx();
    nprocy = Dm.nprocy();
    nprocz = Dm.nprocz();
    size_t count = 0;
    size_t totalGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                size_t n = k*Nx*Ny+j*Nx+i;
                if (id[n] == 2 || id[n] == 1)  count++;
            }
        }
    }
    
    // total Global is the number of nodes in the pore-space
    MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,Dm.Comm);
    double volume=double(nprocx*nprocy*nprocz)*double(Nx-2)*double(Ny-2)*double(Nz-2);
    double porosity=double(totalGlobal)/volume;
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    double dcount = 0;
    double dcountGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                for (int p=0;p<8;p++){
                    int n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                    if (id[n] == 2 || id[n] == 1)  {  dcount += 0.125; }
                }
            }
        }
    }
    MPI_Allreduce(&dcount,&dcountGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    porosity=double(dcountGlobal)/volume;
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
    
    return porosity;
}

inline int populate(int ijk /* regular coordinates */,
                    int m   /* The neighboring solid site */,
                    int dir /* lattice direction */,
                    double * LIBBqA, double * LIBBqBC, double * LIBBqD, double * qDistances,
                    int N   /* total sites for regular layout */
) {
    int count = 0;
    double A,BC,D;
    A = BC = D = 0; // Initialize all weights to zero
    double q = qDistances[ijk+dir*N];
    if (m == 0) { /* if neighbor is a solid site */
        if (q >= 0.5) {
            A = 1.-1./(2*q);
            BC = 1./(2*q);
            D = 0;
        }
        if (q < 0.5) {
            A = 0;
            BC = 2*q;
            D = (1.-2*q);
            
        }
        LIBBqA[ijk+dir*N] = A;
        LIBBqBC[ijk+dir*N] = BC;
        LIBBqD[ijk+dir*N] = D;
        count++;
    }
    else { /* if neighbor is a fluid site */
        LIBBqA[ijk+dir*N] = 0;
        LIBBqBC[ijk+dir*N] = 1;
        LIBBqD[ijk+dir*N] = 0;
        qDistances[ijk+dir*N] = 1;
    }
    return count;
}

inline int LIBB_Init(int * Map,
                     char *id,
                     int strideY,
                     int strideZ,
                     int start,
                     int finish,
                     int N, double * qDistances, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {
    
    int nn,ijk,m;
    size_t count = 0;
    for (int n=start; n<finish; n++){
        // Get the 1D index based on regular data layout
        ijk = Map[n]; // standard coordinates
        //........................................................................
        nn = ijk-1;                            // neighbor index (get convention)
        m = id[nn];                        // get neighbor for id - 2
        count += populate(ijk,m,1,LIBBqA,LIBBqBC,LIBBqD,qDistances,N); // contains the bug
        //........................................................................
        nn = ijk+1;                            // neighbor index (get convention)
        m = id[nn];                        // get neighbor for id - 1
        count += populate(ijk,m,0,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideY;                            // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 4
        count += populate(ijk,m,3,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideY;                            // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 3
        count += populate(ijk,m,2,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideZ;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 6
        count += populate(ijk,m,5,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideZ;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 5
        count += populate(ijk,m,4,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideY-1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 8
        count += populate(ijk,m,7,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideY+1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 7
        count += populate(ijk,m,6,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideY-1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 10
        count += populate(ijk,m,9,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideY+1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 9
        count += populate(ijk,m,8,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideZ-1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 12
        count += populate(ijk,m,11,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideZ+1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 11
        count += populate(ijk,m,10,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideZ-1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 14
        count += populate(ijk,m,13,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideZ+1;                        // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 13
        count += populate(ijk,m,12,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideZ-strideY;                    // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 16
        count += populate(ijk,m,15,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideZ+strideY;                    // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 15
        count += populate(ijk,m,14,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk+strideZ-strideY;                    // neighbor index (get convention)
        m = id[nn];                   // get neighbor for id - 18
        count += populate(ijk,m,17,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        //........................................................................
        nn = ijk-strideZ+strideY;                    // neighbor index (get convention)
        m = id[nn];                    // get neighbor for id - 17
        count += populate(ijk,m,16,LIBBqA,LIBBqBC,LIBBqD,qDistances,N);
        
    }
    return count;
}

ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK,
                                     int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),rhoA(0),rhoB(0),alpha(0),beta(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM){}

ScaLBL_ColorModel::~ScaLBL_ColorModel(){}



void ScaLBL_ColorModel::AssignComponentLabels(double *phase){
    size_t NLABELS=0;
    char VALUE=0;
    double AFFINITY=0.f;
    auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );
    auto LabelList = color_db->getVector<char>( "ComponentLabels" );
    
    
    NLABELS=LabelList.size();
    if (NLABELS != AffinityList.size()){
        ERROR("Error: ComponentLabels and ComponentAffinity must be the same length! \n");
    }
    
    int label_count[NLABELS];
    int label_count_global[NLABELS];
    // Assign the labels
    
    for (int idx=0; idx<NLABELS; idx++) label_count[idx]=0;
    
    for (int k=0;k<Nz;k++){
        for (int j=0;j<Ny;j++){
            for (int i=0;i<Nx;i++){
                int n = k*Nx*Ny+j*Nx+i;
                VALUE=id[n];
                // Assign the affinity from the paired list
                for (unsigned int idx=0; idx < NLABELS; idx++){
                    //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
                    if (VALUE == LabelList[idx]){
                        AFFINITY=AffinityList[idx];
                        label_count[idx]++;
                        idx = NLABELS;
                        Mask->id[n] = 0; // set mask to zero since this is an immobile component
                    }
                }
                // fluid labels are reserved
                if (VALUE == 1) AFFINITY = 1.0;
                else if (VALUE == 2) AFFINITY = -1.0;
                phase[n] = AFFINITY;
            }
        }
    }
    // Set Dm to match Mask
    for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
    
    MPI_Allreduce(&label_count[0],&label_count_global[0],NLABELS,MPI_INT,MPI_SUM,Dm->Comm);
    
    if (rank==0){
        printf("Components labels: %lu \n",NLABELS);
        for (unsigned int idx=0; idx<NLABELS; idx++){
            VALUE=LabelList[idx];
            AFFINITY=AffinityList[idx];
            printf("   label=%i, affinity=%f\n\n",int(VALUE),AFFINITY);
        }
    }
}

void ScaLBL_ColorModel::ReadParams(string filename){
    // read the input database
    db = std::make_shared<Database>( filename );
    domain_db = db->getDatabase( "Domain" );
    color_db = db->getDatabase( "Color" );
    analysis_db = db->getDatabase( "Analysis" );
    
    // Color Model parameters
    timestepMax = color_db->getScalar<int>( "timestepMax" );
    analysis_interval = analysis_db->getScalar<int>( "analysis_interval" );
    tauA = color_db->getScalar<double>( "tauA" );
    tauB = color_db->getScalar<double>( "tauB" );
    rhoA = color_db->getScalar<double>( "rhoA" );
    rhoB = color_db->getScalar<double>( "rhoB" );
    Fx = color_db->getVector<double>( "F" )[0];
    Fy = color_db->getVector<double>( "F" )[1];
    Fz = color_db->getVector<double>( "F" )[2];
    alpha = color_db->getScalar<double>( "alpha" );
    beta = color_db->getScalar<double>( "beta" );
    Restart = color_db->getScalar<bool>( "Restart" );
    din = color_db->getScalar<double>( "din" );
    dout = color_db->getScalar<double>( "dout" );
    flux = color_db->getScalar<double>( "flux" );
    method = color_db->getScalar<int>("method"); // method = 1=BGK, 2=MRT, 3=LIBB
    updateQ = color_db->getScalar<bool>("updateQ");
    inletA=1.f;
    inletB=0.f;
    outletA=0.f;
    outletB=1.f;
    
    // Read domain parameters
    auto L = domain_db->getVector<double>( "L" );
    auto size = domain_db->getVector<int>( "n" );
    auto nproc = domain_db->getVector<int>( "nproc" );
    BoundaryCondition = domain_db->getScalar<int>( "BC" );
    legacySDs = domain_db->getScalar<bool>("legacySDs");
    Nx = size[0];
    Ny = size[1];
    Nz = size[2];
    Lx = L[0];
    Ly = L[1];
    Lz = L[2];
    nprocx = nproc[0];
    nprocy = nproc[1];
    nprocz = nproc[2];
    
    Berea = domain_db->getScalar<bool>( "Berea" );
    
    string susr = domain_db->getScalar<string>( "usr" );
    int fs = domain_db->getScalar<int>( "filesystem" );
    
    if (rank == 0) {
        printf("INPUT PARAMETERS:\n");
        //printf("subdomain size: Nx=%d Ny=%d Nz=%d\n",Nx,Ny,Nz);
        printf("nprocx=%d nprocy=%d nprocz=%d\n", nprocx, nprocy, nprocz);
        printf("Lx=%-5.5E Ly=%-5.5E Lz=%-5.5E\n",Lx,Ly,Lz);
        printf("BC=%d\n",BoundaryCondition);
        printf("tauA=%-5.5E tauB=%-5.5E\n",tauA,tauB); //rhoA=%-5.5E rhoB=%-5.5E\n",tauA,tauB,rhoA,rhoB);
        printf("Fx=%-5.5E Fy=%-5.5E Fz=%-5.5E\n",Fx,Fy,Fz);
        /*
         printf("alpha=%-5.5E beta=%-5.5E\n",alpha,beta);
         printf("din=%-5.5E dout=%-5.5E flux=%-5.5E\n",din,dout,flux);
         printf("filesystem=%d usr=%s\n",fs,susr.c_str());
         */
        printf("timestepMax=%d\n",timestepMax);
        printf("\n");
        printf("LBM METHOD:\n");
        if (method == 1 ) cout << "NON-AA MRT-LIBB" << endl << endl;
        if (method == 2 ) cout << "AA MRT-LIBB" << endl << endl;
    }
    if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)
}

void ScaLBL_ColorModel::SetDomain(){
    
    Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
    Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));    // mask domain removes immobile phases
    Nx+=2; Ny+=2; Nz += 2;
    N = Nx*Ny*Nz;
    id = new char [N];
    for (size_t i=0; i<N; i++) Dm->id[i] = 1;               // initialize this way
    std::array<size_t,3> subdivide = { 1, 1, 1 };
    if ( analysis_db->keyExists( "subdivide" ) ) {
        auto tmp = analysis_db->getVector<size_t>( "subdivide" );
        subdivide = { tmp[0], tmp[1], tmp[2] };
    }
    Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm,subdivide) ); // TwoPhase analysis object
    MPI_Barrier(comm);
    
    Dm->CommInit(); // Sets up communication
    MPI_Barrier(comm);
    rank = Dm->rank();
    
    Averages->method = method;
    
    
    Mask->ReadIDs(); // Read ID files in this function
    for (int i=0; i<N; i++) id[i] = Mask->id[i];  // save what was read
    
}

void ScaLBL_ColorModel::ReadInput() {

    size_t readID;
    
    if (legacySDs == false) {
        if (rank == 0) cout << "Using Ray-traced signed distance" << endl;
        /* Ray traced signed distance */
        char LocalRankFilename[40];
        /* New ray tracing distances: 'signed distance' */
        sprintf(LocalRankFilename,"RTD.%05i",rank);
        SignDist = new double[N];
        FILE *RayTraceDistance = fopen(LocalRankFilename,"rb");
        readID=fread(SignDist,8,N,RayTraceDistance);
        if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading RayTraceDistance \n");
        fclose(RayTraceDistance);
        for (size_t i=0; i < N; i++) Averages->SDs(i) = SignDist[i];  // (1) Save RTD data into Averages->SDs
        
        int * spheresPerRankArray; spheresPerRankArray = new int[nprocx*nprocy*nprocz];
        if (rank == 0) {
            sprintf(LocalRankFilename,"NSpheresPerRankArray.%05i",rank);
            FILE *NSPRA = fopen(LocalRankFilename,"rb");
            readID = fread(spheresPerRankArray,4,nprocx*nprocy*nprocz,NSPRA);
            if (readID != size_t(nprocx*nprocy*nprocz)) printf("Error reading NSPRA \n");
            fclose(NSPRA);
            
        }
        
        int NspheresOnRank = 0;
        MPI_Scatter(spheresPerRankArray, 1,  MPI_INT, &NspheresOnRank, 1, MPI_INT,0, comm);
        
        // std::cout << rank << ") NspheresOnRank=" << NspheresOnRank << std::endl;
        Averages->NSpheresOnRank = NspheresOnRank;
        Averages->InitializeRankSignDistance();
        
        
        double * ranksds; ranksds = new double[N*NspheresOnRank];
        
        sprintf(LocalRankFilename,"RankRTD.%05i",rank);
        FILE *RankRayTraceDistance = fopen(LocalRankFilename,"rb");
        readID=fread(ranksds,8,N*NspheresOnRank,RankRayTraceDistance);
        if (readID != size_t(N*NspheresOnRank)) printf("lbpm_segmented_pp: Error reading RankRayTraceDistance \n");
        fclose(RankRayTraceDistance);
        for (size_t i=0; i < N*NspheresOnRank; i++) Averages->rankSDs(i) = ranksds[i];  // (1) Save RankRTD data into Averages->rankSDs
        delete[] ranksds;
        // std::cout << "Finished reading and passing rankSDs to Averages class" << std::endl;
        // Averages->RankSignDistance();

        
    }
    else {
        if (rank == 0) cout << "Using legacy ID-based signed distance" << endl;
        // Generate the signed distance map
        // Initialize the domain and communication
        Array<char> id_solid(Nx,Ny,Nz);
        int count = 0;
        // Solve for the position of the solid phase
        for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                for (int i=0;i<Nx;i++){
                    int n = k*Nx*Ny+j*Nx+i;
                    // Initialize the solid phase
                    if (Mask->id[n] > 0)    id_solid(i,j,k) = 1;
                    else                    id_solid(i,j,k) = 0;
                }
            }
        }
        // Initialize the signed distance function
        for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    // Initialize distance to +/- 1
                    Averages->SDs( i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
                }
            }
        }
        //    MeanFilter(Averages->SDs);
        if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
        CalcDist(Averages->SDs,id_solid,*Mask);
        
    }
    
    if (Berea == true) {
        char LocalRankString[8];
        char LocalRankFilename1[40];
        sprintf(LocalRankString,"%05d",rank);
        
        char LocalRankString2[8];
        char LocalRankFilename2[40];
        sprintf(LocalRankString2,"%05d",rank);
        
        // Write to disk
        sprintf(LocalRankFilename1,"%s%s","SignDist.",LocalRankString);
        FILE *SDFILE = fopen(LocalRankFilename1,"wb");
        if (SDFILE==NULL) ERROR("Error opening file: SDFILE.xxxxx");
        fwrite(Averages->SDs.data(),8,N,SDFILE);
        fclose(SDFILE);
        
        sprintf(LocalRankFilename2,"%s%s","RTD.",LocalRankString2);
        FILE *SDFILE2 = fopen(LocalRankFilename2,"wb");
        if (SDFILE2==NULL) ERROR("Error opening file: RTD.xxxxx");
        fwrite(Averages->SDs.data(),8,N,SDFILE2);
        fclose(SDFILE2);
        
    }
    
    
} // ray trace distance data is read in here


void ScaLBL_ColorModel::InitializeSolidVolumeFraction() {
    if (Restart == false) {
        Averages->PhaseInit(af); // Computes VFmask (depends on SDs and ID)
        for (size_t n = 0; n < N; n++) PhaseLabel[n] = Averages->VFmask(n);
        ScaLBL_CopyToDevice(vfmask, PhaseLabel, sizeof(double)*N);
        ScaLBL_Comm_Regular->SendHalo(vfmask);
        ScaLBL_Comm_Regular->RecvHalo(vfmask);
        
        
    
        Ni = ScaLBL_Comm->MemoryOptimizedInactiveLayout(InactiveMap, Dm->id, PhaseLabel, Ni);
        ScaLBL_AllocateDeviceMemory((void **) &dvcInactiveMap,     sizeof(int)*Ni);
        
        // Check InactiveMap
        TmpMap=new int[Ni];
        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    int idx=InactiveMap(i,j,k);
                    if (!(idx < 0)) TmpMap[idx] = k*Nx*Ny+j*Nx+i;
                }
            }
        }
        // check that TmpMap is valid
        for (int idx=0; idx<ScaLBL_Comm->LastInactiveExterior(); idx++){
            int n = TmpMap[idx];
            if (n > Nx*Ny*Nz){
                printf("Inactive Site Bad value! idx=%i \n",n);
                TmpMap[idx] = Nx*Ny*Nz-1;
            }
        }
        for (int idx=ScaLBL_Comm->FirstInactiveInterior(); idx<ScaLBL_Comm->LastInactiveInterior(); idx++){
            int n = TmpMap[idx];
            if (n > Nx*Ny*Nz){
                printf("Inactive Site Bad value! idx=%i \n",n);
                TmpMap[idx] = Nx*Ny*Nz-1;
            }
        }
        
        ScaLBL_CopyToDevice(dvcInactiveMap, TmpMap, sizeof(int)*Ni);
        delete[] TmpMap;
        
        if (rank == 0 ) std::cout << " dvcInactiveMap on Device" << std::endl;
        
        char LocalRankString[8];
        char LocalRankFilename1[40];
        sprintf(LocalRankString,"%05d",rank);
        
        // Write to disk
        sprintf(LocalRankFilename1,"%s%s","VFmask.",LocalRankString);
        FILE *VFFILE = fopen(LocalRankFilename1,"wb");
        if (VFFILE==NULL) ERROR("Error opening file: VFmask.xxxxx");
        fwrite(PhaseLabel,8,N,VFFILE);
        fclose(VFFILE);
        /* Save VFmask to disk*/
    } else {
        /* Read VFmask from disk */
        char LocalRankString[8];
        char LocalRankFilename1[40];
        sprintf(LocalRankString,"%05d",rank);
        size_t read_file;
        
        // Read from disk
        sprintf(LocalRankFilename1,"VFmask.%05i",rank);
        FILE *VFFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(PhaseLabel,8,N,VFFILE);
        if (read_file != size_t(N)) printf("Error reading VFFILE \n");
        fclose(VFFILE);
    }
    ScaLBL_CopyToDevice(vfmask, PhaseLabel, sizeof(double)*N);
    for (size_t n = 0; n < N; n++) Averages->VFmask(n) = PhaseLabel[n];
}

void ScaLBL_ColorModel::InitializeInterpolationArrays() {
    if (Restart == false ) {
        // Populate the q distances
        
        if (Berea == false) {
            sprintf(LocalRankFilename,"q.%05i",rank);
            FILE *DIST = fopen(LocalRankFilename,"rb");
            size_t ReadSignDist;
            ReadSignDist=fread(qDistances,8,NqN,DIST);
            if (ReadSignDist != size_t(NqN)) fclose(DIST);
            
            NormalizeQdistances(Dm->id, N, qDistances, libbqA,libbqBC, libbqD);
            
            // Populates the libbq values based upon qDistances
            size_t count = LIBB_Init(TmpMap, Dm->id, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(),
                                     ScaLBL_Comm->LastInterior(), N, qDistances, libbqA, libbqBC, libbqD);
            size_t count2 = LIBB_Init(TmpMap, Dm->id, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(),
                                      N, qDistances, libbqA, libbqBC, libbqD);
            
            for (size_t n = 0; n < NqN; n++) Averages->qDistances(n) = qDistances[n];
        
        }
        else {
            
            for (int dir = 0; dir < 18; dir++) {
                for (int ijk = 0; ijk < N; ijk++) {
                    LIBBqA[ijk+dir*N] = 0;
                    LIBBqBC[ijk+dir*N] = 1;
                    LIBBqD[ijk+dir*N] = 0;
                }
            }
            
        }
        
        char LocalRankString[8];
        char LocalRankFilename1[40];
        sprintf(LocalRankString,"%05d",rank);
        
        // Write to disk
        sprintf(LocalRankFilename1,"%s%s","libbqA.",LocalRankString);
        FILE *AFILE = fopen(LocalRankFilename1,"wb");
        if (AFILE==NULL) ERROR("Error opening file: libbqA.xxxxx");
        fwrite(libbqA,8,18*N,AFILE);
        fclose(AFILE);
        
        sprintf(LocalRankFilename1,"%s%s","libbqBC.",LocalRankString);
        FILE *BCFILE = fopen(LocalRankFilename1,"wb");
        if (BCFILE==NULL) ERROR("Error opening file: libbqBC.xxxxx");
        fwrite(libbqBC,8,18*N,BCFILE);
        fclose(BCFILE);
        
        sprintf(LocalRankFilename1,"%s%s","libbqD.",LocalRankString);
        FILE *DFILE = fopen(LocalRankFilename1,"wb");
        if (DFILE==NULL) ERROR("Error opening file: libbqD.xxxxx");
        fwrite(libbqD,8,18*N,DFILE);
        fclose(DFILE);
        
        sprintf(LocalRankFilename1,"%s%s","q2.",LocalRankString);
        FILE *QFILE = fopen(LocalRankFilename1,"wb");
        if (QFILE==NULL) ERROR("Error opening file: Q.xxxxx");
        fwrite(qDistances,8,18*N,QFILE);
        fclose(QFILE);
        
                
        
        
        
        
    } else {
        char LocalRankString[8];
        char LocalRankFilename1[40];
        sprintf(LocalRankString,"%05d",rank);
        
        size_t read_file;
        
        // Read from disk
        sprintf(LocalRankFilename1,"libbqA.%05i",rank);
        FILE *AFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(libbqA,8,18*N,AFILE);
        if (read_file != size_t(18*N)) printf("Error reading AFILE \n");
        fclose(AFILE);
        
        sprintf(LocalRankFilename1,"libbqBC.%05i",rank);
        FILE *BCFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(libbqBC,8,18*N,BCFILE);
        if (read_file != size_t(18*N)) printf("Error reading BCFILE \n");
        fclose(BCFILE);
        
        sprintf(LocalRankFilename1,"libbqD.%05i",rank);
        FILE *DFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(libbqD,8,18*N,DFILE);
        if (read_file != size_t(18*N)) printf("Error reading DFILE \n");
        fclose(DFILE);
        
        sprintf(LocalRankFilename1,"q2.%05i",rank);
        FILE *QFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(qDistances,8,18*N,QFILE);
        if (read_file != size_t(18*N)) printf("Error reading QFILE \n");
        fclose(QFILE);
    }
    
    
    
    ScaLBL_CopyToDevice(LIBBqA, libbqA, libbSize);
    ScaLBL_CopyToDevice(LIBBqBC, libbqBC, libbSize);
    ScaLBL_CopyToDevice(LIBBqD, libbqD, libbSize);
    for (size_t n = 0; n < NqN; n++) Averages->qDistances(n) = qDistances[n];
}





void ScaLBL_ColorModel::RestartDistances_pp() {

        char LocalRankString[8];
        char LocalRankFilename1[40];
        sprintf(LocalRankString,"%05d",rank);
    
        size_t read_file;
    
        // Available from RayTrace
        sprintf(LocalRankFilename1,"q.%05i",rank);
        FILE *DIST = fopen(LocalRankFilename1,"rb");
        read_file = fread(qDistances,8,18*N,DIST);
        if (read_file != size_t(18*N)) printf("Error reading DIST \n");
        fclose(DIST);
        std::cout << "Finished reading q distances" << std::endl;
        MPI_Barrier(comm);
    
        sprintf(LocalRankFilename1,"libbqA.%05i",rank);
        FILE *AFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(libbqA,8,18*N,AFILE);
        if (read_file != size_t(18*N)) printf("Error reading AFILE \n");
        fclose(AFILE);
        std::cout << "Finished reading A file" << std::endl;
        MPI_Barrier(comm);
    
        sprintf(LocalRankFilename1,"libbqBC.%05i",rank);
        FILE *BCFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(libbqBC,8,18*N,BCFILE);
        if (read_file != size_t(18*N)) printf("Error reading BCFILE \n");
        fclose(BCFILE);
        std::cout << "Finished reading BC file" << std::endl;
        MPI_Barrier(comm);
    
    
        sprintf(LocalRankFilename1,"libbqD.%05i",rank);
        FILE *DFILE = fopen(LocalRankFilename1,"rb");
        read_file = fread(libbqD,8,18*N,DFILE);
        if (read_file != size_t(18*N)) printf("Error reading DFILE \n");
        fclose(DFILE);
        std::cout << "Finished reading D file" << std::endl;
        MPI_Barrier(comm);
    
    
        NormalizeQdistances(Dm->id, N, qDistances, libbqA,libbqBC, libbqD);
        std::cout << "Finished NormalizeQdistances" << std::endl;
        MPI_Barrier(comm);
        
        // Populates the libbq values based upon qDistances
        size_t count = LIBB_Init(TmpMap, Dm->id, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(),
                                 ScaLBL_Comm->LastInterior(), N, qDistances, libbqA, libbqBC, libbqD);
        size_t count2 = LIBB_Init(TmpMap, Dm->id, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(),
                                  N, qDistances, libbqA, libbqBC, libbqD);
        std::cout << "Finished LIBB_Init" << std::endl;
        MPI_Barrier(comm);
    
    
        
    
        sprintf(LocalRankFilename1,"%s%s","q2.",LocalRankString);
        FILE *QFILE = fopen(LocalRankFilename1,"wb");
        if (QFILE==NULL) ERROR("Error opening file: Q.xxxxx");
        fwrite(qDistances,8,18*N,QFILE);
        fclose(QFILE);
        std::cout << "Finished writing q2 files" << std::endl;
    
    
}


void ScaLBL_ColorModel::ReadRestartFiles() {
    /* Overwrite things if restart is true */
    if (Restart == true){
        if (rank==0){
            printf("Reading restart file! \n");
            ifstream restart("Restart.txt");
            if (restart.is_open()){
                restart  >> timestep;
                printf("Restarting from timestep %i \n",timestep);
            }
            else{
                printf("WARNING:No Restart.txt file, setting timestep=0 \n");
                timestep=0;
            }
        }
        
        auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );
        affinity = AffinityList[0]; 
        MPI_Bcast(&affinity,1,MPI_DOUBLE,0,comm);
        
        Fx = color_db->getVector<double>( "F" )[0];
        Fy = color_db->getVector<double>( "F" )[1];
        Fz = color_db->getVector<double>( "F" )[2];
        MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
        MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
        MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
        
        MPI_Bcast(&timestep,1,MPI_INT,0,comm);
        // Read in the restart file to CPU buffers
        double *cPhi, *cfq, *cDenA, *cDenB, *cAq, *cBq;
        
        cAq = new double[7*Np];
        cBq = new double[7*Np];
        cfq = new double[19*Np];
        cDenA = new double[N];
        cDenB = new double[N];
        cPhi = new double[N];
        
        char LocalRestartFile[40];
        sprintf(LocalRestartFile,"Restart.%05i",rank);
        ifstream File(LocalRestartFile,ios::binary);
        double value;
        
        for (int n=0; n < Np; n++){
            for (int q=0; q<7; q++) {
                size_t element = q*Np+n;
                File.read((char*) &value, sizeof(value));
                cAq[element] = value;
            }
        }
        for (int n=0; n < Np; n++){
            for (int q=0; q<7; q++) {
                size_t element = q*Np+n;
                File.read((char*) &value, sizeof(value));
                cBq[element] = value;
            }
        }
        for (int n=0; n<Np; n++){
            for (int q=0; q<19; q++){
                size_t element = q*Np+n;
                File.read((char*) &value, sizeof(value));
                cfq[element] = value;
            }
        }
        for (int n=0; n < N; n++){
            File.read((char*) &value, sizeof(value));
            cDenA[n] = value;
        }
        for (int n=0; n < N; n++){
            File.read((char*) &value, sizeof(value));
            cDenB[n] = value;
        }
        for (int n=0; n<N; n++){
            File.read((char*) &value, sizeof(value));
            cPhi[n] = value;
        }
        File.close();
        
        
        
    
        for (size_t n = 0; n < N; n++) {
            if (Dm->id[n] == 0) cPhi[n] = affinity; // Updates the phase affinity
        }
        
        if (rank==0){
            printf("Solid phase affinity=%f\n",affinity);
        }
        
        if (rank==0){
            printf("Fx=%f Fy=%f Fz=%f\n",Fx,Fy,Fz);
        }
            
        
        
        // Copy the restart data to the GPU
        ScaLBL_CopyToDevice(DenA,cDenA,N*sizeof(double));
        ScaLBL_CopyToDevice(DenB,cDenB,N*sizeof(double));
        
        ScaLBL_CopyToDevice(fq,cfq,19*Np*sizeof(double));
        ScaLBL_CopyToDevice(Phi,cPhi,N*sizeof(double));
        ScaLBL_CopyToDevice(Aq,cAq,7*Np*sizeof(double));
        ScaLBL_CopyToDevice(Bq,cBq,7*Np*sizeof(double));
        ScaLBL_DeviceBarrier();
        
        MPI_Barrier(comm);
        delete[] cDenA;
        delete[] cDenB;
        delete[] cfq;
        delete[] cPhi;
        delete[] cAq;
        delete[] cBq;
    }
    
    if (BoundaryCondition >0 ){
        if (Dm->kproc()==0){
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
        }
        if (Dm->kproc() == nprocz-1){
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
        }
        
        
    // Correct the coefficient arrays
        
        int i,j,k;
        
        if (Dm->kproc() == 0)
        {
            for (int dir = 0; dir < 18; dir++) {
                for (int ijk = 0; ijk < N; ijk++) {
                    k = ijk/(Nx*Ny);
                    j = (ijk-Nx*Ny*k)/Nx;
                    i = ijk-Nx*Ny*k-Nx*j;
                    if (k < 3) {
                        LIBBqA[ijk+dir*N] = 0;
                        LIBBqBC[ijk+dir*N] = 1;
                        LIBBqD[ijk+dir*N] = 0;
                    }
                }
            }
        }
        
        if (Dm->kproc() == nprocz-1)
        {
            for (int dir = 0; dir < 18; dir++) {
                for (int ijk = 0; ijk < N; ijk++) {
                    k = ijk/(Nx*Ny);
                    j = (ijk-Nx*Ny*k)/Nx;
                    i = ijk-Nx*Ny*k-Nx*j;
                    if (k > Nz-4) {
                        LIBBqA[ijk+dir*N] = 0;
                        LIBBqBC[ijk+dir*N] = 1;
                        LIBBqD[ijk+dir*N] = 0;
                    }
                }
            }
        }
        
        
        
    } // BC > 0
}

void ScaLBL_ColorModel::UpdateMapsAndLists() {
    // Check Map:
    TmpMap=new int[Np];
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                int idx=Map(i,j,k);
                if (!(idx < 0)) TmpMap[idx] = k*Nx*Ny+j*Nx+i;
            }
        }
    }
    // check that TmpMap is valid
    for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
        int n = TmpMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Bad value! idx=%i \n",n);
            TmpMap[idx] = Nx*Ny*Nz-1;
        }
    }
    for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
        int n = TmpMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Bad value! idx=%i \n",n);
            TmpMap[idx] = Nx*Ny*Nz-1;
        }
    }
    ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    ScaLBL_CopyToDevice(InterpolationList, interpolationList, neighborSize);
}

void ScaLBL_ColorModel::Create() {
    //.........................................................
    // Initialize communication structures in averaging domain
    for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
    Mask->CommInit();
    Np=Mask->PoreCount();
    ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
    ScaLBL_Comm_Regular  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
    
    size_t Npad=(Np/16 + 2)*16;
    Map.resize(Nx,Ny,Nz); Map.fill(-2);
    InactiveMap.resize(Nx,Ny,Nz); InactiveMap.fill(-2);
    
    
    /* Fill in the qDistances arrays */
    Nq = 18;
    NqN= Nq*N;
    NqNpad = Nq*Npad;
    
    
    /* local CPU registers */
    qDistances= new double[NqN];
    neighborList= new int[NqNpad];
    interpolationList= new int[NqNpad];
    libbqA= new double[NqN];
    libbqBC= new double[NqN];
    libbqD= new double[NqN];
    
    
    for (size_t i = 0; i < NqN; i++) { libbqA[i] = 0;  libbqBC[i] = 0;  libbqD[i] = 0;  }
    
    // Updated Np
    Np = ScaLBL_Comm->MemoryOptimizedLayoutAA_LIBB(Map,neighborList,interpolationList,Mask->id,Np);
    
    NqNp = Nq*Np;
    dist_mem_size = Np*sizeof(double);
    neighborSize=NqNp*sizeof(int);
    libbSize=NqN*sizeof(double);
    signSize = N*sizeof(double);
    qsize = NqN*sizeof(double);
    
    
    gradphix = new double[N];
    gradphiy = new double[N];
    gradphiz = new double[N];
    for (size_t i = 0; i < N; i++) { gradphix[i] = 0;  gradphiy[i] = 0;  gradphiz[i] = 0;  }
    
    //......... GPU STATE SPACE ............................................
    ScaLBL_AllocateDeviceMemory((void **) &NeighborList,      neighborSize);
    ScaLBL_AllocateDeviceMemory((void **) &InterpolationList, neighborSize);
    ScaLBL_AllocateDeviceMemory((void **) &dvcMap,     sizeof(int)*Np);
    ScaLBL_AllocateDeviceMemory((void **) &fq,         19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &fq2,        19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &savedfq,    19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Aq,         7*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Bq,         7*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Phi,          sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &dvcSignDist,  sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Pressure,     sizeof(double)*Np);
    ScaLBL_AllocateDeviceMemory((void **) &Velocity,   3*sizeof(double)*Np);
    ScaLBL_AllocateDeviceMemory((void **) &LIBBqA,     libbSize);
    ScaLBL_AllocateDeviceMemory((void **) &LIBBqBC,    libbSize);
    ScaLBL_AllocateDeviceMemory((void **) &LIBBqD,     libbSize);
    ScaLBL_AllocateDeviceMemory((void **) &DenA,       sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &DenB,       sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &vfmask,     sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &ID,         sizeof(char)*N);
    
    ScaLBL_AllocateDeviceMemory((void **) &GradPhiX,       sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &GradPhiY,       sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &GradPhiZ,       sizeof(double)*N);
    

    //.....................................................................
    ScaLBL_CopyToDevice(GradPhiX, gradphix, N*sizeof(double));
    ScaLBL_CopyToDevice(GradPhiY, gradphiy, N*sizeof(double));
    ScaLBL_CopyToDevice(GradPhiZ, gradphiz, N*sizeof(double));
    
    
    ScaLBL_CopyToDevice(ID, Dm->id, N*sizeof(char));
    for (int i = 0; i < N; i++) Averages->ID(i) = Dm->id[i];
    fflush(stdout);
    
    UpdateMapsAndLists();
    
    if (!updateQ) {
        InitializeInterpolationArrays();
    }
    
    if (updateQ) {
        RestartDistances_pp();
    }
    
    /* Memory management */
    delete[] qDistances;
    delete[] neighborList;
    delete[] interpolationList;
    delete[] libbqA;
    delete[] libbqBC;
    delete[] libbqD;
    delete[] TmpMap;
}


void ScaLBL_ColorModel::InitializePhaseAndDensityFields() {
    if (Restart == false) {
        for (size_t n = 0; n < N; n++) {
            if (Dm->id[n] == 0) PhaseLabel[n] = af;
            else if (Dm->id[n] == 1) PhaseLabel[n] = 1.0;
            else if (Dm->id[n] == 2) PhaseLabel[n] = -1.0;
        }
        
        /* Save for visualization of initial state */
        for (size_t n = 0; n < N; n++) Averages->Phase(n)  = PhaseLabel[n];
        ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
        
        ScaLBL_PhaseField_Init_D3Q7(dvcMap,Phi, DenA, DenB, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np, N);
        ScaLBL_PhaseField_Init_D3Q7(dvcMap,Phi, DenA, DenB, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N);
    }
    
    
}

void ScaLBL_ColorModel::ReadSignedDistanceField() {
    sprintf(LocalRankFilename,"SignDist.%05i",rank);
    SignDist = new double[N];
    
    FILE *SDist = fopen(LocalRankFilename,"rb");
    size_t readID=fread(SignDist,8,N,SDist);
    if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading SDist \n");
    fclose(SDist);
    for (size_t i=0; i < N; i++) Averages->SDs(i) = SignDist[i];
}

void ScaLBL_ColorModel::Initialize() {
    auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );
    af = AffinityList[0]; Averages->affinity = af;
    
    ScaLBL_D3Q19_Init(fq, Np);
    PhaseLabel = new double[N];
    InitializeSolidVolumeFraction();
    InitializePhaseAndDensityFields();
    ReadSignedDistanceField();
    
    ReadRestartFiles();
    
    
    if (BoundaryCondition >0 ){
        if (Dm->kproc()==0){
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
        }
        if (Dm->kproc() == nprocz-1){
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
        }
    }

    delete[] PhaseLabel;
}


std::string getLastLine(std::ifstream& in){
    std::string line;
    while (in >> std::ws && std::getline(in, line));
    
    return line;
}

void ScaLBL_ColorModel::Analyze() {
    
    
    if ( timestep % analysis_interval == 0) {
        auto CaThresh = analysis_db->getScalar<size_t>( "CaThreshold" );
        
        double progress = 0;
       
        double Ca = Averages->Ca_global;
        double Re = Averages->Re_global;
        double pw = Averages->pw_global;
        double pn = Averages->pn_global;
        double poros = Averages->porosity_global;
        double png = Averages->pn_global;
        double ca = Averages->efawns_global;
        double total_mass = Averages->mass_global;
        double awn = Averages->awn_global;
        double Jwn = Averages->Jwn_global;
        double Kwn = Averages->Kwn_global;
        double euler = Averages->euler_global;
        double nw_phase_vol = Averages->global_nw_phase_volume;
        double vwndnw = Averages->vnwdnw_global;
        double an = Averages->an_global;
        double Jwn2 = Averages->Jwn2_global;
        double Kwn2 = Averages->Kwn2_global;
        double nwOwnOx = Averages->nwOwnO_global(0);
        double nwOwnOy = Averages->nwOwnO_global(1);
        double nwOwnOz = Averages->nwOwnO_global(2);
        double nwOwnOx2 = Averages->nwOwnO_global2(0);
        double nwOwnOy2 = Averages->nwOwnO_global2(1);
        double nwOwnOz2 = Averages->nwOwnO_global2(2);
        double awn2 = Averages->awn2_global;
        double Eta = Averages->Eta_global;
        double igl_Jwn = Averages->igl_mean_curvature;
        double igl_Kwn = Averages->igl_Gaussian_curvature;
        
        

        
        
        progress = double(timestep) / double(timestepMax) * 100.0;
        if (rank == 0) {
            cout <<
            progress << std::setprecision(12) <<
            "," << Ca <<
            "," << Re <<
            "," << ca <<
            "," << awn <<
            "," << Jwn <<
            "," << Kwn <<
//            "," << vwndnw <<
            "," << an <<
//            "," << pn <
//            "," << pw <<
            "," << Jwn2 <<
            "," << Kwn2 <<
            "," << Eta <<
            "," << nwOwnOx <<
            "," << nwOwnOy <<
            "," << nwOwnOz <<
            "," << awn2 <<
            "," << nwOwnOx2 <<
            "," << nwOwnOy2 <<
            "," << nwOwnOz2 <<
            "," << euler <<
            "\r" << std::endl;
        }
    }
}

void ScaLBL_ColorModel::LIBBNonAAPhaseField_1(){
    ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq);
    ScaLBL_D3Q7_PhaseField_LIBB(NeighborList, dvcMap, Aq, Bq, DenA, DenB, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N, Nx, Nx*Ny);
    ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq);
    
    ScaLBL_DeviceBarrier();
    ScaLBL_D3Q7_PhaseField_LIBB(NeighborList, dvcMap, Aq, Bq, DenA, DenB, Phi, 0, ScaLBL_Comm->LastExterior(), Np, N, Nx, Nx*Ny);
    MPI_Barrier(comm);  ScaLBL_DeviceBarrier();
}

void ScaLBL_ColorModel::LIBBNonAAPhaseField_2(){
    ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq);
    ScaLBL_D3Q7_PhaseField_LIBB(NeighborList,dvcMap, Aq, Bq, DenA, DenB, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N, Nx, Nx*Ny);
    ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq);
    
    ScaLBL_DeviceBarrier();
    ScaLBL_D3Q7_PhaseField_LIBB(NeighborList,dvcMap, Aq, Bq, DenA, DenB, Phi, 0, ScaLBL_Comm->LastExterior(), Np, N, Nx, Nx*Ny);
    ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
}

void ScaLBL_ColorModel::LIBBNonAAStreamAndCollide_1(int timestep){
    ScaLBL_Comm_Regular->SendHalo(Phi);
    if (BoundaryCondition > 0) {
        ScaLBL_Comm->Color_BC_z_2(dvcMap, Phi, DenA, DenB, inletA, inletB);
        ScaLBL_Comm->Color_BC_Z_2(dvcMap, Phi, DenA, DenB, outletA, outletB);
    }
    ScaLBL_Comm->SendD3Q19AA(fq);
    ScaLBL_D3Q19_Color_LIBB(InterpolationList, NeighborList, dvcMap, fq, fq2, savedfq, Aq, Bq, DenA, DenB, Phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ);
    ScaLBL_Comm_Regular->RecvHalo(Phi);
    ScaLBL_Comm->RecvD3Q19AA(fq);
    
    if (BoundaryCondition == 3){
        ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
        ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
    }
    if (BoundaryCondition == 4){
        din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
        ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
    }
    
    ScaLBL_DeviceBarrier();
    ScaLBL_D3Q19_Color_LIBB(InterpolationList, NeighborList, dvcMap, fq, fq2, savedfq, Aq, Bq, DenA, DenB, Phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ);
    ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
}

void ScaLBL_ColorModel::LIBBNonAAStreamAndCollide_2(int timestep){
    ScaLBL_Comm_Regular->SendHalo(Phi);
    if (BoundaryCondition > 0){
        ScaLBL_Comm->Color_BC_z_2(dvcMap, Phi, DenA, DenB, inletA, inletB);
        ScaLBL_Comm->Color_BC_Z_2(dvcMap, Phi, DenA, DenB, outletA, outletB);
    }
    ScaLBL_Comm->SendD3Q19AA(fq2);
    ScaLBL_D3Q19_Color_LIBB(InterpolationList, NeighborList, dvcMap, fq2, fq, savedfq, Aq, Bq, DenA, DenB, Phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD,GradPhiX, GradPhiY, GradPhiZ);
    ScaLBL_Comm_Regular->RecvHalo(Phi);
    ScaLBL_Comm->RecvD3Q19AA(fq2);
    
    if (BoundaryCondition == 3){
        ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq2, din, timestep);
        ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq2, dout, timestep);
    }
    if (BoundaryCondition == 4){
        din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq2, flux, timestep);
        ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq2, dout, timestep);
    }
    ScaLBL_DeviceBarrier();
    ScaLBL_D3Q19_Color_LIBB(InterpolationList, NeighborList, dvcMap, fq2, fq, savedfq, Aq, Bq, DenA, DenB, Phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ);
    ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
}

void ScaLBL_ColorModel::SaveMomentumDistribution_1() {
    save_state(fq,savedfq,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),19,Np);
    save_state(fq,savedfq,0, ScaLBL_Comm->LastExterior(),19,Np);
    ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
}

void ScaLBL_ColorModel::SaveMomentumDistribution_2() {
    save_state(fq2,savedfq,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),19,Np);
    save_state(fq2,savedfq,0, ScaLBL_Comm->LastExterior(),19,Np);
    ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
}

void ScaLBL_ColorModel::LIBBwithAAPhaseField_Odd() {
}

void ScaLBL_ColorModel::LIBBwithAAStreamAndCollideAndStream_Odd() {
}

void ScaLBL_ColorModel::LIBBwithAAPhaseField_Even() {
}

void ScaLBL_ColorModel::LIBBwithAACollide_Even() {
}

void ScaLBL_ColorModel::FillInactiveDensities() {
}






void ReadSpherePacking2(int nspheres,
                       double *List_cx, double *List_cy, double *List_cz, double *List_rad) {
    FILE *fid = fopen("pack.out","rb");
    char line[100];
    int count = 0;
    while ( !feof(fid) && fgets(line,100,fid)!=NULL ) {
        char* line2 = line;
        List_cx[count] = strtod(line2,&line2);
        List_cy[count] = strtod(line2,&line2);
        List_cz[count] = strtod(line2,&line2);
        List_rad[count] = strtod(line2,&line2);
        count++;
    }
    INSIST( count==nspheres, "Specified number of spheres is probably incorrect!" );
    // .............................................................
}



// Run Function

void ScaLBL_ColorModel::Run(){
    
    
    
  
    
    int nprocs=nprocx*nprocy*nprocz;
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
    
    int nspheres = domain_db->getScalar<int>( "nspheres");
    double *ccx,*ccy,*ccz,*r;
    ccx = new double[nspheres];
    ccy = new double[nspheres];
    ccz = new double[nspheres];
    r = new double[nspheres];
    if (Berea == false) {
        ReadSpherePacking2(nspheres,
                          ccx,ccy,ccz,r);
        double identical_sphere_radius = r[0];
        double sphere_diameter = 2.0*identical_sphere_radius;
        double avgMu = (0.5*(tauA + tauB)-0.5)/3.0;
        double avgRho = 0.5*(rhoA + rhoB);
        Averages->alpha = alpha;
        Averages->sphere_diameter = sphere_diameter;
        Averages->average_rho_input = avgRho;
        Averages->average_mu = avgMu;
        Averages->muA = (tauA-0.5)/3.0;
        Averages->muB = (tauB-0.5)/3.0;
    }
    
    
    PROFILE_START("Loop");
    bool Regular = false;
    runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, beta, Map );
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

    
    if (rank == 0)
    {
        cout << "\nLONG-TERM ITERATIONS:"<< endl;
        cout << "%" << ",Ca" << ",Re" << ",ca" << ",awn" << ",Jwn" << ",Kwn" << ",an"<< ",Jwn2" << ",Kwn2" << ",Eta" << ",nwOwnOx" << ",nwOwnOy" << ",nwOwnOz" << ",awn2" <<  ",nwOwnOx2" << ",nwOwnOy2" << ",nwOwnOz2" << ",euler" << "\r" << std::endl;
    }
    
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    while (timestep < timestepMax) {
        LIBBNonAAPhaseField_1();
        SaveMomentumDistribution_1();
        LIBBNonAAStreamAndCollide_1(timestep);
        timestep++;
        analysis.run4(timestep,*Averages,Phi,Pressure,Velocity,fq2,Aq,Bq,GradPhiX, GradPhiY, GradPhiZ, DenA, DenB, Np,Fx,Fy,Fz);  analysis.finish();
        Analyze();
        
        LIBBNonAAPhaseField_2();
        SaveMomentumDistribution_2();
        LIBBNonAAStreamAndCollide_2(timestep);
        timestep++;
        analysis.run4(timestep,*Averages,Phi,Pressure,Velocity,fq,Aq,Bq,GradPhiX, GradPhiY, GradPhiZ, DenA, DenB ,Np,Fx,Fy,Fz);  analysis.finish();
        Analyze();
    }
    
    
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    if (rank == 0) cout << "Finished lattice Boltzmann iterations.  Exiting simulator." << endl;
}
