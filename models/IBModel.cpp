///*
// Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
// 
// This file is part of the Open Porous Media project (OPM).
// OPM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// OPM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with OPM.  If not, see <http://www.gnu.org/licenses/>.
// */
///*
// color lattice boltzmann model
// */
//#include "models/IBModel.h"
//#include "analysis/distance.h"
//#include <vector>
//#include <string>
//
//
//inline void MeanFilter(Array<double> & tmp, DoubleArray &Mesh){
//
//    
//      double a = 1./3.;
//      double b = 1./18.;
//      double c = 1./36.;
//      double sum = 0;
//      size_t n, nn;
//      int i,j,k;
//        int Nx = Mesh.size(0);
//        int Ny = Mesh.size(1);
//        int Nz = Mesh.size(2);
//      size_t N = Nx*Ny*Nz;
//      for (n=0; n<N; n++){
//          sum = 0;
//          sum = a*tmp(n);
//
//
//          //.......Back out the 3-D indices for node n..............
//          k = n/(Nx*Ny);
//          j = (n-Nx*Ny*k)/Nx;
//          i = n-Nx*Ny*k-Nx*j;
//          //........................................................................
//          //........Get 1-D index for this thread....................
//          //        n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
//          //........................................................................
//          //                    COMPUTE THE COLOR GRADIENT
//          //........................................................................
//          //.................Read Phase Indicator Values............................
//          //........................................................................
//          nn = n-1;                            // neighbor index (get convention)
//          if (i-1<0)        nn += Nx;            // periodic BC along the x-boundary
//          sum += b*tmp(nn);                      // get neighbor for phi - 1
//          //........................................................................
//          nn = n+1;                            // neighbor index (get convention)
//          if (!(i+1<Nx))    nn -= Nx;            // periodic BC along the x-boundary
//          sum += b*tmp(nn);                     // get neighbor for phi - 2
//          //........................................................................
//          nn = n-Nx;                            // neighbor index (get convention)
//          if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
//          sum += b*tmp(nn);                   // get neighbor for phi - 3
//          //........................................................................
//          nn = n+Nx;                            // neighbor index (get convention)
//          if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
//          sum += b*tmp(nn);                   // get neighbor for phi - 4
//          //........................................................................
//          nn = n-Nx*Ny;                        // neighbor index (get convention)
//          if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += b*tmp(nn);                   // get neighbor for phi - 5
//          //........................................................................
//          nn = n+Nx*Ny;                        // neighbor index (get convention)
//          if (!(k+1<Nz))    nn -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += b*tmp(nn);                 // get neighbor for phi - 6
//          //........................................................................
//          nn = n-Nx-1;                        // neighbor index (get convention)
//          if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
//          if (j-1<0)            nn += Nx*Ny;    // Perioidic BC along the y-boundary
//          sum += c*tmp(nn);                    // get neighbor for phi - 7
//          //........................................................................
//          nn = n+Nx+1;                        // neighbor index (get convention)
//          if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
//          if (!(j+1<Ny))        nn -= Nx*Ny;    // Perioidic BC along the y-boundary
//          sum += c*tmp(nn);                   // get neighbor for phi - 8
//          //........................................................................
//          nn = n+Nx-1;                        // neighbor index (get convention)
//          if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
//          if (!(j+1<Ny))        nn -= Nx*Ny;    // Perioidic BC along the y-boundary
//          sum += c*tmp(nn);           // get neighbor for phi - 9
//          //........................................................................
//          nn = n-Nx+1;                        // neighbor index (get convention)
//          if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
//          if (j-1<0)            nn += Nx*Ny;    // Perioidic BC along the y-boundary
//          sum += c*tmp(nn);               // get neighbor for phi - 10
//          //........................................................................
//          nn = n-Nx*Ny-1;                        // neighbor index (get convention)
//          if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
//          if (k-1<0)            nn += Nx*Ny*Nz;    // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);              // get neighbor for phi - 11
//          //........................................................................
//          nn = n+Nx*Ny+1;                        // neighbor index (get convention)
//          if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
//          if (!(k+1<Nz))        nn -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);                 // get neighbor for phi - 12
//          //........................................................................
//          nn = n+Nx*Ny-1;                        // neighbor index (get convention)
//          if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
//          if (!(k+1<Nz))        nn -= Nx*Ny*Nz;    // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);            // get neighbor for phi - 13
//          //........................................................................
//          nn = n-Nx*Ny+1;                        // neighbor index (get convention)
//          if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
//          if (k-1<0)            nn += Nx*Ny*Nz;    // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);                  // get neighbor for phi - 14
//          //........................................................................
//          nn = n-Nx*Ny-Nx;                    // neighbor index (get convention)
//          if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
//          if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);                // get neighbor for phi - 15
//          //........................................................................
//          nn = n+Nx*Ny+Nx;                    // neighbor index (get convention)
//          if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
//          if (!(k+1<Nz))    nn -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);                   // get neighbor for phi - 16
//          //........................................................................
//          nn = n+Nx*Ny-Nx;                    // neighbor index (get convention)
//          if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
//          if (!(k+1<Nz))    nn -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);                    // get neighbor for phi - 17
//          //........................................................................
//          nn = n-Nx*Ny+Nx;                    // neighbor index (get convention)
//          if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
//          if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
//          sum += c*tmp(nn);                 // get neighbor for phi - 18
//    
//                     //........................................................................
//          Mesh(n) = sum;
//                     
//                     
//                      // Original method
//       // Mesh(n) = tmp(n);
//          }
//    
//    
//    
//    
//}
//
//ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK, int NP, MPI_Comm COMM):
//rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),rhoA(0),rhoB(0),alpha(0),beta(0),
//Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
//Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
//{
//    
//}
//ScaLBL_ColorModel::~ScaLBL_ColorModel(){
//    
//}
//
///*void ScaLBL_ColorModel::WriteCheckpoint(const char *FILENAME, const double *cPhi, const double *cfq, int Np)
// {
// int q,n;
// double value;
// ofstream File(FILENAME,ios::binary);
// for (n=0; n<Np; n++){
// // Write the two density values
// value = cPhi[n];
// File.write((char*) &value, sizeof(value));
// // Write the even distributions
// for (q=0; q<19; q++){
// value = cfq[q*Np+n];
// File.write((char*) &value, sizeof(value));
// }
// }
// File.close();
// 
// }
// 
// void ScaLBL_ColorModel::ReadCheckpoint(char *FILENAME, double *cPhi, double *cfq, int Np)
// {
// int q=0, n=0;
// double value=0;
// ifstream File(FILENAME,ios::binary);
// for (n=0; n<Np; n++){
// File.read((char*) &value, sizeof(value));
// cPhi[n] = value;
// // Read the distributions
// for (q=0; q<19; q++){
// File.read((char*) &value, sizeof(value));
// cfq[q*Np+n] = value;
// }
// }
// File.close();
// }
// */
//
//
//void ScaLBL_ColorModel::ReadParams(string filename){
//    // read the input database
//    db = std::make_shared<Database>( filename );
//    domain_db = db->getDatabase( "Domain" );
//    color_db = db->getDatabase( "Color" );
//    analysis_db = db->getDatabase( "Analysis" );
//    
//    // Color Model parameters
//    timestepMax = color_db->getScalar<int>( "timestepMax" );
//    tauA = color_db->getScalar<double>( "tauA" );
//    tauB = color_db->getScalar<double>( "tauB" );
//    rhoA = color_db->getScalar<double>( "rhoA" );
//    rhoB = color_db->getScalar<double>( "rhoB" );
//    Fx = color_db->getVector<double>( "F" )[0];
//    Fy = color_db->getVector<double>( "F" )[1];
//    Fz = color_db->getVector<double>( "F" )[2];
//    alpha = color_db->getScalar<double>( "alpha" );
//    beta = color_db->getScalar<double>( "beta" );
//    Restart = color_db->getScalar<bool>( "Restart" );
//    din = color_db->getScalar<double>( "din" );
//    dout = color_db->getScalar<double>( "dout" );
//    flux = color_db->getScalar<double>( "flux" );
//    inletA=1.f;
//    inletB=0.f;
//    outletA=0.f;
//    outletB=1.f;
//    
//    // Read domain parameters
//    auto L = domain_db->getVector<double>( "L" );
//    auto size = domain_db->getVector<int>( "n" );
//    auto nproc = domain_db->getVector<int>( "nproc" );
//    BoundaryCondition = domain_db->getScalar<int>( "BC" );
//    Nx = size[0];
//    Ny = size[1];
//    Nz = size[2];
//    Lx = L[0];
//    Ly = L[1];
//    Lz = L[2];
//    nprocx = nproc[0];
//    nprocy = nproc[1];
//    nprocz = nproc[2];
//    
//    string susr = domain_db->getScalar<string>( "usr" );
//    int fs = domain_db->getScalar<int>( "filesystem" );
//    
//    if (rank == 0) {
//        printf("============ INPUT PARAMETERS =================\n");
//        printf("subdomain size: Nx=%d Ny=%d Nz=%d\n",Nx,Ny,Nz);
//        printf("nprocx=%d nprocy=%d nprocz=%d\n", nprocx, nprocy, nprocz);
//        printf("Lx=%-5.5E Ly=%-5.5E Lz=%-5.5E\n",Lx,Ly,Lz);
//        printf("BC=%d\n",BoundaryCondition);
//        printf("tauA=%-5.5E tauB=%-5.5E rhoA=%-5.5E rhoB=%-5.5E\n",tauA,tauB,rhoA,rhoB);
//        printf("Fx=%-5.5E Fy=%-5.5E Fz=%-5.5E\n",Fx,Fy,Fz);
//        printf("alpha=%-5.5E beta=%-5.5E\n",alpha,beta);
//        printf("din=%-5.5E dout=%-5.5E flux=%-5.5E\n",din,dout,flux);
//        printf("filesystem=%d usr=%s\n",fs,susr.c_str());
//    }
//    
//    if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)
//    
//}
//void ScaLBL_ColorModel::SetDomain(){
//    Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
//    Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));    // mask domain removes immobile phases
//    Nx+=2; Ny+=2; Nz += 2;
//    N = Nx*Ny*Nz;
//    id = new char [N];
//    for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
//    std::array<size_t,3> subdivide = { 1, 1, 1 };
//    if ( analysis_db->keyExists( "subdivide" ) ) {
//        auto tmp = analysis_db->getVector<size_t>( "subdivide" );
//        subdivide = { tmp[0], tmp[1], tmp[2] };
//    }
//    Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm,subdivide) ); // TwoPhase analysis object
//    MPI_Barrier(comm);
//    Dm->CommInit();
//    MPI_Barrier(comm);
//    rank = Dm->rank();
//}
//
//
//
//bool ScaLBL_ColorModel::GridIndependence() {
//    
// if (Dm->rank()==0) printf("vector Grid Indep code\n"); 
//    
//    char *ReadID;
//    ReadID = new char[N];
//    
//    // Read the signed distance from file
//    sprintf(LocalRankFilename,"ID.%05i",rank);
//    FILE *DIST = fopen(LocalRankFilename,"rb");
//    size_t ReadSignDist;
//    ReadSignDist=fread(ReadID,1,N,DIST);
//    if (ReadSignDist != size_t(N)) printf("lbpm_grid_independence_pp: Error reading ID file (rank=%i)\n",rank);
//    fclose(DIST);
//    
//    int resolve;
//    if ( analysis_db->keyExists( "resolution" ) ) {
//        auto tmp = analysis_db->getScalar<int>( "resolution" );
//        resolve = int(tmp);
//    }
//    else {
//        resolve = 1;
//        
//    }
//    int resolution = resolve;
//    if (rank == 0) printf("resolution level=%d\n",resolution);
//    if (rank == 0) printf("nprocs=%d\n",nprocs);
//    double sum,sum_sat;
//    double porosity,saturation;
//    double sum_local=0.0;
//    double sum_sat_local=0.0;
//    double iVol_global = 0.0;
//    int refinedN = 0;
//    int Nx2, Ny2, Nz2;
//    
//    // Step 2
//    int Nx3, Ny3, Nz3, FinN;
//    
//    
//  //  resolution = 3;
//  //  for (int resolution = 1; resolution < 6; resolution++) {
//        
//        sum = 0;
//        sum_sat = 0;
//        sum_local = 0.0;
//        sum_sat_local = 0.0;
//        
//     
//        
//        // Step 1: Just increase the resolution and determine saturation and porosity: This gives the correct answer
//        // but the arrays are far too large - need to trim
//        Nx2 = resolution * (Nx) ;
//        Ny2 = resolution * (Ny) ;
//        Nz2 = resolution * (Nz) ;
//        //printf("%dx New size = %d\n",resolution,Nx2);
//        refinedN = Nx2*Ny2*Nz2;
//        //printf("%dx Media size = %d | ",resolution,Nx2);
//        
//        // Step 2
//        Nx3 = resolution * (Nx - 2) + 2;
//        Ny3 = resolution * (Ny - 2) + 2;
//        Nz3 = resolution * (Nz - 2) + 2;
//        FinN = Nx3*Ny3*Nz3;
//        //printf("%dx Fin size = %d\n",resolution,Nx3);
//   
//    vector<char> NewID(refinedN);
//    vector<char> TmpID(refinedN);
//    vector<char> FinID(FinN);
////        char *NewID;
////        NewID=new char[refinedN];
////        char *TmpID;
////        TmpID=new char[refinedN];
////        char *FinID;
////        FinID=new char[FinN];
//    
//    char *WriteID;
//    int WriteN = FinN;
//    WriteID=new char[WriteN];
//    
//        int k2, j2, i2;
//        // Step 1: Just increase the resolution and determine saturation and porosity: This gives the correct answer
//        // but the arrays are far too large - need to trim
//        for (int k=0;k<Nz;k++){
//            k2 = resolution*k;
//            for (int j=0;j<Ny;j++){
//                j2 = resolution*j;
//                for (int i=0;i<Nx;i++){
//                    i2 = resolution*i;
//                    int n = k*Nx*Ny+j*Nx+i;
//                    char VALUE=ReadID[n];
//                    for (int z = 0; z < resolution; z++) {
//                        for (int y = 0; y < resolution; y++) {
//                            for (int x = 0; x < resolution; x++) {
//                                int o = (k2+z)*(Nx2)*(Ny2) + (j2+y)*(Nx2) + (i2+x); // root lattice site of local block
//                                NewID.at(o) = VALUE;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        
//        // Step 2: Offset the NewID data into TmpID field
//        for (int k=0;k<Nz2;k++){
//            for (int j=0;j<Ny2;j++){
//                for (int i=0;i<Nx2;i++){
//                    int n = k*Nx2*Ny2+j*Nx2+i;
//                    int m = (k-(2*resolution-2)+resolution-1)*Nx2*Ny2 + (j-(2*resolution-2)+resolution-1)*Nx2+(i-(2*resolution-2)+resolution-1);
//                    if ( (k-(2*resolution-2)+resolution-1)*Nx2*Ny2      < 0 ) {}
//                    else if ( (j-(2*resolution-2)+resolution-1)*Nx2     < 0 ) {}
//                    else if  ( (i-(2*resolution-2)+resolution-1)        < 0 ) {}
//                    else { TmpID.at(m) = NewID.at(n); }
//                }
//            }
//        }
//
//        // Causes memory issues
//        for (int k=0;k<Nz2-(2*resolution-2);k++){
//            for (int j=0;j<Ny2-(2*resolution-2);j++){
//                for (int i=0;i<Nx2-(2*resolution-2);i++){
//                    int n = k*Nx2*Ny2+j*Nx2+i;
//                    char VALUE=TmpID.at(n);
//                    int m = k*Nx3*Ny3+j*Nx3+i;
//                    FinID.at(m) = VALUE;
//                    WriteID[m] = FinID.at(m);
//                }
//            }
//        }
//        
//        
//
//        
//        
//
//        int pcount = 0; 
//	int scount = 0;
//        sum_local = 0.0;
//        sum_sat_local = 0.0;
//        sum = 0;
//        sum_sat = 0.0;
//        
//        for (int k=1;k<Nz3-1;k++){
//            for (int j=1;j<Ny3-1;j++){
//                for (int i=1;i<Nx3-1;i++){
//                    int n = k*Nx3*Ny3+j*Nx3+i;
//                    if (FinID.at(n) == 1 || FinID.at(n) == 2){ pcount++; }
//                    if (FinID.at(n) == 1) { scount++; }
//                   // WriteID[n] = FinID.at(n);
//                }
//            }
//        }
//        
//        iVol_global = 1.0/(1.0*(Nx3-2)*(Ny3-2)*(Nz3-2)*nprocs);
//        
//	sum_local = 1.0*pcount;
//	sum_sat_local = 1.0*scount;
//
//        MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
//        porosity = porosity*iVol_global;
//
//
//	MPI_Barrier(comm);
//
//	MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
//
//	MPI_Barrier(comm);
//
//        if (Dm->rank()==0) printf("%dx Media porosity (internal)=%-15.15E ",resolution,porosity);
//        MPI_Allreduce(&sum_sat_local,&sum_sat,1,MPI_DOUBLE,MPI_SUM,comm);
//        saturation = 1.-sum_sat/sum;
//        if (Dm->rank()==0) printf(" WP saturation (internal)=%-15.15E\n",saturation);
//        
//	MPI_Barrier(comm);
//
//
//        sum_local = 0.0;
//        sum_sat_local = 0.0;
//        sum = 0;
//        sum_sat = 0.0;
//        
//         const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
//        
//        for (int k=1;k<Nz3-1;k++){
//            for (int j=1;j<Ny3-1;j++){
//                for (int i=1;i<Nx3-1;i++){
//                    for (int p=0; p<8; p++) {
//                    int n = i+cube[p][0] + (j+cube[p][1])*Nx3 + (k+cube[p][2])*Nx3*Ny3;
//                        if (FinID.at(n) == 1 || FinID.at(n) == 2){ sum_local+=0.125; }
//                        if (FinID.at(n) == 1) { sum_sat_local+=0.125; }
//                        
//
//                    }
//                }
//            }
//        }
//        
//    double porosity_halo, saturation_halo, porosity_error, saturation_error;
//        
//       
//        
//        
//        MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
//        porosity_halo = sum*iVol_global;
//        if (Dm->rank()==0) printf("%dx Media porosity (w/  halo)=%-15.15E ",resolution,porosity_halo);
//        
//	MPI_Barrier(comm);
//
//
//	MPI_Allreduce(&sum_sat_local,&sum_sat,1,MPI_DOUBLE,MPI_SUM,comm);
//        saturation_halo = 1.-sum_sat/sum;
//        if (Dm->rank()==0) printf(" WP saturation (w/  halo)=%-15.15E\n",saturation_halo);
//        if (Dm->rank()==0) { if(resolution == 1) printf("-------------------------------------------\n\n"); }
//        
//    
//    
//    
//    
//    porosity_error = abs(porosity - porosity_halo);
//    saturation_error = abs(saturation - saturation_halo);
//    if (Dm->rank()==0) printf("porosity error = %-15.15E\nsaturation error = %-15.15E\n\n",porosity_error,saturation_error);
//    
//    bool pass = true;
//    if (porosity_error > 0) {
//        if (Dm->rank()==0) printf("Error exists in porosity calculations, exiting...\n");
//        pass = false;
//    }
//    
//    if (saturation_error > 0) {
//        if (Dm->rank()==0) printf("Error exists in saturation calculations, exiting...\n");
//        pass = false;
//    }
//    
//        
////        }
////        }
//        
//        
//        //        for (int k=1;k<Nz3-1;k++){
//        //            for (int j=1;j<Ny3-1;j++){
//        //                for (int i=1;i<Nx3-1;i++){
//        //                    int n = k*Nx3*Ny3+j*Nx3+i;
//        //                    if (FinID[n] > 0){ sum_local+=1.0; }
//        //                    if (FinID[n] == 1) { sum_sat_local+=1.0; }
//        //                }
//        //            }
//        //        }
//        //               iVol_global = 1.0/(1.0*double((Nx3-2)*(Ny3-2)*(Nz3-2)*nprocs));
//        
//        
//        
//        
//        
//    
////        NewID = NULL;
////        TmpID = NULL;
////        FinID = NULL;
////        MPI_Barrier(comm);
//   // }
//    
//    
//    char LocalRankFilename[40];
//    sprintf(LocalRankFilename,"%s%s","newID.",LocalRankString);
//
//    // Write new ID files out
//    sprintf(LocalRankFilename,"%s%s","newID.",LocalRankString);
//    FILE *IDFILE = fopen(LocalRankFilename,"wb");
//    if (IDFILE==NULL) ERROR("Error opening file: newID.xxxxx");
//    fwrite(WriteID,1,WriteN,IDFILE);
//    fclose(IDFILE);
//
//    int numOfFiles = nprocx*nprocy*nprocz;
//    int subdomain = Nz-2;
//    int newsubdomain = subdomain*resolution;
//
//    if (rank == 0) {
//    printf("input subdomain size (z dir) =%d\n",subdomain);
//    printf("%d newID files written with resolution=%d\n",numOfFiles,resolution);
//    printf("output subdomain size (zidr) =%d\n\n",newsubdomain);
//    printf("---------------------------\nPASSED GRID RESOLUTION TEST\n---------------------------\n");
//    
//    }
//   delete[] ReadID;
//   delete[] WriteID;
//   
//    return pass;
//}
//
//
//bool ScaLBL_ColorModel::ReadInputAndValidate(){
//    size_t readID;
//    bool pass = true;
//    pass = Mask->ReadIDsAndValidationIDs();
//    for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read
//    
//    sprintf(LocalRankString,"%05d",rank);
//    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
//    sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
//    
//    // Generate the signed distance map
//    // Initialize the domain and communication
//    Array<char> id_solid(Nx,Ny,Nz);
//    int count = 0;
//    // Solve for the position of the solid phase
//    for (int k=0;k<Nz;k++){
//        for (int j=0;j<Ny;j++){
//            for (int i=0;i<Nx;i++){
//                int n = k*Nx*Ny+j*Nx+i;
//                // Initialize the solid phase
//                if (Mask->id[n] > 0)    id_solid(i,j,k) = 1;
//                else                       id_solid(i,j,k) = 0;
//            }
//        }
//    }
//    // Initialize the signed distance function
//    for (int k=0;k<Nz;k++){
//        for (int j=0;j<Ny;j++){
//            for (int i=0;i<Nx;i++){
//                int n=k*Nx*Ny+j*Nx+i;
//                // Initialize distance to +/- 1
//                Averages->SDs(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
//            }
//        }
//    }
//    //    MeanFilter(Averages->SDs);
//    if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
//    CalcDist(Averages->SDs,id_solid,*Mask);
//    
//    if (rank == 0) cout << "Domain set." << endl;
//    
//    //    Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha);
//   delete[] id;
//    return pass;
//}
//
//
//
//
//
//void ScaLBL_ColorModel::ReadInput(){
//    size_t readID;
//  
//    Mask->ReadIDs();
//    for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read
//
// 
//    sprintf(LocalRankString,"%05d",rank);
//    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
//    sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
//    
//    // Generate the signed distance map
//    // Initialize the domain and communication
//    Array<char> id_solid(Nx,Ny,Nz);
//    int count = 0;
//    // Solve for the position of the solid phase
//    for (int k=0;k<Nz;k++){
//        for (int j=0;j<Ny;j++){
//            for (int i=0;i<Nx;i++){
//                int n = k*Nx*Ny+j*Nx+i;
//                // Initialize the solid phase
//                if (Mask->id[n] > 0)    id_solid(i,j,k) = 1;
//                else                    id_solid(i,j,k) = 0;
//            }
//        }
//    }
//    
//    // Original Method
////    for (int k=0;k<Nz;k++){
////        for (int j=0;j<Ny;j++){
////            for (int i=0;i<Nx;i++){
////                int n=k*Nx*Ny+j*Nx+i;
////                // Initialize distance to +/- 1
////                Averages->SDs(i,j,k) = 2.0*double(id_solid(i,j,k)) -1.0;
////            }
////        }
////    }
////    CalcDist(Averages->SDs,id_solid,*Mask,{true,true,true},{1,1,1});
//
//    
//    
//    
//    // METHOD SORT OF WORKS
////        Array<double> TempDouble(Nx,Ny,Nz);
////        Array<double> TempDouble2(Nx,Ny,Nz);
//      for (int k=0;k<Nz;k++){
//          for (int j=0;j<Ny;j++){
//              for (int i=0;i<Nx;i++){
//                  int n=k*Nx*Ny+j*Nx+i;
//                  // Initialize distance to +/- 1
//                 Averages->SDs(i,j,k) = 2.0*double(id_solid(i,j,k)) -1;
//              }
//          }
//      }
//
//   CalcDist(Averages->SDs,id_solid,*Mask,{true,true,true},{1,1,1});
//    
//    
//  //  Dm->CommunicateMeshHalo(Averages->SDs);
//    
//    /*
//    MeanFilter(Averages->SDs,TempDouble);
//    Dm->CommunicateMeshHalo(TempDouble);
//    MeanFilter(TempDouble,Averages->SDs);
//    Dm->CommunicateMeshHalo(Averages->SDs);
//    MeanFilter(Averages->SDs,TempDouble);
//    Dm->CommunicateMeshHalo(TempDouble);
//    MeanFilter(TempDouble,Averages->SDs);
//    Dm->CommunicateMeshHalo(Averages->SDs);
//     */
//     
//
//
//    
//}
//
//
//
//
//
//
//
//
//
//
//
//void ScaLBL_ColorModel::AssignComponentLabels(double *phase)
//{
//    size_t NLABELS=0;
//    char VALUE=0;
//    double AFFINITY=0.f;
//    
//    auto LabelList = color_db->getVector<char>( "ComponentLabels" );
//    auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );
//    
//    NLABELS=LabelList.size();
//    if (NLABELS != AffinityList.size()){
//        ERROR("Error: ComponentLabels and ComponentAffinity must be the same length! \n");
//    }
//    
//    int label_count[NLABELS];
//    int label_count_global[NLABELS];
//    // Assign the labels
//    
//    for (int idx=0; idx<NLABELS; idx++) label_count[idx]=0;
//    
//    for (int k=0;k<Nz;k++){
//        for (int j=0;j<Ny;j++){
//            for (int i=0;i<Nx;i++){
//                int n = k*Nx*Ny+j*Nx+i;
//                VALUE=id[n];
//                // Assign the affinity from the paired list
//                for (unsigned int idx=0; idx < NLABELS; idx++){
//                    //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
//                    if (VALUE == LabelList[idx]){
//                        AFFINITY=AffinityList[idx];
//                        label_count[idx]++;
//                        idx = NLABELS;
//                        Mask->id[n] = 0; // set mask to zero since this is an immobile component
//                    }
//                }
//                // fluid labels are reserved
//                if (VALUE == 1) AFFINITY=1.0;
//                else if (VALUE == 2) AFFINITY=-1.0;
//                phase[n] = AFFINITY;
//            }
//        }
//    }
//    // Set Dm to match Mask
//    for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
//    
//    MPI_Allreduce(&label_count[0],&label_count_global[0],NLABELS,MPI_INT,MPI_SUM,Dm->Comm);
//    
//    if (rank==0){
//        printf("Components labels: %lu \n",NLABELS);
//        for (unsigned int idx=0; idx<NLABELS; idx++){
//            VALUE=LabelList[idx];
//            AFFINITY=AffinityList[idx];
//            double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocz);
//            printf("   label=%i, affinity=%f, volume fraction==%f\n",int(VALUE),AFFINITY,volume_fraction);
//        }
//    }
//    
//}
//
//
//void ScaLBL_ColorModel::Create(){
///*
//	 *  This function creates the variables needed to run a LBM 
//	 */
//	//.........................................................
//	// don't perform computations at the eight corners
//	//id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
//	//id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
//
//	//.........................................................
//	// Initialize communication structures in averaging domain
//	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
//	Mask->CommInit();
//	Np=Mask->PoreCount();
//	//...........................................................................
//	if (rank==0)    printf ("Create ScaLBL_Communicator \n");
//	// Create a communicator for the device (will use optimized layout)
//	// ScaLBL_Communicator ScaLBL_Comm(Mask); // original
//	ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
//	ScaLBL_Comm_Regular  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
//
//	int Npad=(Np/16 + 2)*16;
//	if (rank==0)    printf ("Set up memory efficient layout, %i | %i | %i \n", Np, Npad, N);
//	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
//	auto neighborList= new int[18*Npad];
//	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id,Np);
//	MPI_Barrier(comm);
//   
//	//...........................................................................
//	//                MAIN  VARIABLES ALLOCATED HERE
//	//...........................................................................
//	// LBM variables
//	if (rank==0)    printf ("Allocating distributions \n");
//	//......................device distributions.................................
//	dist_mem_size = Np*sizeof(double);
//	neighborSize=18*(Np*sizeof(int));
//	//...........................................................................
//	ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
//	ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
//	ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
//	ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*dist_mem_size);
//	ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*dist_mem_size);
//	ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
//	ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nx*Ny*Nz);		
//	ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
//	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
//	ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*sizeof(double)*Np);
//    ScaLBL_AllocateDeviceMemory((void **) &force, 19*sizeof(double)*Np);
//    ScaLBL_AllocateDeviceMemory((void **) &Xi, sizeof(int)*Np);
//    ScaLBL_AllocateDeviceMemory((void **) &TimeAveragedVelocity, sizeof(double)*Np);
//	//...........................................................................
//    
//     cout << "Generating a solid sphere ?" << endl;
//	// Update GPU data structures
//    
//    for (int i = 0; i < Nx*Ny*Nz; i++) Xi[i] = 0;
//    
//    double rad = double(Nx-2)*0.4;
//    double center = double(Nx-2)*0.5;
//    
//    for (int k=1; k<Nz-1; k++){
//        for (int j=1; j<Ny-1; j++){
//            for (int i=1; i<Nx-1; i++){
//                size_t n = k*Nx*Ny+j*Nx+i;
//                if ( (i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center) < rad*rad) Xi[n] = 1;
//            }
//        }
//    }
//
//        double *TAVel;
//        TAVel = new double[Np];
//    
//
//    for (int n = 0; n < Np; n++) TAVel[n] = n;
//    ScaLBL_CopyToDevice(TimeAveragedVelocity, TAVel, Np*sizeof(double));
//    delete[] TAVel;
//
//	if (rank==0)	printf ("Setting up device map and neighbor list \n");
//	fflush(stdout);
//	int *TmpMap;
//	TmpMap=new int[Np];
//	for (int k=1; k<Nz-1; k++){
//		for (int j=1; j<Ny-1; j++){
//			for (int i=1; i<Nx-1; i++){
//				int idx=Map(i,j,k);
//				if (!(idx < 0))
//					TmpMap[idx] = k*Nx*Ny+j*Nx+i;
//			}
//		}
//	}
//	// check that TmpMap is valid
//	for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
//		int n = TmpMap[idx];
//		if (n > Nx*Ny*Nz){
//			printf("Bad value! idx=%i \n");
//			TmpMap[idx] = Nx*Ny*Nz-1;
//		}
//	}
//	for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
//		int n = TmpMap[idx];
//		if (n > Nx*Ny*Nz){
//			printf("Bad value! idx=%i \n");
//			TmpMap[idx] = Nx*Ny*Nz-1;
//		}
//	}
//	ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
//	ScaLBL_DeviceBarrier();
//	delete [] TmpMap;
//	
//	// copy the neighbor list 
//	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
//	// initialize phi based on PhaseLabel (include solid component labels)
//	double *PhaseLabel;
//	PhaseLabel = new double[N];
//	AssignComponentLabels(PhaseLabel);
//	ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
//
//}        
//
///********************************************************
// * AssignComponentLabels                                 *
// ********************************************************/
//
//void ScaLBL_ColorModel::Initialize(){
//    
//    
//    if (rank==0)    printf ("Initializing distributions \n");
//    ScaLBL_D3Q19_Init(fq, Np);
//    
//    if (rank==0)    printf ("Initializing phase field \n");
////    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
////    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//    
//  
//    
//    /* Overwrite things if restart is true */
//    if (Restart == true){
//        
//        auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );
//        
//        if (rank==0){
//            printf("Reading restart file! \n");
//            ifstream restart("Restart.txt");
//            if (restart.is_open()){
//                restart  >> timestep;
//                printf("Restarting from timestep %i \n",timestep);
//            }
//            else{
//                printf("WARNING:No Restart.txt file, setting timestep=0 \n");
//                timestep=0;
//            }
//        }
//      
//        MPI_Bcast(&timestep,1,MPI_INT,0,comm);
//        // Read in the restart file to CPU buffers
//        double *cPhi, *cfq, *cDen, *cAq, *cBq;
//        cPhi = new double[N];
//        cDen = new double[2*Np];
//        cAq = new double[7*Np];
//        cBq = new double[7*Np];
//        cfq = new double[19*Np];
//        
//        
//        ifstream File(LocalRestartFile,ios::binary);
//        double value, valueAq, valueBq, va,vb,vphase;
//        
//        for (int n=0; n < Np; n++){
//            for (int q=0; q<7; q++) {
//                File.read((char*) &valueAq, sizeof(value));
//                File.read((char*) &valueBq, sizeof(value));
//                cAq[q*Np+n] = valueAq;
//                cBq[q*Np+n] = valueBq;
//            }
//        }
//        
//        for (int n=0; n<Np; n++){
//            // Read the distributions
//            for (int q=0; q<19; q++){
//                File.read((char*) &value, sizeof(value));
//                cfq[q*Np+n] = value;
//            }
//        }
//        
//        for (int n=0; n < Np; n++){
//            File.read((char*) &va, sizeof(va));
//            File.read((char*) &vb, sizeof(vb));
//            cDen[n]    = va;
//            cDen[Np+n] = vb;
//        }
//        
//        for (int n=0; n<N; n++){
//            File.read((char*) &vphase, sizeof(value));
//            cPhi[n] = vphase;
//            if (id[n] == 0) cPhi[n] = AffinityList[0];
//        }
//        
//        File.close();
//        
//        
//        // Copy the restart data to the GPU
//        ScaLBL_CopyToDevice(Den,cDen,2*Np*sizeof(double));
//        ScaLBL_CopyToDevice(fq,cfq,19*Np*sizeof(double));
//        ScaLBL_CopyToDevice(Phi,cPhi,N*sizeof(double));
//        ScaLBL_CopyToDevice(Aq,cAq,7*Np*sizeof(double));
//        ScaLBL_CopyToDevice(Bq,cBq,7*Np*sizeof(double));
//        ScaLBL_DeviceBarrier();
//        
//        MPI_Barrier(comm);
//      delete[] cDen;
//      delete[] cfq;
//      delete[] cPhi;
//      delete[]  cAq;
//      delete[] cBq;
//    }
//    
//    if (BoundaryCondition >0 ){
//        if (Dm->kproc()==0){
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
//        }
//        /* Change from -1 to 1 depending on if single or two phase*/
//        if (Dm->kproc() == nprocz-1){
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,Nz-1);
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,Nz-2);
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,Nz-3);
//        }
//    }
//    
//}
//
//std::string getLastLine(std::ifstream& in)
//{
//    std::string line;
//    while (in >> std::ws && std::getline(in, line)) // skip empty lines
//        ;
//
//    return line;
//}
//
//void ScaLBL_ColorModel::Run(){
//    int nprocs=nprocx*nprocy*nprocz;
//    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
//    
//    
//    int n_smooth = 2000;
//    double deltaX_star = 1./64.;
//    double beta_star = 0.35;
//    double exponent = 0.4;
//    double lambda = 1.0;
//    double deltaX = 1./(static_cast<double>(Nx*nprocx));
//    
//    
//   // int modtimestep = 0;
//    
//    /* Initial smoothing method */
////    if (Restart == false) {
////        timestep = -1.0 * n_smooth;
////    }
//    
//    double deltaT = deltaX*deltaX;
//    if (rank == 0) printf("deltaT value = %f\n",deltaT);
//    
//  //  double numerical_time = -static_cast<double>(n_smooth) * deltaX*deltaX;
//    
//    
//    int analysis_interval = analysis_db->getScalar<int>( "analysis_interval" );
//    
//    if (rank==0){
//        printf("********************************************************\n");
//        printf("No. of timesteps: %i \n", timestepMax);
//        fflush(stdout);
//    }
//    //.......create and start timer............
//    double starttime,stoptime,cputime;
//    ScaLBL_DeviceBarrier();
//    MPI_Barrier(comm);
//    starttime = MPI_Wtime();
//    //.........................................
//    double init_Fz = Fz;
//    
//    double png = 0;
//    double pwg = 0;
//    double diffp = 0;
//    double iftJwn = 0;
//
//    int halfTimestepMax = 100000;
//
//    double Ca_previous = 0.0;
//    
//    double Ca, muA, vA_z, relError;
//    Ca = muA = vA_z = 0.0;
//    relError = 1.0;
//    //************ MAIN ITERATION LOOP ***************************************/
//    PROFILE_START("Loop");
//    //std::shared_ptr<Database> analysis_db;
//    bool Regular = false;
//    runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, beta, Map );
////    timestep += analysis_interval-4;
////     analysis.run2( timestep, *Averages, Phi, Pressure, Velocity, fq, Aq, Bq, Den, Np);
//    
//   auto CaThresh = analysis_db->getScalar<size_t>( "CaThreshold" );
//    
//    double even_analysis_interval = static_cast<double>(analysis_interval)*0.5;
//    
//    auto count = 0;
//   
//	if (rank == 0) {
//    if (Restart == true) {
//     std::string line;
//    std::ifstream file("force.csv");
//
//    if (file)
//    {
//        line = getLastLine(file);
//        std::cout << "Restarting with force:" << line << '\n';
//    }
//    else{ std::cout << "Unable to open force file.\n"; }
//        
//    Fz = std::stod(line);
//    }
//   }
//    
//    double currentAffinity;
//    if (Restart == true) {
////        currentAffinity = -0.7;
////        if (rank == 0) {printf("Adjusting phase affinity, new affinity= %f \n",currentAffinity);}
////        currentAffinity = AdjustPhaseAffinity(dvcMap, Phi, N,currentAffinity,Dm->id);
//       
//    }
//
// 
//    MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
//    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
//
//    double epsilon = 0;
//    
//    cout << " We in the right function?" << endl;
//
//    //analysis.createThreads( analysis_method, 4 );
////    while (timestep < timestepMax) {
////        /* initial smoothing phase */
//////        if (Restart == false) {
//////            if (numerical_time < 0) {
//////                 epsilon = 1.0-(-0.1+(1.0*std::abs(numerical_time))/(static_cast<double>(n_smooth)*deltaT) - ( (1.0*std::abs(numerical_time))/(static_cast<double>(n_smooth)*deltaT) -1.0  )*beta_star*std::pow((deltaX/deltaX_star),exponent));
//////            }
//////        }
////        //beta = epsilon;
////
////
////        // if (rank == 0) { if (timestep%analysis_interval + 4 == analysis_interval){ printf("numerical_time=%f  beta=%f\n",numerical_time,beta); } }
////
////        //if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
////        PROFILE_START("Update");
////        // *************ODD TIMESTEP*************
////        timestep++; //modtimestep++;
////        //numerical_time = numerical_time + deltaT;
////        // Compute the Phase indicator field
////        // Read for Aq, Bq happens in this routine (requires communication)
////        ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
////        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
////        ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
////        ScaLBL_DeviceBarrier();
////        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
////
////        // Perform the collision operation
////        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
////        if (BoundaryCondition > 0){
////            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
////            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
////        }
////        // Halo exchange for phase field
////        ScaLBL_Comm_Regular->SendHalo(Phi);
////
////        ScaLBL_D3Q19_AAodd_Color_VP(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
////                                 alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np,Xi,force);
////        ScaLBL_Comm_Regular->RecvHalo(Phi);
////        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
////        ScaLBL_DeviceBarrier();
////        // Set BCs
////        if (BoundaryCondition == 3){
////            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
////            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
////        }
////        if (BoundaryCondition == 4){
////            din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
////            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
////        }
////        ScaLBL_D3Q19_AAodd_Color_VP(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
////                                 alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np,Xi,force);
////        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
////
////        // *************EVEN TIMESTEP*************
////        timestep++;// modtimestep++;
////        //numerical_time = numerical_time + deltaT;
////        // Compute the Phase indicator field
////        ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
////        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
////        ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
////        ScaLBL_DeviceBarrier();
////        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
////
////        // Perform the collision operation
////        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
////        // Halo exchange for phase field
////        if (BoundaryCondition > 0){
////            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
////            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
////        }
////        ScaLBL_Comm_Regular->SendHalo(Phi);
////        ScaLBL_D3Q19_AAeven_Color_VP(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
////                                  alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np,Xi,force);
////        ScaLBL_Comm_Regular->RecvHalo(Phi);
////        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
////        ScaLBL_DeviceBarrier();
////        // Set boundary conditions
////        if (BoundaryCondition == 3){
////            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
////            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
////        }
////        else if (BoundaryCondition == 4){
////            din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
////            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
////        }
////        /* Includes accumulated velocity field */
////        ScaLBL_D3Q19_AAeven_Color_VP(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
////                                     alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np,Xi,force);
////        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
////        //************************************************************************
////
////        MPI_Barrier(comm);
////        PROFILE_STOP("Update");
////
////
////        if (timestep%analysis_interval + 4 == analysis_interval){
////            Divide_Accumulated_Momentum(TimeAveragedVelocity,even_analysis_interval,Np);
////           // count = 0;
////        }
////
////        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
////
////        // Run the analysis
////        analysis.run3( timestep, *Averages, Phi, Pressure, Velocity, fq, Aq, Bq, Den, Np, TimeAveragedVelocity,Fx,Fy,Fz);
////
////
//////        if (rank == 0 ) {contactAngle = Averages->efawns_global;
//////
//////            printf("contact angle = %.10f \n",contactAngle);
//////
//////            if (contactAngle == 0 || contactAngle > 1.0) {
//////                // CPU:
//////                // Set phase field to new value thru function
//////
//////
//////                // GPU:
//////                // Set phase field to new value thru kernal
//////            }
//////
//////        }
////
////       // deltaT += 2*deltaT;
////      //  if (rank == 0) printf("t=%f\n",deltaT);
////
////        if (timestep%analysis_interval + 4 == analysis_interval){  Clear_Accumulated_Momentum(TimeAveragedVelocity,Np);  }
////
////
////
////
/////* Ca number work */
////        if ( timestep > analysis_interval && (timestep%analysis_interval == analysis_interval-20)){
////            vA_z = Averages->vn_global(2);
////            muA = rhoA*(tauA-0.5)/3.f;
////
////            png = Averages->pn_global;
////            pwg = Averages->pw_global;
////            diffp = png-pwg;
////            iftJwn = 5.796*alpha*Averages->Jwn_global;
////            double maxmagvel = Averages->maxMagVel_global;
////
////
////
////			Ca = vA_z*muA/(5.796*alpha*1.0);
////			if (Restart == true) { Ca = 100; Ca_previous = 200; Restart = false; }
////			if (Ca == 0) Ca = 1;
////                relError = fabs((Ca - Ca_previous)/Ca);
////                    if (rank==0){
////                        printf("timestep=%d Ca = %-15.15E, (previous = %-15.15E, rel error = %-15.15E) \n",timestep,Ca,Ca_previous,relError);
////                        FILE * kr_log_file = fopen("relperm.csv","a+");
////                        fprintf(kr_log_file,"%i,%-15.15E,%-15.15E,%-15.15E,%-15.15E,%-15.15E,%-15.15E\n",timestep,Ca,relError,Fz,beta,maxmagvel);
////                        fclose(kr_log_file);
////                        FILE * pressdiff = fopen("pressdiff.csv","a+");
////                        fprintf(pressdiff,"%i,%-15.15E,%-15.15E,%-15.15E,%-15.15E\n",timestep,png,pwg,diffp,iftJwn);
////                        fclose(pressdiff);
////                        FILE * force_file = fopen("force.csv","a+");
////                        fprintf(force_file,"%-15.15E\n",Fz);
////                        fclose(force_file);
////
////                     }
////                     Ca_previous = Ca;
////                 }
////          analysis.finish();
////
/////*
////    if (timestep >= 80000 && timestep <= 130000) {
////       if (timestep %  4000 == 0 && Ca > 1E-5) {
////
////	if (Fz > 1E-8) { Fz *= 0.8; }
////           printf("Fz=%-15.15E\n",Fz);
////       }
////    }
////
////
////if (timestep >= 180000 && timestep <= 230000) {
////       if (timestep %  4000 == 0 && Ca > 1E-5) {
////
////        if (Fz > 1E-8) { Fz *= 0.8; }
////           printf("Fz=%-15.15E\n",Fz);
////       }
////    }
////
////
////if (timestep >= 280000 && timestep <= 330000) {
////       if (timestep %  4000 == 0 && Ca > 1E-5) {
////
////        if (Fz > 1E-8) { Fz *= 0.8; }
////           printf("Fz=%-15.15E\n",Fz);
////       }
////    }
////*/
////
////
////
////
//// ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
////}
////
////    analysis.finish();
//    PROFILE_STOP("Loop");
//    PROFILE_SAVE("lbpm_color_simulator",1);
//    //************************************************************************
//    ScaLBL_DeviceBarrier();
//    MPI_Barrier(comm);
//    stoptime = MPI_Wtime();
//    if (rank==0) printf("-------------------------------------------------------------------\n");
//    // Compute the walltime per timestep
//    cputime = (stoptime - starttime)/timestep;
//    // Performance obtained from each node
//    double MLUPS = double(Np)/cputime/1000000;
//    
//    if (rank==0) printf("********************************************************\n");
//    if (rank==0) printf("CPU time = %f \n", cputime);
//    if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
//    MLUPS *= nprocs;
//    if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
//    if (rank==0) printf("********************************************************\n");
//    
//    // ************************************************************************
//}
//
//double ScaLBL_ColorModel::MorphInit(const double beta, const double morph_delta){
//    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
//    
//    double vF = 0.f;
//    double vS = 0.f;
//    
//    DoubleArray phase(Nx,Ny,Nz);
//    IntArray phase_label(Nx,Ny,Nz);;
//    DoubleArray phase_distance(Nx,Ny,Nz);
//    Array<char> phase_id(Nx,Ny,Nz);
//    
//    // Basic algorithm to
//    // 1. Copy phase field to CPU
//    ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));
//    
//    double count,count_global,volume_initial,volume_final;
//    count = 0.f;
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
//            }
//        }
//    }
//    MPI_Allreduce(&count,&count_global,1,MPI_DOUBLE,MPI_SUM,comm);
//    volume_initial = count_global;
//    
//    // 2. Identify connected components of phase field -> phase_label
//    BlobIDstruct new_index;
//    ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,phase,Averages->SDs,vF,vS,phase_label,comm);
//    MPI_Barrier(comm);
//    // only operate on component "0"
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                int label = phase_label(i,j,k);
//                if (label == 0 )     phase_id(i,j,k) = 0;
//                else 		     phase_id(i,j,k) = 1;
//            }
//        }
//    }
//    // 3. Generate a distance map to the largest object -> phase_distance
//    CalcDist(phase_distance,phase_id,*Dm);
//    
//    double temp,value;
//    double factor=0.5/beta;
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                if (phase_distance(i,j,k) < 3.f ){
//                    value = phase(i,j,k);
//                    if (value > 1.f)   value=1.f;
//                    if (value < -1.f)  value=-1.f;
//                    // temp -- distance based on analytical form McClure, Prins et al, Comp. Phys. Comm.
//                    temp = -factor*log((1.0+value)/(1.0-value));
//                    /// use this approximation close to the object
//                    if (fabs(value) < 0.8 && Averages->SDs(i,j,k) > 1.f ){
//                        phase_distance(i,j,k) = temp;
//                    }
//                }
//            }
//        }
//    }
//    
//    // 4. Apply erosion / dilation operation to phase_distance
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                double walldist=Averages->SDs(i,j,k);
//                double wallweight = 1.f / (1+exp(-5.f*(walldist-1.f)));
//                phase_distance(i,j,k) -= wallweight*morph_delta;
//            }
//        }
//    }
//    
//    // 5. Update phase indicator field based on new distnace
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                int n = k*Nx*Ny + j*Nx + i;
//                double d = phase_distance(i,j,k);
//                if (Averages->SDs(i,j,k) > 0.f){
//                    // only update phase field in immediate proximity of largest component
//                    if (d < 3.f){
//                        phase(i,j,k) = (2.f*(exp(-2.f*beta*d))/(1.f+exp(-2.f*beta*d))-1.f);
//                    }
//                }
//            }
//        }
//    }
//    
//    count = 0.f;
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
//            }
//        }
//    }
//    MPI_Allreduce(&count,&count_global,1,MPI_DOUBLE,MPI_SUM,comm);
//    volume_final=count_global;
//    
//    double delta_volume = (volume_final-volume_initial);
//    if (rank == 0)  printf("MorphInit: change fluid volume fraction by %f \n", delta_volume/volume_initial);
//    
//    // 6. copy back to the device
//    //if (rank==0)  printf("MorphInit: copy data  back to device\n");
//    ScaLBL_CopyToDevice(Phi,phase.data(),N*sizeof(double));
//    
//    // 7. Re-initialize phase field and density
////    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
////    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//    if (BoundaryCondition >0 ){
//        if (Dm->kproc()==0){
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
//            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
//        }
//        if (Dm->kproc() == nprocz-1){
//            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
//            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
//            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
//        }
//    }
//    return delta_volume;
//}
//
//void ScaLBL_ColorModel::WriteDebug(){
//    // Copy back final phase indicator field and convert to regular layout
//    DoubleArray PhaseField(Nx,Ny,Nz);
//    //ScaLBL_Comm->RegularLayout(Map,Phi,PhaseField);
//    ScaLBL_CopyToHost(PhaseField.data(), Phi, sizeof(double)*N);
//    
//    FILE *OUTFILE;
//    sprintf(LocalRankFilename,"Phase.%05i.raw",rank);
//    OUTFILE = fopen(LocalRankFilename,"wb");
//    fwrite(PhaseField.data(),8,N,OUTFILE);
//    fclose(OUTFILE);
//}
