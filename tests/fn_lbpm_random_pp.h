/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>


//#include "common/pmmc.h"
#include "common/Domain.h"
#include "common/SpherePack.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

inline void PopulateHalo_Char(char *id,
                              Domain &Dm, int nx, int ny, int nz, int rank) {
    double GlobalNumber = 1.f;
    double LocalNumber=0.f;
    int Nx = nx;
    int Ny = ny;
    int Nz = nz;
    int i,j,k,n;
    double count,countGlobal,totalGlobal;
    count = 0.f;
    double maxdist=0.f;
    double maxdistGlobal;
    int nprocx=Dm.nprocx();
    int nprocy=Dm.nprocy();
    int nprocz=Dm.nprocz();
    
    
    // Communication buffers
    char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
    char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
    char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
    char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
    char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
    char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
    // send buffers
    sendID_x = new char [Dm.sendCount_x];
    sendID_y = new char [Dm.sendCount_y];
    sendID_z = new char [Dm.sendCount_z];
    sendID_X = new char [Dm.sendCount_X];
    sendID_Y = new char [Dm.sendCount_Y];
    sendID_Z = new char [Dm.sendCount_Z];
    sendID_xy = new char [Dm.sendCount_xy];
    sendID_yz = new char [Dm.sendCount_yz];
    sendID_xz = new char [Dm.sendCount_xz];
    sendID_Xy = new char [Dm.sendCount_Xy];
    sendID_Yz = new char [Dm.sendCount_Yz];
    sendID_xZ = new char [Dm.sendCount_xZ];
    sendID_xY = new char [Dm.sendCount_xY];
    sendID_yZ = new char [Dm.sendCount_yZ];
    sendID_Xz = new char [Dm.sendCount_Xz];
    sendID_XY = new char [Dm.sendCount_XY];
    sendID_YZ = new char [Dm.sendCount_YZ];
    sendID_XZ = new char [Dm.sendCount_XZ];
    //......................................................................................
    // recv buffers
    recvID_x = new char [Dm.recvCount_x];
    recvID_y = new char [Dm.recvCount_y];
    recvID_z = new char [Dm.recvCount_z];
    recvID_X = new char [Dm.recvCount_X];
    recvID_Y = new char [Dm.recvCount_Y];
    recvID_Z = new char [Dm.recvCount_Z];
    recvID_xy = new char [Dm.recvCount_xy];
    recvID_yz = new char [Dm.recvCount_yz];
    recvID_xz = new char [Dm.recvCount_xz];
    recvID_Xy = new char [Dm.recvCount_Xy];
    recvID_xZ = new char [Dm.recvCount_xZ];
    recvID_xY = new char [Dm.recvCount_xY];
    recvID_yZ = new char [Dm.recvCount_yZ];
    recvID_Yz = new char [Dm.recvCount_Yz];
    recvID_Xz = new char [Dm.recvCount_Xz];
    recvID_XY = new char [Dm.recvCount_XY];
    recvID_YZ = new char [Dm.recvCount_YZ];
    recvID_XZ = new char [Dm.recvCount_XZ];
    //......................................................................................
    int sendtag,recvtag;
    sendtag = recvtag = 7;
    
    // Pack and send the updated ID values
    PackID(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
    PackID(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
    PackID(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
    PackID(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
    PackID(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
    PackID(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
    PackID(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
    PackID(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
    PackID(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
    PackID(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
    PackID(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
    PackID(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
    PackID(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
    PackID(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
    PackID(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
    PackID(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
    PackID(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
    PackID(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
    //......................................................................................
    MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_CHAR,Dm.rank_x(),sendtag,
                 recvID_X,Dm.recvCount_X,MPI_CHAR,Dm.rank_X(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_CHAR,Dm.rank_X(),sendtag,
                 recvID_x,Dm.recvCount_x,MPI_CHAR,Dm.rank_x(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_CHAR,Dm.rank_y(),sendtag,
                 recvID_Y,Dm.recvCount_Y,MPI_CHAR,Dm.rank_Y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_CHAR,Dm.rank_Y(),sendtag,
                 recvID_y,Dm.recvCount_y,MPI_CHAR,Dm.rank_y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_CHAR,Dm.rank_z(),sendtag,
                 recvID_Z,Dm.recvCount_Z,MPI_CHAR,Dm.rank_Z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_CHAR,Dm.rank_Z(),sendtag,
                 recvID_z,Dm.recvCount_z,MPI_CHAR,Dm.rank_z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_CHAR,Dm.rank_xy(),sendtag,
                 recvID_XY,Dm.recvCount_XY,MPI_CHAR,Dm.rank_XY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_CHAR,Dm.rank_XY(),sendtag,
                 recvID_xy,Dm.recvCount_xy,MPI_CHAR,Dm.rank_xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_CHAR,Dm.rank_Xy(),sendtag,
                 recvID_xY,Dm.recvCount_xY,MPI_CHAR,Dm.rank_xY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_CHAR,Dm.rank_xY(),sendtag,
                 recvID_Xy,Dm.recvCount_Xy,MPI_CHAR,Dm.rank_Xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_CHAR,Dm.rank_xz(),sendtag,
                 recvID_XZ,Dm.recvCount_XZ,MPI_CHAR,Dm.rank_XZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_CHAR,Dm.rank_XZ(),sendtag,
                 recvID_xz,Dm.recvCount_xz,MPI_CHAR,Dm.rank_xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_CHAR,Dm.rank_Xz(),sendtag,
                 recvID_xZ,Dm.recvCount_xZ,MPI_CHAR,Dm.rank_xZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_CHAR,Dm.rank_xZ(),sendtag,
                 recvID_Xz,Dm.recvCount_Xz,MPI_CHAR,Dm.rank_Xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_CHAR,Dm.rank_yz(),sendtag,
                 recvID_YZ,Dm.recvCount_YZ,MPI_CHAR,Dm.rank_YZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_CHAR,Dm.rank_YZ(),sendtag,
                 recvID_yz,Dm.recvCount_yz,MPI_CHAR,Dm.rank_yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_CHAR,Dm.rank_Yz(),sendtag,
                 recvID_yZ,Dm.recvCount_yZ,MPI_CHAR,Dm.rank_yZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_CHAR,Dm.rank_yZ(),sendtag,
                 recvID_Yz,Dm.recvCount_Yz,MPI_CHAR,Dm.rank_Yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    //......................................................................................
    UnpackID(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
    UnpackID(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
    UnpackID(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
    UnpackID(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
    UnpackID(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
    UnpackID(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
    UnpackID(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
    UnpackID(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
    UnpackID(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
    UnpackID(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
    UnpackID(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
    UnpackID(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
    UnpackID(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
    UnpackID(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
    UnpackID(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
    UnpackID(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
    UnpackID(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
    UnpackID(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
    //......................................................................................
    for (int k=0;k<2;k++){
        for (int j=0;j<2;j++){
            for (int i=0;i<2;i++){
                int n = k*Nx*Ny+j*Nx+i;
                if (i == 1 && j == 1 && k == 1) {int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);    id[m] = id[n]; }
            }
        }
    }
}


//inline double calculate_porosity(char * id,
//                                 int Nx, int Ny, int Nz, Domain &Dm, int rank) {
//    int nprocx,nprocy,nprocz;
//
//    nprocx = Dm.nprocx();
//    nprocy = Dm.nprocy();
//    nprocz = Dm.nprocz();
//    size_t count = 0;
//    size_t totalGlobal = 0;
//    for (int k=1; k<Nz-1; k++){
//        for (int j=1; j<Ny-1; j++){
//            for (int i=1; i<Nx-1; i++){
//                size_t n = k*Nx*Ny+j*Nx+i;
//                if (id[n] == 2 || id[n] == 1)  count++;
//            }
//        }
//    }
//
//    // total Global is the number of nodes in the pore-space
//    MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,Dm.Comm);
//    double volume=double(nprocx*nprocy*nprocz)*double(Nx-2)*double(Ny-2)*double(Nz-2);
//    double porosity=double(totalGlobal)/volume;
//    // if (rank == 0) printf("porosity = %.12f\n",porosity);
//
//    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
//    double dcount = 0;
//    double dcountGlobal = 0;
//    for (int k=1; k<Nz-1; k++){
//        for (int j=1; j<Ny-1; j++){
//            for (int i=1; i<Nx-1; i++){
//                for (int p=0;p<8;p++){
//                    int n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
//                    if (id[n] == 2 || id[n] == 1)  {  dcount += 0.125; }
//                }
//            }
//        }
//    }
//    MPI_Allreduce(&dcount,&dcountGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
//   // porosity=double(dcountGlobal)/volume;
//    // if (rank == 0) printf("porosity = %.12f\n",porosity);
//
//    return porosity;
//}


void lbpm_random_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    
    string filename = argv[1];
    
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );
    auto color_db = db->getDatabase( "Color" );
    auto analysis_db = db->getDatabase( "Analysis" );
    
    
    bool Restart = color_db->getScalar<bool>( "Restart" );
    
    // if (Restart == 0) {
    
    if (rank == 0) { printf("\nUsing legacy pre-processor: LBPM_RANDOM_PP \n");}
    
    
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, Nx, Ny, Nz, nspheres;
    double Lx, Ly, Lz;
    int i,j,k,n;
    int BC=0;
    
    
    //  auto db = std::make_shared<Database>( filename );
    //  auto domain_db = db->getDatabase( "Domain" );
    
    auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;
    int iproc = Dm->iproc();
    int jproc = Dm->jproc();
    int kproc = Dm->kproc();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
    nspheres = domain_db->getScalar<int>( "nspheres");
    int InitialWetting = domain_db->getScalar<int>("InitialSaturation");
    double Saturation = domain_db->getScalar<double>("RandomSaturation");
    /* Fact check these print statements... */
    // if (rank==0){ printf("Initializing wetting phase saturation of %f \n",Saturation);
    // if (InitialWetting == 1) {printf("Initial connected phase labeled (1) \n");} else  {printf("Initial connected phase labeled (2) \n");} }
    if (InitialWetting == 1) Saturation=1.0-Saturation;
    int BoundaryCondition=1;
    
    Nx += 2;
    Nx = Ny = Nz;
    size_t N = Nx*Ny*Nz;
    
    /* Initialize domain class communication */
    for (n=0; n<N; n++) Dm->id[n]=1;  Dm->CommInit();  for (n=0; n<N; n++) Dm->id[n]=2; // Set wetting
    
    char LocalRankFilename[40];
    
    
    char *ReadID;
    ReadID = new char[N];
    
    // Read the signed distance from file
    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"rb");
    size_t ReadSignDist;
    ReadSignDist=fread(ReadID,1,N,DIST);
    if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading ID file (rank=%i)\n",rank);
    fclose(DIST);
    
    // double porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
    // if (rank == 0) cout << "Read ID file porosity=" << porosity << endl;
    
    int count = 0;
    double denom = 0;
    double sum_local = 0;
    
    
    /* Update the Dm->id */
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                n = k*Nx*Ny+j*Nx+i;
                if (!(ReadID[n] == 0))  { count++;}
                
            }
        }
    }
    sum_local = 1.0*double(count); MPI_Allreduce(&sum_local,&denom,1,MPI_DOUBLE,MPI_SUM,comm);
    
    
    
    
    
    
    int bin, binCount;
    ifstream Dist("BlobSize.in");
    Dist >> binCount;
    int *SizeX, *SizeY, *SizeZ;
    //    printf("Number of blob sizes: %i \n",binCount);
    SizeX = new int [binCount];
    SizeY = new int [binCount];
    SizeZ = new int [binCount];
    for (bin=0; bin<binCount; bin++){
        Dist >> SizeX[bin];
        Dist >> SizeY[bin];
        Dist >> SizeZ[bin];
        //  printf("Blob %i dimension: %i x %i x %i \n",bin, SizeX[bin], SizeY[bin], SizeZ[bin]);
    }
    Dist.close();
    
    // Generate the residual NWP
    if (rank==0) printf("Initializing with saturation (phase 1) = %f \n",Saturation);
    //    GenerateResidual(id,Nx,Ny,Nz,Saturation);
    
    
    
    
    
    double sum_sat_local = 0;
    double sat = 0;
    
    
    count = 0;
    int x,y,z;
    int sizeX,sizeY,sizeZ;
    int ii,jj,kk;
    int Number = 0;        // number of features
    while (sat < Saturation) {
        if (rank==0) {
            Number++;
            // Randomly generate a point in the domain
            x = (Nx-2)*nprocx*float(rand())/float(RAND_MAX);
            y = (Ny-2)*nprocy*float(rand())/float(RAND_MAX);
            z = (Nz-2)*nprocz*float(rand())/float(RAND_MAX);
            
            bin = int(floor(binCount*float(rand())/float(RAND_MAX)));
            sizeX = SizeX[bin];
            sizeY = SizeY[bin];
            sizeZ = SizeZ[bin];
        }
        MPI_Bcast(&x,1,MPI_INT,0,comm);
        MPI_Bcast(&y,1,MPI_INT,0,comm);
        MPI_Bcast(&z,1,MPI_INT,0,comm);
        MPI_Bcast(&sizeX,1,MPI_INT,0,comm);
        MPI_Bcast(&sizeY,1,MPI_INT,0,comm);
        MPI_Bcast(&sizeZ,1,MPI_INT,0,comm);
        
        //if (rank==0) printf("Broadcast block at %i,%i,%i \n",x,y,z);
        
        for (k=z;k<z+sizeZ;k++){
            for (j=y;j<y+sizeY;j++){
                for (i=x;i<x+sizeX;i++){
                    
                    // Identify nodes in the domain (periodic BC)
                    ii = i;
                    jj = j;
                    kk = k;
                    
                    if (ii>nprocx*(Nx-2)) ii-=nprocx*(Nx-2);
                    if (jj>nprocy*(Ny-2)) jj-=nprocy*(Ny-2);
                    if (kk>nprocz*(Nz-2)) kk-=nprocz*(Nz-2);
                    
                    // Check if this is in the subdomain
                    if (ii < (iproc+1)*(Nx-2)+1 && jj < (jproc+1)*(Ny-2)+1 && kk < (kproc+1)*(Nz-2)+1 &&
                        ii  > iproc*(Nx-2) && jj > jproc*(Ny-2) && kk > kproc*(Nz-2) ){
                        
                        // Map from global to local coordinates
                        ii -= iproc*(Nx-2);
                        jj -= jproc*(Ny-2);
                        kk -= kproc*(Nz-2);
                        
                        n = kk*Nx*Ny+jj*Nx+ii;
                        
                        if (ReadID[n] == 2) { ReadID[n] = 1; }
                        
                    }
                }
            }
        }
        
        
        count = 0;
        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    n=k*Nx*Ny+j*Nx+i;
                    if (ReadID[n] == 1){
                        count++;
                    }
                }
            }
        }
        
        sat = 0;
        sum_sat_local=1.0*double(count);
        MPI_Allreduce(&sum_sat_local,&sat,1,MPI_DOUBLE,MPI_SUM,comm);
        
        
        sat = sat/denom;
        //if (rank==0) printf("New count=%i\n",countGlobal);
        if (rank==0) printf("New saturation=%f\n",sat);
    }
    
    
//    for ( k=0;k<2;k++){
//        for ( j=0;j<2;j++){
//            for ( i=0;i<2;i++){
//                n = k*Nx*Ny+j*Nx+i;
//                if (i == 0 && j == 0 && k == 0) ReadID[n] = 0;
//                 if (i == 1 && j == 1 && k == 1) {
////                     int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);
////                     ReadID[m] = ReadID[n] = 2;
//                     
//                 }
//            }
//        }
//    }
//    
//    for ( k=Nz-2;k<Nz;k++){
//            for ( j=Ny-2;j<Ny;j++){
//                for ( i=Nx-2;i<Nx;i++){
//                    n = k*Nx*Ny+j*Nx+i;
//                    if (i == Nx-1 && j == Ny-1 && k == Nz-1) ReadID[n] = 1;
//    //                     int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);
//    //                     ReadID[m] = ReadID[n] = 2;
//
//                     }
//                }
//            }
        
    
    //  if (InitialWetting == 1)    FlipID(ReadID,Nx*Ny*Nz);
    
    
    PopulateHalo_Char(ReadID, *Dm, Nx, Ny, Nz, rank);
    
    
    //double porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
    // if (rank == 0) cout << "Writing new ID files with porosity=" << porosity << endl;
    
    
    //sprintf(LocalRankFilename,"%s%s/ID.%05i","/mnt/bb/",usr,rank);
    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *ID = fopen(LocalRankFilename,"wb");
    fwrite(ReadID,1,N,ID);
    fclose(ID);
    
    if (rank == 0 ) cout << "End of random saturation pp." << endl << endl;
    //    }
    //    else {
    //
    //    }
}



void lbpm_partial_saturation_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    
    string filename = argv[1];
    
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );
    auto color_db = db->getDatabase( "Color" );
    auto analysis_db = db->getDatabase( "Analysis" );
    
    
    bool Restart = color_db->getScalar<bool>( "Restart" );
    
    // if (Restart == 0) {
    
    if (rank == 0) { printf("\nUsing legacy pre-processor: LBPM_PARTIAL_SATURATION_PP \n");}
    
    
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, Nx, Ny, Nz, nspheres;
    double Lx, Ly, Lz;
    int i,j,k,n;
    int BC=0;
    
    
    //  auto db = std::make_shared<Database>( filename );
    //  auto domain_db = db->getDatabase( "Domain" );
    
    auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;
    int iproc = Dm->iproc();
    int jproc = Dm->jproc();
    int kproc = Dm->kproc();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
    nspheres = domain_db->getScalar<int>( "nspheres");
    int InitialWetting = domain_db->getScalar<int>("InitialSaturation");
    double Saturation = domain_db->getScalar<double>("RandomSaturation");
    /* Fact check these print statements... */
    // if (rank==0){ printf("Initializing wetting phase saturation of %f \n",Saturation);
    // if (InitialWetting == 1) {printf("Initial connected phase labeled (1) \n");} else  {printf("Initial connected phase labeled (2) \n");} }
    if (InitialWetting == 1) Saturation=1.0-Saturation;
    int BoundaryCondition=1;
    
    Nx += 2;
    Nx = Ny = Nz;
    int N = Nx*Ny*Nz;
    
    /* Initialize domain class communication */
    for (n=0; n<N; n++) Dm->id[n]=1;  Dm->CommInit();  for (n=0; n<N; n++) Dm->id[n]=2; // Set wetting
    
    char LocalRankFilename[40];
    
    
    char *ReadID;
    ReadID = new char[N];
    
    // Read the signed distance from file
    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"rb");
    size_t ReadSignDist;
    ReadSignDist=fread(ReadID,1,N,DIST);
    if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading ID file (rank=%i)\n",rank);
    fclose(DIST);
    
    //double porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
    // if (rank == 0) cout << "Read ID file porosity=" << porosity << endl;
    
    int count = 0;
    double denom = 0;
    double sum_local = 0;
    
    
    /* Update the Dm->id */
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                n = k*Nx*Ny+j*Nx+i;
                if (!(ReadID[n] == 0))  { count++;}
            }
        }
    }
    
    // Single sphere in center of domain

        
        int vol_count = 0;
        
  double sphere_radius = 8;
      double sphere_center_x = 10;
    double sphere_center_y = 10;
      double sphere_center_z = 10;

        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    size_t n = k*Nx*Ny+j*Nx+i;
                    if (ReadID[n] == 2) {
                        // Box of fluid (cube to sphere or used in parallel plates)
//                                            if (i >= 10 && i <= 40 && j >= 10 && j <= 40 && k >= 10 && k <= 40  ){
//                                                ReadID[n] = 1;
//                                                vol_count++;
//                                            }
                        //                    if (i >= 6 && i <= 27 && j >= 6 && j <= 27 && k >= 6 && k <= 27  ){
                        //                        ReadID[n] = 1;
                        //                        vol_count++;
                        //                    }

                        // Sphere of fluid
                        //                    if ((i-16)*(i-16) + (j - 16)*(j - 16) + (k - 16)*(k - 16) <= double(12*12)   ){
                        //                        ReadID[n] = 1;
                        //                    }

    //                                        if ( ((i + (Nx-2)*iproc)-16)*((i + (Nx-2)*iproc)-16) + ((j + (Ny-2)*jproc)-16)*((j + (Ny-2)*jproc)-16) + ((k + (Nz-2)*kproc)-16)*((k + (Nz-2)*kproc)-16) <= sphere_radius*sphere_radius ) { ReadID[n] = 1; }



                        if  ( ((i + (Nx-2)*iproc)-sphere_center_x)*((i + (Nx-2)*iproc)-sphere_center_x)
                            + ((j + (Ny-2)*jproc)-sphere_center_y)*((j + (Ny-2)*jproc)-sphere_center_y)
                            + ((k + (Nz-2)*kproc)-sphere_center_z)*((k + (Nz-2)*kproc)-sphere_center_z)
                            <= sphere_radius*sphere_radius ) { ReadID[n] = 1; }

                        // Cube of domain filled with fluid
//                        if ((k + (Nz-2)*kproc) <  2*int(double(nprocz*(Nz-2))/double(3)) &&
//                            (k + (Nz-2)*kproc) > int(double(nprocz*(Nz-2))/double(3))
//
//                            && (j + (Ny-2)*jproc) <  2*int(double(nprocy*(Ny-2))/double(3)) &&
//                            (j + (Ny-2)*jproc) > int(double(nprocy*(Ny-2))/double(3))
//
//                            && (i + (Nx-2)*iproc) <  2*int(double(nprocx*(Nx-2))/double(3)) &&
//                            (i + (Nx-2)*iproc) > int(double(nprocx*(Nx-2))/double(3))     ) {
//
//                            ReadID[n] = 1;
//                        }

                        // square tube of fluid (used in parallel plates)
                        //                    if ((i + (Nx-2)*iproc) > int(double(nprocx*(Nx-2))/double(4)+2) && (i + (Nx-2)*iproc) < 3*int(double(nprocx*(Nx-2))/double(4))   ) {
                        //                        ReadID[n] = 1;
                        //                    }
                    }
                }
            }
        }
        
//        sphere_radius = 12;
//        sphere_center_x = 25;
//        sphere_center_y = 25;
//        sphere_center_z = 15;  // 0 for the periodic sphere
//
//        vol_count = 0;
//
//        // other side of sphere
//            for (int k=1; k<Nz-1; k++){
//                for (int j=1; j<Ny-1; j++){
//                    for (int i=1; i<Nx-1; i++){
//                        size_t n = k*Nx*Ny+j*Nx+i;
//                        if (ReadID[n] == 2) {
//                            // Box of fluid (cube to sphere or used in parallel plates)
////                                                if (i >= 10 && i <= 40 && j >= 10 && j <= 40 && k >= 10 && k <= 40  ){
////                                                    ReadID[n] = 1;
////                                                    vol_count++;
////                                                }
//                            //                    if (i >= 6 && i <= 27 && j >= 6 && j <= 27 && k >= 6 && k <= 27  ){
//                            //                        ReadID[n] = 1;
//                            //                        vol_count++;
//                            //                    }
//
//                            // Sphere of fluid
//                            //                    if ((i-16)*(i-16) + (j - 16)*(j - 16) + (k - 16)*(k - 16) <= double(12*12)   ){
//                            //                        ReadID[n] = 1;
//                            //                    }
//
//        //                                        if ( ((i + (Nx-2)*iproc)-16)*((i + (Nx-2)*iproc)-16) + ((j + (Ny-2)*jproc)-16)*((j + (Ny-2)*jproc)-16) + ((k + (Nz-2)*kproc)-16)*((k + (Nz-2)*kproc)-16) <= sphere_radius*sphere_radius ) { ReadID[n] = 1; }
//
//
//
//
//
//                            // Sphere
//                            if ( ((i + (Nx-2)*iproc)-sphere_center_x)*((i + (Nx-2)*iproc)-sphere_center_x) + ((j + (Ny-2)*jproc)-sphere_center_y)*((j + (Ny-2)*jproc)-sphere_center_y) + ((k + (Nz-2)*kproc)-sphere_center_z)*((k + (Nz-2)*kproc)-sphere_center_z) <= sphere_radius*sphere_radius ) { ReadID[n] = 1; }
//
//
//
//
//
//
//
////                             Slab of domain filled with fluid
////    if ((k + (Nz-2)*kproc) < 3*int(double(nprocz*(Nz-2))/double(4)) && (k + (Nz-2)*kproc) > int(double(nprocz*(Nz-2))/double(4))     ) { ReadID[n] = 1;    }
//
//                            // square tube of fluid (used in parallel plates)
//                            //                    if ((i + (Nx-2)*iproc) > int(double(nprocx*(Nx-2))/double(4)+2) && (i + (Nx-2)*iproc) < 3*int(double(nprocx*(Nx-2))/double(4))   ) {
//                            //                        ReadID[n] = 1;
//                            //                    }
//                        }
//                    }
//                }
//            }
//
//
//
//    sphere_radius = 8;
//    sphere_center_x = 25;
//    sphere_center_y = 25;
//    sphere_center_z = 38;
//
//    vol_count = 0;
//
//
//    // one side of sphere
//    for (int k=1; k<Nz-1; k++){
//        for (int j=1; j<Ny-1; j++){
//            for (int i=1; i<Nx-1; i++){
//                size_t n = k*Nx*Ny+j*Nx+i;
//                if (ReadID[n] == 2) {
//                    // Box of fluid (cube to sphere or used in parallel plates)
//                    //                    if (i >= 12 && i <= 55 && j >= 12 && j <= 55 && k >= 12 && k <= 55  ){
//                    //                        ReadID[n] = 1;
//                    //                        vol_count++;
//                    //                    }
//                    //                    if (i >= 6 && i <= 27 && j >= 6 && j <= 27 && k >= 6 && k <= 27  ){
//                    //                        ReadID[n] = 1;
//                    //                        vol_count++;
//                    //                    }
//
//                    // Sphere of fluid
//                    //                    if ((i-16)*(i-16) + (j - 16)*(j - 16) + (k - 16)*(k - 16) <= double(12*12)   ){
//                    //                        ReadID[n] = 1;
//                    //                    }
//
////                                        if ( ((i + (Nx-2)*iproc)-16)*((i + (Nx-2)*iproc)-16) + ((j + (Ny-2)*jproc)-16)*((j + (Ny-2)*jproc)-16) + ((k + (Nz-2)*kproc)-16)*((k + (Nz-2)*kproc)-16) <= sphere_radius*sphere_radius ) { ReadID[n] = 1; }
//
//
//
//                    if ( ((i + (Nx-2)*iproc)-sphere_center_x)*((i + (Nx-2)*iproc)-sphere_center_x) + ((j + (Ny-2)*jproc)-sphere_center_y)*((j + (Ny-2)*jproc)-sphere_center_y) + ((k + (Nz-2)*kproc)-sphere_center_z)*((k + (Nz-2)*kproc)-sphere_center_z) <= sphere_radius*sphere_radius ) { ReadID[n] = 1; }
//
//                    // Half of domain filled with fluid
//                    //                    if ((k + (Nz-2)*kproc) < int(double(nprocz*(Nz-2))/double(2))     ) {
//                    //                        ReadID[n] = 1;
//                    //                    }
//
//                    // square tube of fluid (used in parallel plates)
//                    //                    if ((i + (Nx-2)*iproc) > int(double(nprocx*(Nx-2))/double(4)+2) && (i + (Nx-2)*iproc) < 3*int(double(nprocx*(Nx-2))/double(4))   ) {
//                    //                        ReadID[n] = 1;
//                    //                    }
//                }
//            }
//        }
//    }


    
    // std::cout << "solid volume=" << vol_count << std::endl;
    sum_local = 1.0*double(count); MPI_Allreduce(&sum_local,&denom,1,MPI_DOUBLE,MPI_SUM,comm);
    
    
    //  PopulateHalo_Char(ReadID, *Dm, Nx, Ny, Nz, rank);
    
    
    //    porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
    //    if (rank == 0) cout << "Writing new ID files with porosity=" << porosity << endl;
    
    
    //sprintf(LocalRankFilename,"%s%s/ID.%05i","/mnt/bb/",usr,rank);
    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *ID = fopen(LocalRankFilename,"wb");
    fwrite(ReadID,1,N,ID);
    fclose(ID);
    
    if (rank == 0 ) cout << "End of random saturation pp." << endl << endl;
    //    }
    //    else {
    //
    //    }
}

void lbpm_morphopen_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
    double Lx, Ly, Lz;
    int i,j,k,n;
    int BC=0;
    //  char fluidValue,solidValue;
    int MAXTIME=1000;
    int READ_FROM_BLOCK=0;
    
    char LocalRankString[8];
    char LocalRankFilename[40];
    
    string filename;
    double Rcrit_new, SW;
    if (argc > 1){
        filename=argv[1];
        Rcrit_new=0.f;
        //SW=strtod(argv[2],NULL);
    }
    else ERROR("No input database provided\n");
    // read the input database
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );
    
    // Read domain parameters
    auto L = domain_db->getVector<double>( "L" );
    auto size = domain_db->getVector<int>( "n" );
    auto nproc = domain_db->getVector<int>( "nproc" );
    auto ReadValues = domain_db->getVector<char>( "ReadValues" );
    auto WriteValues = domain_db->getVector<char>( "WriteValues" );
    SW = domain_db->getScalar<double>("RandomSaturation");
    
    if (rank == 0) { printf("\nUsing legacy pre-processor: LBPM_MORPHOPEN_PP \n");}
    
    if (rank==0)    printf("Target saturation %f \n",SW);
    
    nx = size[0];
    ny = size[1];
    nz = size[2];
    nprocx = nproc[0];
    nprocy = nproc[1];
    nprocz = nproc[2];
    
    int N = (nx+2)*(ny+2)*(nz+2);
    
    std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
    //        std::shared_ptr<Domain> Dm (new Domain(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC));
    for (n=0; n<N; n++) Dm->id[n]=1;
    Dm->CommInit();
    
   
     
    
    char *ReadID;
      ReadID = new char[N];
      
      // Read the signed distance from file
      sprintf(LocalRankFilename,"ID.%05i",rank);
      FILE *DIST = fopen(LocalRankFilename,"rb");
      size_t ReadSignDist;
      ReadSignDist=fread(ReadID,1,N,DIST);
      if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading ID file (rank=%i)\n",rank);
      fclose(DIST);
    
    
    
    
    nx+=2; ny+=2; nz+=2;
    // Generate the signed distance map
    // Initialize the domain and communication
    Array<char> id_solid(nx,ny,nz);
    DoubleArray SignDist(nx,ny,nz);
    
  
    sprintf(LocalRankFilename,"SignDist.%05i",rank);
    double *SDs; SDs = new double[N];

    FILE *SDist = fopen(LocalRankFilename,"rb");
    size_t readID=fread(SDs,8,N,SDist);
    if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading SDist \n");
    fclose(SDist);
    for (size_t i=0; i < N; i++) SignDist(i) = SDs[i];
    
    MPI_Barrier(comm);
    double count,countGlobal,totalGlobal;
    count = 0.f;
    double maxdist=-200.f;
    double maxdistGlobal;
    for (int k=1; k<nz-1; k++){
        for (int j=1; j<ny-1; j++){
            for (int i=1; i<nx-1; i++){
                n = k*nx*ny+j*nx+i;
                // extract distance for critical radius
                if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
            }
        }
    }
    for (int k=0; k<nz; k++){
        for (int j=0; j<ny; j++){
            for (int i=0; i<nx; i++){
                n = k*nx*ny+j*nx+i;
                if (SignDist(i,j,k) < 0.f){
                    // don't do anything
                }
                else{
                    // initially saturated with wetting phase
                    ReadID[n] = 2;
                    count+=1.0;
                }
                // don't let halo be the maximum dist
                if ( SignDist(i,j,k) > maxdist) SignDist(i,j,k) = maxdist;
            }
        }
    }
    MPI_Barrier(comm);
    // total Global is the number of nodes in the pore-space
    MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
    MPI_Allreduce(&maxdist,&maxdistGlobal,1,MPI_DOUBLE,MPI_MAX,comm);
    double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
    double porosity=totalGlobal/volume;
    if (rank==0) printf("Media Porosity: %f \n",porosity);
    if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);
    
    Dm->CommInit();
    int iproc = Dm->iproc();
    int jproc = Dm->jproc();
    int kproc = Dm->kproc();
    
    // Generate the NWP configuration
    //if (rank==0) printf("Initializing morphological distribution with critical radius %f \n", Rcrit);
    if (rank==0) printf("Performing morphological opening with target saturation %f \n", SW);
    //    GenerateResidual(id,nx,ny,nz,Saturation);
    
    // Communication buffers
    char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
    char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
    char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
    char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
    char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
    char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
    // send buffers
    sendID_x = new char [Dm->sendCount_x];
    sendID_y = new char [Dm->sendCount_y];
    sendID_z = new char [Dm->sendCount_z];
    sendID_X = new char [Dm->sendCount_X];
    sendID_Y = new char [Dm->sendCount_Y];
    sendID_Z = new char [Dm->sendCount_Z];
    sendID_xy = new char [Dm->sendCount_xy];
    sendID_yz = new char [Dm->sendCount_yz];
    sendID_xz = new char [Dm->sendCount_xz];
    sendID_Xy = new char [Dm->sendCount_Xy];
    sendID_Yz = new char [Dm->sendCount_Yz];
    sendID_xZ = new char [Dm->sendCount_xZ];
    sendID_xY = new char [Dm->sendCount_xY];
    sendID_yZ = new char [Dm->sendCount_yZ];
    sendID_Xz = new char [Dm->sendCount_Xz];
    sendID_XY = new char [Dm->sendCount_XY];
    sendID_YZ = new char [Dm->sendCount_YZ];
    sendID_XZ = new char [Dm->sendCount_XZ];
    //......................................................................................
    // recv buffers
    recvID_x = new char [Dm->recvCount_x];
    recvID_y = new char [Dm->recvCount_y];
    recvID_z = new char [Dm->recvCount_z];
    recvID_X = new char [Dm->recvCount_X];
    recvID_Y = new char [Dm->recvCount_Y];
    recvID_Z = new char [Dm->recvCount_Z];
    recvID_xy = new char [Dm->recvCount_xy];
    recvID_yz = new char [Dm->recvCount_yz];
    recvID_xz = new char [Dm->recvCount_xz];
    recvID_Xy = new char [Dm->recvCount_Xy];
    recvID_xZ = new char [Dm->recvCount_xZ];
    recvID_xY = new char [Dm->recvCount_xY];
    recvID_yZ = new char [Dm->recvCount_yZ];
    recvID_Yz = new char [Dm->recvCount_Yz];
    recvID_Xz = new char [Dm->recvCount_Xz];
    recvID_XY = new char [Dm->recvCount_XY];
    recvID_YZ = new char [Dm->recvCount_YZ];
    recvID_XZ = new char [Dm->recvCount_XZ];
    //......................................................................................
    int sendtag,recvtag;
    sendtag = recvtag = 7;
    
    int x,y,z;
    int ii,jj,kk;
    int Nx = nx;
    int Ny = ny;
    int Nz = nz;
    
    double sw_old=1.0;
    double sw_new=1.0;
    double sw_diff_old = 1.0;
    double sw_diff_new = 1.0;
    
    // Increase the critical radius until the target saturation is met
    double deltaR=0.05; // amount to change the radius in voxel units
    double Rcrit_old;
    
    double GlobalNumber = 1.f;
    int imin,jmin,kmin,imax,jmax,kmax;
    
    Rcrit_new = maxdistGlobal;
    //if (argc>2){
    //    Rcrit_new = strtod(argv[2],NULL);
    //    if (rank==0) printf("Max. distance =%f, Initial critical radius = %f \n",maxdistGlobal,Rcrit_new);
    //}
    while (sw_new > SW)
    {
        sw_diff_old = sw_diff_new;
        sw_old = sw_new;
        Rcrit_old = Rcrit_new;
        Rcrit_new -= deltaR*Rcrit_old;
        int Window=round(Rcrit_new);
        if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
        // and sw<Sw will be immediately broken
        double LocalNumber=0.f;
        for(k=0; k<Nz; k++){
            for(j=0; j<Ny; j++){
                for(i=0; i<Nx; i++){
                    n = k*nx*ny + j*nx+i;
                    if (SignDist(i,j,k) > Rcrit_new){
                        // loop over the window and update
                        imin=max(1,i-Window);
                        jmin=max(1,j-Window);
                        kmin=max(1,k-Window);
                        imax=min(Nx-1,i+Window);
                        jmax=min(Ny-1,j+Window);
                        kmax=min(Nz-1,k+Window);
                        for (kk=kmin; kk<kmax; kk++){
                            for (jj=jmin; jj<jmax; jj++){
                                for (ii=imin; ii<imax; ii++){
                                    int nn = kk*nx*ny+jj*nx+ii;
                                    double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
                                    if (ReadID[nn] == 2 && dsq <= Rcrit_new*Rcrit_new){
                                        LocalNumber+=1.0;
                                        ReadID[nn]=1;
                                    }
                                }
                            }
                        }
                        
                    }
                    // move on
                }
            }
        }
        // Pack and send the updated ID values
        PackID(Dm->sendList_x, Dm->sendCount_x ,sendID_x, ReadID);
        PackID(Dm->sendList_X, Dm->sendCount_X ,sendID_X, ReadID);
        PackID(Dm->sendList_y, Dm->sendCount_y ,sendID_y, ReadID);
        PackID(Dm->sendList_Y, Dm->sendCount_Y ,sendID_Y, ReadID);
        PackID(Dm->sendList_z, Dm->sendCount_z ,sendID_z, ReadID);
        PackID(Dm->sendList_Z, Dm->sendCount_Z ,sendID_Z, ReadID);
        PackID(Dm->sendList_xy, Dm->sendCount_xy ,sendID_xy, ReadID);
        PackID(Dm->sendList_Xy, Dm->sendCount_Xy ,sendID_Xy, ReadID);
        PackID(Dm->sendList_xY, Dm->sendCount_xY ,sendID_xY, ReadID);
        PackID(Dm->sendList_XY, Dm->sendCount_XY ,sendID_XY, ReadID);
        PackID(Dm->sendList_xz, Dm->sendCount_xz ,sendID_xz, ReadID);
        PackID(Dm->sendList_Xz, Dm->sendCount_Xz ,sendID_Xz, ReadID);
        PackID(Dm->sendList_xZ, Dm->sendCount_xZ ,sendID_xZ, ReadID);
        PackID(Dm->sendList_XZ, Dm->sendCount_XZ ,sendID_XZ, ReadID);
        PackID(Dm->sendList_yz, Dm->sendCount_yz ,sendID_yz, ReadID);
        PackID(Dm->sendList_Yz, Dm->sendCount_Yz ,sendID_Yz, ReadID);
        PackID(Dm->sendList_yZ, Dm->sendCount_yZ ,sendID_yZ, ReadID);
        PackID(Dm->sendList_YZ, Dm->sendCount_YZ ,sendID_YZ, ReadID);
        //......................................................................................
        MPI_Sendrecv(sendID_x,Dm->sendCount_x,MPI_CHAR,Dm->rank_x(),sendtag,
                     recvID_X,Dm->recvCount_X,MPI_CHAR,Dm->rank_X(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_X,Dm->sendCount_X,MPI_CHAR,Dm->rank_X(),sendtag,
                     recvID_x,Dm->recvCount_x,MPI_CHAR,Dm->rank_x(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_y,Dm->sendCount_y,MPI_CHAR,Dm->rank_y(),sendtag,
                     recvID_Y,Dm->recvCount_Y,MPI_CHAR,Dm->rank_Y(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Y,Dm->sendCount_Y,MPI_CHAR,Dm->rank_Y(),sendtag,
                     recvID_y,Dm->recvCount_y,MPI_CHAR,Dm->rank_y(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_z,Dm->sendCount_z,MPI_CHAR,Dm->rank_z(),sendtag,
                     recvID_Z,Dm->recvCount_Z,MPI_CHAR,Dm->rank_Z(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Z,Dm->sendCount_Z,MPI_CHAR,Dm->rank_Z(),sendtag,
                     recvID_z,Dm->recvCount_z,MPI_CHAR,Dm->rank_z(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xy,Dm->sendCount_xy,MPI_CHAR,Dm->rank_xy(),sendtag,
                     recvID_XY,Dm->recvCount_XY,MPI_CHAR,Dm->rank_XY(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_XY,Dm->sendCount_XY,MPI_CHAR,Dm->rank_XY(),sendtag,
                     recvID_xy,Dm->recvCount_xy,MPI_CHAR,Dm->rank_xy(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Xy,Dm->sendCount_Xy,MPI_CHAR,Dm->rank_Xy(),sendtag,
                     recvID_xY,Dm->recvCount_xY,MPI_CHAR,Dm->rank_xY(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xY,Dm->sendCount_xY,MPI_CHAR,Dm->rank_xY(),sendtag,
                     recvID_Xy,Dm->recvCount_Xy,MPI_CHAR,Dm->rank_Xy(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xz,Dm->sendCount_xz,MPI_CHAR,Dm->rank_xz(),sendtag,
                     recvID_XZ,Dm->recvCount_XZ,MPI_CHAR,Dm->rank_XZ(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_XZ,Dm->sendCount_XZ,MPI_CHAR,Dm->rank_XZ(),sendtag,
                     recvID_xz,Dm->recvCount_xz,MPI_CHAR,Dm->rank_xz(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Xz,Dm->sendCount_Xz,MPI_CHAR,Dm->rank_Xz(),sendtag,
                     recvID_xZ,Dm->recvCount_xZ,MPI_CHAR,Dm->rank_xZ(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xZ,Dm->sendCount_xZ,MPI_CHAR,Dm->rank_xZ(),sendtag,
                     recvID_Xz,Dm->recvCount_Xz,MPI_CHAR,Dm->rank_Xz(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_yz,Dm->sendCount_yz,MPI_CHAR,Dm->rank_yz(),sendtag,
                     recvID_YZ,Dm->recvCount_YZ,MPI_CHAR,Dm->rank_YZ(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_YZ,Dm->sendCount_YZ,MPI_CHAR,Dm->rank_YZ(),sendtag,
                     recvID_yz,Dm->recvCount_yz,MPI_CHAR,Dm->rank_yz(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Yz,Dm->sendCount_Yz,MPI_CHAR,Dm->rank_Yz(),sendtag,
                     recvID_yZ,Dm->recvCount_yZ,MPI_CHAR,Dm->rank_yZ(),recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_yZ,Dm->sendCount_yZ,MPI_CHAR,Dm->rank_yZ(),sendtag,
                     recvID_Yz,Dm->recvCount_Yz,MPI_CHAR,Dm->rank_Yz(),recvtag,comm,MPI_STATUS_IGNORE);
        //......................................................................................
        UnpackID(Dm->recvList_x, Dm->recvCount_x ,recvID_x, ReadID);
        UnpackID(Dm->recvList_X, Dm->recvCount_X ,recvID_X, ReadID);
        UnpackID(Dm->recvList_y, Dm->recvCount_y ,recvID_y, ReadID);
        UnpackID(Dm->recvList_Y, Dm->recvCount_Y ,recvID_Y, ReadID);
        UnpackID(Dm->recvList_z, Dm->recvCount_z ,recvID_z, ReadID);
        UnpackID(Dm->recvList_Z, Dm->recvCount_Z ,recvID_Z, ReadID);
        UnpackID(Dm->recvList_xy, Dm->recvCount_xy ,recvID_xy, ReadID);
        UnpackID(Dm->recvList_Xy, Dm->recvCount_Xy ,recvID_Xy, ReadID);
        UnpackID(Dm->recvList_xY, Dm->recvCount_xY ,recvID_xY, ReadID);
        UnpackID(Dm->recvList_XY, Dm->recvCount_XY ,recvID_XY, ReadID);
        UnpackID(Dm->recvList_xz, Dm->recvCount_xz ,recvID_xz, ReadID);
        UnpackID(Dm->recvList_Xz, Dm->recvCount_Xz ,recvID_Xz, ReadID);
        UnpackID(Dm->recvList_xZ, Dm->recvCount_xZ ,recvID_xZ, ReadID);
        UnpackID(Dm->recvList_XZ, Dm->recvCount_XZ ,recvID_XZ, ReadID);
        UnpackID(Dm->recvList_yz, Dm->recvCount_yz ,recvID_yz, ReadID);
        UnpackID(Dm->recvList_Yz, Dm->recvCount_Yz ,recvID_Yz, ReadID);
        UnpackID(Dm->recvList_yZ, Dm->recvCount_yZ ,recvID_yZ, ReadID);
        UnpackID(Dm->recvList_YZ, Dm->recvCount_YZ ,recvID_YZ, ReadID);
        //......................................................................................
        
        MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_DOUBLE,MPI_SUM,comm);
        
        count = 0.f;
        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    n=k*Nx*Ny+j*Nx+i;
                    if (ReadID[n] == 2){
                        count+=1.0;
                    }
                }
            }
        }
        MPI_Allreduce(&count,&countGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
        sw_new = countGlobal/totalGlobal;
        sw_diff_new = abs(sw_new-SW);
        if (rank==0){
            printf("     %f ",sw_new);
            printf("     %f\n",Rcrit_new);
        }
    }
    
    if (sw_diff_new<sw_diff_old){
        if (rank==0){
            printf("Final saturation=%f\n",sw_new);
            printf("Final critical radius=%f\n",Rcrit_new);
        }
    }
    else{
        if (rank==0){
            printf("Final saturation=%f\n",sw_old);
            printf("Final critical radius=%f\n",Rcrit_old);
        }
    }
    
    if (rank==0) printf("Writing ID file \n");
    sprintf(LocalRankFilename,"ID.%05i",rank);
    
    char *id;
    id = new char[N];
    
    FILE *ID = fopen(LocalRankFilename,"wb");
    fwrite(ReadID,1,N,ID);
    fclose(ID);
    
}


