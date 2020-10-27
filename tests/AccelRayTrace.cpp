#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

//#include "common/pmmc.h"
#include "common/Domain.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

#include "AccelRayTrace.h"

using namespace std;

int RT_DOUBLE = 8;
int RT_INT = 4;
int RT_CHAR = 1;

inline double calculate_porosity(char * id,
                              int Nx, int Ny, int Nz, Domain &Dm, int rank) {
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
                if (id[n] == 2)  count++;
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
                    if (id[n] == 2)  {  dcount += 0.125; }
                }
            }
        }
    }
    MPI_Allreduce(&dcount,&dcountGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    porosity=double(dcountGlobal)/volume;
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
 
 return porosity;
}


inline void fillBoolfield(bool* boolfield,
                          int Nx, int Ny, int Nz) {
    
    
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                size_t n = k*Nx*Ny+j*Nx+i;
                boolfield[n] = false;
            }
        }
    }
    
    
    
} // Sets interior to false, keeps halo true

inline void PackID(int *list,
                   int count, char *sendbuf, char *ID) {
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}

inline void UnpackID(int *list,
                     int count, char *recvbuf, char *ID) {
    // Fill in the phase ID values from neighboring processors
    // This unpacks the values once they have been recieved from neighbors
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

inline void PackID_Double(int *list,
                          int count, double *sendbuf, double *ID) {
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}

inline void UnpackID_Double(int *list,
                            int count, double *recvbuf, double *ID) {
    // Fill in the phase ID values from neighboring processors
    // This unpacks the values once they have been recieved from neighbors
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

inline void PopulateHalo_Double(double *id,
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
    double *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
    double *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
    double *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
    double *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
    double *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
    double *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
    // send buffers
    sendID_x = new double [Dm.sendCount_x];
    sendID_y = new double [Dm.sendCount_y];
    sendID_z = new double [Dm.sendCount_z];
    sendID_X = new double [Dm.sendCount_X];
    sendID_Y = new double [Dm.sendCount_Y];
    sendID_Z = new double [Dm.sendCount_Z];
    sendID_xy = new double [Dm.sendCount_xy];
    sendID_yz = new double [Dm.sendCount_yz];
    sendID_xz = new double [Dm.sendCount_xz];
    sendID_Xy = new double [Dm.sendCount_Xy];
    sendID_Yz = new double [Dm.sendCount_Yz];
    sendID_xZ = new double [Dm.sendCount_xZ];
    sendID_xY = new double [Dm.sendCount_xY];
    sendID_yZ = new double [Dm.sendCount_yZ];
    sendID_Xz = new double [Dm.sendCount_Xz];
    sendID_XY = new double [Dm.sendCount_XY];
    sendID_YZ = new double [Dm.sendCount_YZ];
    sendID_XZ = new double [Dm.sendCount_XZ];
    //......................................................................................
    // recv buffers
    recvID_x = new double [Dm.recvCount_x];
    recvID_y = new double [Dm.recvCount_y];
    recvID_z = new double [Dm.recvCount_z];
    recvID_X = new double [Dm.recvCount_X];
    recvID_Y = new double [Dm.recvCount_Y];
    recvID_Z = new double [Dm.recvCount_Z];
    recvID_xy = new double [Dm.recvCount_xy];
    recvID_yz = new double [Dm.recvCount_yz];
    recvID_xz = new double [Dm.recvCount_xz];
    recvID_Xy = new double [Dm.recvCount_Xy];
    recvID_xZ = new double [Dm.recvCount_xZ];
    recvID_xY = new double [Dm.recvCount_xY];
    recvID_yZ = new double [Dm.recvCount_yZ];
    recvID_Yz = new double [Dm.recvCount_Yz];
    recvID_Xz = new double [Dm.recvCount_Xz];
    recvID_XY = new double [Dm.recvCount_XY];
    recvID_YZ = new double [Dm.recvCount_YZ];
    recvID_XZ = new double [Dm.recvCount_XZ];
    //......................................................................................
    int sendtag,recvtag;
    sendtag = recvtag = 7;
    
    // Pack and send the updated ID values
    PackID_Double(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
    PackID_Double(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
    PackID_Double(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
    PackID_Double(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
    PackID_Double(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
    PackID_Double(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
    PackID_Double(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
    PackID_Double(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
    PackID_Double(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
    PackID_Double(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
    PackID_Double(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
    PackID_Double(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
    PackID_Double(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
    PackID_Double(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
    PackID_Double(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
    PackID_Double(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
    PackID_Double(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
    PackID_Double(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
    //......................................................................................
    MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_DOUBLE,Dm.rank_x(),sendtag,
                 recvID_X,Dm.recvCount_X,MPI_DOUBLE,Dm.rank_X(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_DOUBLE,Dm.rank_X(),sendtag,
                 recvID_x,Dm.recvCount_x,MPI_DOUBLE,Dm.rank_x(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_DOUBLE,Dm.rank_y(),sendtag,
                 recvID_Y,Dm.recvCount_Y,MPI_DOUBLE,Dm.rank_Y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_DOUBLE,Dm.rank_Y(),sendtag,
                 recvID_y,Dm.recvCount_y,MPI_DOUBLE,Dm.rank_y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_DOUBLE,Dm.rank_z(),sendtag,
                 recvID_Z,Dm.recvCount_Z,MPI_DOUBLE,Dm.rank_Z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_DOUBLE,Dm.rank_Z(),sendtag,
                 recvID_z,Dm.recvCount_z,MPI_DOUBLE,Dm.rank_z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_DOUBLE,Dm.rank_xy(),sendtag,
                 recvID_XY,Dm.recvCount_XY,MPI_DOUBLE,Dm.rank_XY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_DOUBLE,Dm.rank_XY(),sendtag,
                 recvID_xy,Dm.recvCount_xy,MPI_DOUBLE,Dm.rank_xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_DOUBLE,Dm.rank_Xy(),sendtag,
                 recvID_xY,Dm.recvCount_xY,MPI_DOUBLE,Dm.rank_xY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_DOUBLE,Dm.rank_xY(),sendtag,
                 recvID_Xy,Dm.recvCount_Xy,MPI_DOUBLE,Dm.rank_Xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_DOUBLE,Dm.rank_xz(),sendtag,
                 recvID_XZ,Dm.recvCount_XZ,MPI_DOUBLE,Dm.rank_XZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_DOUBLE,Dm.rank_XZ(),sendtag,
                 recvID_xz,Dm.recvCount_xz,MPI_DOUBLE,Dm.rank_xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_DOUBLE,Dm.rank_Xz(),sendtag,
                 recvID_xZ,Dm.recvCount_xZ,MPI_DOUBLE,Dm.rank_xZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_DOUBLE,Dm.rank_xZ(),sendtag,
                 recvID_Xz,Dm.recvCount_Xz,MPI_DOUBLE,Dm.rank_Xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_DOUBLE,Dm.rank_yz(),sendtag,
                 recvID_YZ,Dm.recvCount_YZ,MPI_DOUBLE,Dm.rank_YZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_DOUBLE,Dm.rank_YZ(),sendtag,
                 recvID_yz,Dm.recvCount_yz,MPI_DOUBLE,Dm.rank_yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_DOUBLE,Dm.rank_Yz(),sendtag,
                 recvID_yZ,Dm.recvCount_yZ,MPI_DOUBLE,Dm.rank_yZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_DOUBLE,Dm.rank_yZ(),sendtag,
                 recvID_Yz,Dm.recvCount_Yz,MPI_DOUBLE,Dm.rank_Yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    //......................................................................................
    UnpackID_Double(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
    UnpackID_Double(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
    UnpackID_Double(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
    UnpackID_Double(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
    UnpackID_Double(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
    UnpackID_Double(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
    UnpackID_Double(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
    UnpackID_Double(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
    UnpackID_Double(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
    UnpackID_Double(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
    UnpackID_Double(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
    UnpackID_Double(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
    UnpackID_Double(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
    UnpackID_Double(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
    UnpackID_Double(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
    UnpackID_Double(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
    UnpackID_Double(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
    UnpackID_Double(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
    //......................................................................................
    
    /* Correct corner distance map */
    for (int k=0;k<2;k++){
       for (int j=0;j<2;j++){
           for (int i=0;i<2;i++){
               int n = k*Nx*Ny+j*Nx+i;
               if (i == 1 && j == 1 && k == 1) {int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);    id[m] = id[n]; }
           }
       }
    }
}

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

int main(int argc, char **argv)
{
    // Initialize MPI
    int rank,nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    {
        if (rank == 0) cout << "========================================\nRAY TRACING PRE-PROCESSOR\n========================================\n";
        size_t m,n;
        size_t count = 0;
        
        int Nx,Ny,Nz;
        double Lx,Ly,Lz;
        int iproc,jproc,kproc;
        int nprocx,nprocy,nprocz;

        // Filenames used
        char LocalRankString[8];
        char LocalRankFilename1[40];
        char LocalRankFilename2[40];

        char LocalRankFilename3[40];
        char LocalRankFilename4[40];
        sprintf(LocalRankString,"%05d",rank);
        
        /* Read in input.db information */
        string filename;
        if (argc > 1) filename=argv[1];
        else ERROR("No input database provided\n");
        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));
        Nx = Dm->Nx;  Ny = Dm->Ny;  Nz = Dm->Nz;
        Lx = Dm->Lx;  Ly = Dm->Ly;  Lz = Dm->Lz;
        iproc = Dm->iproc();  jproc = Dm->jproc(); kproc = Dm->kproc();
        nprocx = Dm->nprocx();  nprocy = Dm->nprocy();  nprocz = Dm->nprocz();
        
        int nspheres = domain_db->getScalar<int>( "nspheres");
        // if (rank==0) printf("Applying ray tracing to );
        auto ReadValues = domain_db->getVector<char>( "ReadValues" );  // Should this be int??
        auto WriteValues = domain_db->getVector<char>( "WriteValues" );
        
        auto rt_db = db->getDatabase( "RayTrace" );
        int geometry = rt_db->getScalar<int>( "geometry" );
        double offset = rt_db->getScalar<double>( "offset" );
        
        if (rank == 0) { cout << "Geometry:";
            if (geometry == 1 ) cout << "periodic sphere pack";
            if (geometry == 2 ) cout << "slab" << endl << endl;
            if (geometry == 3 ) cout << "parallel plates";
            if (geometry == 4 ) cout << "square tube" << endl << endl;
            if (geometry == 5 ) cout << "capillary tube" << endl << endl;
            /* Parameters for particular methods */
            if (geometry == 1 ) cout << " with " << nspheres << " spheres" << endl << endl;
            if (geometry == 3 ) cout << " with offset=" << offset << endl << endl;
        }
        
        /* Set up domain sizes with halo */
        size_t Nq = 18;
        Nx += 2;
        Nx = Ny = Nz;
        size_t N = Nx*Ny*Nz;
       

        size_t NqN = Nq*N;

        /* Initialize domain class communication */
        for (n=0; n<N; n++) Dm->id[n] = 1;  Dm->CommInit();  for (n=0; n<N; n++) Dm->id[n] = 2; // Set wetting phase (id=2)

        std::cout << " After CommInit " << std::endl;



        /* Declare ray tracing data arrays */
        auto boolfield = new bool[N];  for (n=0; n<N; n++) boolfield[n]=true; fillBoolfield(boolfield, Nx, Ny, Nz);
        
        /* Initialize ray-trace data */
        char* LIBB_id;  LIBB_id = new char[N];  for (size_t n=0; n<N; n++) LIBB_id[n] = 2; // Set domain to wetting fluid
        double*  qList;   qList = new double[NqN]; for (size_t n = 0; n < NqN; n++) qList[n] = -1;
        double* rtDistance; rtDistance = new double[N];  for (size_t n = 0; n < N; n++) rtDistance[n] = 1000.0;
        double* SDs; SDs = new double[N];  for (size_t n = 0; n < N; n++) SDs[n] = 1000.0; 
        

        // Sphere pack
        if ( geometry == 1)
        {
            run_sp_raytrace(LIBB_id,
                            qList,
                            rtDistance,SDs,
                            Nx,Ny,Nz,
                            Lx,
                            nspheres,
                            iproc,jproc,kproc,
                            nprocx,nprocy,nprocz,
                            rank,
                            N,NqN);
            
//            ComputeLIBBqDistances_sp(LIBB_id,
//                                     qList, rtDistance, SDs,
//                                     Nx, Ny, Nz,
//                                     Lx,
//                                     nspheres,
//                                     iproc, jproc, kproc,
//                                     nprocx, nprocy, nprocz,
//                                     rank);
        }
        
        /* Clear halo and then populate halo */
//        for (int n = 0; n < N; n++) if (boolfield[n] == true) LIBB_id[n] = 0;
//        for (int n = 0; n < N; n++) if (boolfield[n] == true) rtDistance[n] = 1;
//        PopulateHalo_Char(LIBB_id, *Dm, Nx, Ny, Nz, rank);
//        PopulateHalo_Double(rtDistance, *Dm, Nx, Ny, Nz, rank);
        
  //      MPI_Barrier(comm);
        
        // ComputeHWBBqDistances(LIBB_id, qList, Nx, Ny, Nz, count);
        
         MPI_Barrier(comm); 
        
        
        
        /* Check the porosity */
        double porosity = calculate_porosity(LIBB_id, Nx, Ny, Nz, *Dm, rank);
        if (rank == 0) cout << "Actual id porosity=" << porosity << endl;
        



        MPI_Barrier(comm);

        sprintf(LocalRankFilename1,"%s%s","RTD.",LocalRankString);
        FILE *RTFILE = fopen(LocalRankFilename1,"wb");
        if (RTFILE==NULL) ERROR("Error opening file: rt.xxxxx");
        fwrite(rtDistance,8,N,RTFILE);
        fclose(RTFILE);

        
       MPI_Barrier(comm);       
        
 
        sprintf(LocalRankFilename2,"%s%s","ID.",LocalRankString);
        FILE *IDFILE = fopen(LocalRankFilename2,"wb");
        if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
        fwrite(LIBB_id,1,N,IDFILE);
        fclose(IDFILE);

        MPI_Barrier(comm);
        
        sprintf(LocalRankFilename3,"%s%s","q.",LocalRankString);
        FILE *QFILE = fopen(LocalRankFilename3,"wb");
        if (QFILE==NULL) ERROR("Error opening file: q.xxxxx");
        fwrite(qList,8,NqN,QFILE);
        fclose(QFILE);



        MPI_Barrier(comm);


        sprintf(LocalRankFilename4,"%s%s","SignDist.",LocalRankString);
        FILE *SDFILE = fopen(LocalRankFilename4,"wb");
        if (SDFILE==NULL) ERROR("Error opening file: sd.xxxxx");
        fwrite(SDs,8,N,SDFILE);
        fclose(SDFILE);

        
        /* Clean up memory */
        delete[] LIBB_id;
        delete[] qList;
        delete[] rtDistance;
        delete[] SDs;
        if (rank == 0) cout << "End of ray tracing pre-processor. " << endl;
    }
    // ****************************************************
    MPI_Barrier(comm);
    MPI_Finalize();
    // ****************************************************
}


void IDFunction(int imax, int imin,
                int jmax, int jmin,
                int kmax, int kmin,
                int Nx, int Ny, int Nz,
                double sphere_radius,
                double sphere_center_x,
                double sphere_center_y,
                double sphere_center_z,
                char * LIBB_id)
{
    double* sphere_center; sphere_center = new double[3];
    sphere_center[0] = sphere_center_x;
    sphere_center[1] = sphere_center_y;
    sphere_center[2] = sphere_center_z;
    
    int di = imax-imin;
    int dj = jmax-jmin;
    int dk = kmax-kmin;
    int idx_max = di*dj*dk;
    for (int idx=0;idx<idx_max;idx++)
    {
        int i = idx / (dj*dk);
        int j = (idx - i*dj*dk) / dk;
        int k = idx - j*dk - i*dj*dk;

        i += imin;
        j += jmin;
        k += kmin;

        int m = k*Nx*Ny+j*Nx+i;

        /* ray origin: LBM units */
        double * ray_origin; ray_origin = new double[3];
        ray_origin[0] = double(i);
        ray_origin[1] = double(j);
        ray_origin[2] = double(k);

        double* ray_direction; ray_direction = new double[3];
        ray_direction[0] = 0.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 0.0;
        
        //code for getRayTraceDistance_sp
        double* o_minus_c; o_minus_c = new double[3];
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        double tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        double dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        double omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        double under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        double d1 = 0;
        double d2 = 0;
        double pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        double dat = pos;
        //
        
        if (dat > 0) LIBB_id[m] = 0;
    }//end of idx loop
}

void SignedDistanceFunction(int p,
                            double hx,
                            int Nx, int Ny, int Nz,
                            double min_x, double min_y, double min_z,
                            double* ccx, double* ccy, double* ccz,
                            double* r,
                            double* signedDistance)
{
//    std::cout << "Inside SignedDistanceFunction... \n";
    
    int di = Nx-2;
    int dj = Ny-2;
    int dk = Nz-2;
    int idx_max = di*dj*dk;
    for (int idx=0;idx<idx_max;idx++)
    {
        int i = idx / (dj*dk);
        int j = (idx - i*dj*dk) / dk;
        int k = idx - j*dk - i*dj*dk;

        i += 1; //imin=1
        j += 1; //jmin=1
        k += 1; //kmin=1
        
//        std::cout << "Index (" << i << ", " << j << ", " << k << "), now accessing arrays... \n";
        
        double cx = ccx[p] - min_x;
        double cy = ccy[p] - min_y;
        double cz = ccz[p] - min_z;
        double rad = r[p];

//        std::cout << "Converting sphere center to LBM units \n";
        
        /* converting sphere center to LBM units */
        double* sphere_center; sphere_center = new double[3];
        sphere_center[0] = cx/hx+0.5+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5+0.5;
        sphere_center[2] = cz/hx+0.5+0.5;
        double sphere_radius = rad/hx;

//        std::cout << "Calculating signed distance... \n";
        
        int m = k*Nx*Ny+j*Nx+i;
        double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
        if (dist < signedDistance[m]) signedDistance[m] = dist;
    }//end of idx loop
}

void SignedDistanceFunction2(int p,
                             double hx,
                             int Nx, int Ny, int Nz,
                             double min_x, double min_y, double min_z,
                             double* ccx, double* ccy, double* ccz,
                             double* r,
                             double* signedDistance2)
{
    int di = Nx-2;
    int dj = Ny-2;
    int dk = Nz-2;
    int idx_max = di*dj*dk;
    for (int idx=0;idx<idx_max;idx++)
    {
        int i = idx / (dj*dk);
        int j = (idx - i*dj*dk) / dk;
        int k = idx - j*dk - i*dj*dk;

        i += 1; //imin=1
        j += 1; //jmin=1
        k += 1; //kmin=1

        double cx = ccx[p] - min_x;
        double cy = ccy[p] - min_y;
        double cz = ccz[p] - min_z;
        double rad = r[p];

        /* converting sphere center to LBM units */
        double* sphere_center; sphere_center = new double[3];
        sphere_center[0] = cx/hx+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        double sphere_radius = rad/hx;

        int m = k*Nx*Ny+j*Nx+i;
        double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
        if (dist < signedDistance2[m]) signedDistance2[m] = dist;
    }
}


void QDistancesFunction(int imax, int imin,
                        int jmax, int jmin,
                        int kmax, int kmin,
                        int Nx, int Ny, int Nz,
                        double sphere_radius,
                        double sphere_center_x,
                        double sphere_center_y,
                        double sphere_center_z,
                        double* qDistances)
{
    double* sphere_center; sphere_center = new double[3];
    sphere_center[0] = sphere_center_x;
    sphere_center[1] = sphere_center_y;
    sphere_center[2] = sphere_center_z;
    
    int di = imax-imin;
    int dj = jmax-jmin;
    int dk = kmax-kmin;
    int idx_max = di*dj*dk;
    for (int idx=0;idx<idx_max;idx++){
        int i = idx / (dj*dk);
        int j = (idx - i*dj*dk) / dk;
        int k = idx - j*dk - i*dj*dk;

        i += imin;
        j += jmin;
        k += kmin;
        
        int m = k*Nx*Ny+j*Nx+i;

        /* ray origin: LBM units */
        double* ray_origin; ray_origin = new double[3];
        ray_origin[0] = double(i);
        ray_origin[1] = double(j);
        ray_origin[2] = double(k);

        //initialize, to be used for each direction
        double* ray_direction; ray_direction = new double[3];
        double* o_minus_c; o_minus_c = new double[3];
        double tmp, dp, omc_magSq, under_sqrt, d1, d2, pos, dat;
        
        /* Generate ray trace distance map */
        // 1
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= 1 && dat > 0) qDistances[m] = dat;

        // 2
        ray_direction[0] = -1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= 1 && dat > 0) qDistances[m+Nx*Ny*Nz] = dat;

        // 3
        ray_direction[0] = 0.0;
        ray_direction[1] = 1.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= 1 && dat > 0) qDistances[m+2*Nx*Ny*Nz] = dat;

        // 4
        ray_direction[0] = 0.0;
        ray_direction[1] = -1.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= 1 && dat > 0) qDistances[m+3*Nx*Ny*Nz] = dat;

        // 5
        ray_direction[0] = 0.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= 1 && dat > 0) qDistances[m+4*Nx*Ny*Nz] = dat;

        // 6
        ray_direction[0] = 0.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = -1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= 1 && dat > 0) qDistances[m+5*Nx*Ny*Nz] = dat;

        // 7
        ray_direction[0] = 1.0;
        ray_direction[1] = 1.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+6*Nx*Ny*Nz] = dat;

        // 8
        ray_direction[0] = -1.0;
        ray_direction[1] = -1.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+7*Nx*Ny*Nz] = dat;

        // 9
        ray_direction[0] = 1.0;
        ray_direction[1] = -1.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+8*Nx*Ny*Nz] = dat;

        // 10
        ray_direction[0] = -1.0;
        ray_direction[1] = 1.0;
        ray_direction[2] = 0.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+9*Nx*Ny*Nz] = dat;

        // 11
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+10*Nx*Ny*Nz] = dat;

        // 12
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+11*Nx*Ny*Nz] = dat;

        // 13
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+12*Nx*Ny*Nz] = dat;

        // 14
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+13*Nx*Ny*Nz] = dat;

        // 15
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+14*Nx*Ny*Nz] = dat;

        // 16
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        dat = abs(dat);
        if (dat <= sqrt(2) && dat > 0) qDistances[m+15*Nx*Ny*Nz] = dat;

        // 17
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+16*Nx*Ny*Nz] = dat;

        // 18
        ray_direction[0] = 1.0;
        ray_direction[1] = 0.0;
        ray_direction[2] = 1.0;
        o_minus_c[0] = ray_origin[0] - sphere_center[0];
        o_minus_c[1] = ray_origin[1] - sphere_center[1];
        o_minus_c[2] = ray_origin[2] - sphere_center[2];
        tmp = sqrt(ray_direction[0]*ray_direction[0] + ray_direction[1]*ray_direction[1] + ray_direction[2]*ray_direction[2]);
        if (tmp == 0.0) tmp = 1.0;
        ray_direction[0] /= tmp; ray_direction[1] /= tmp; ray_direction[2] /= tmp;
        dp = ray_direction[0]*o_minus_c[0] + ray_direction[1]*o_minus_c[1] + ray_direction[2]*o_minus_c[2];
        omc_magSq = o_minus_c[0]*o_minus_c[0] + o_minus_c[1]*o_minus_c[1] + o_minus_c[2]*o_minus_c[2];
        under_sqrt = dp*dp - (omc_magSq - sphere_radius*sphere_radius);
        d1 = 0;
        d2 = 0;
        pos = 0;
        if (under_sqrt < 0) pos = 0;
        if (under_sqrt == 0) {
            d1 = -dp;
            pos = d1;
        }
        if (under_sqrt > 0) {
            d1 = -dp + sqrt(under_sqrt);
            d2 = -dp - sqrt(under_sqrt);
            if (d2 >= d1) std::swap(d2,d1);
            if (d2 < 0) d2 = d1;
            pos = d2;
        }
        dat = pos;
        if (dat <= sqrt(2) && dat > 0) qDistances[m+17*Nx*Ny*Nz] = dat;
    }//idx loop
}
