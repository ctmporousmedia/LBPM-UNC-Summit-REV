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
#include "common/Array.h"
#include "common/Domain.h"


inline void PackID(int *list, int count, char *sendbuf, char *ID){
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}
//***************************************************************************************

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
    // Fill in the phase ID values from neighboring processors
    // This unpacks the values once they have been recieved from neighbors
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}



inline void Saturate(char *id, Domain &Dm, int nx, int ny, int nz, int rank)
{
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
    
     /* Walls */
    
    /*
    if (0==Dm.rank()) printf("Adding walls...\n");
    for (int k=0; k<nz; k++){
        for (int j=0; j<ny; j++){
            for (int i=0; i<nx; i++){
                size_t n = k*nx*ny+j*nx+i;
                if (SignDist(i,j,k) < 0.0)  id[n] = 0;
                else{
                    // initially saturated with wetting phase
                    id[n] = 2;
                    if (j + (Ny-2)*(Dm.jproc()) <= 1)  id[n] = 0;
                    if (i + (Nx-2)*(Dm.iproc()) <= 1)  id[n] = 0;
                    if (j + (Ny-2)*(Dm.jproc()) >= (Ny-2)*(nprocy))  id[n] = 0;
                    if (i + (Nx-2)*(Dm.iproc()) >= (Nx-2)*(nprocx))  id[n] = 0;
                }
            }
        }
    }
    
    
    // Add reservoirs for pressure and or flux boundary conditions
    if (0==Dm.rank()) printf("Adding reservoirs...\n");
    for (int k=0; k<nz; k++){
       for (int j=0; j<ny; j++){
           for (int i=0; i<nx; i++){
               n = k*nx*ny+j*nx+i;
               // initially saturated with wetting phase
               if (k + (Nz-2)*(Dm.kproc()) < 4)  id[n] = 2;
               if (k + (Nz-2)*(Dm.kproc()) > (Nz-2)*nprocz-3   ) { id[n] = 2; }
           }
       }
    }
    */
    
    
    for (int k=1; k<nz-1; k++){
        for (int j=1; j<ny-1; j++){
            for (int i=1; i<nx-1; i++){
                size_t n = k*nx*ny+j*nx+i;
                if (id[n] == 2 )  count += 1.0;
            }
        }
    }
    
    
    MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    
    
    // total Global is the number of nodes in the pore-space
    MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    MPI_Allreduce(&maxdist,&maxdistGlobal,1,MPI_DOUBLE,MPI_MAX,Dm.Comm);
    double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
    double porosity=totalGlobal/volume;
    if (rank==0) printf("Media Porosity: %.12f \n",porosity);
    
    
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
    
    MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    double sw_new, sw_diff_new;
    
    porosity = 0;
    
    count = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                n=k*Nx*Ny+j*Nx+i;
                if (id[n] == 2){
                    count+=1.0;
                }
            }
        }
    }
    MPI_Allreduce(&count,&countGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    sw_new = countGlobal/totalGlobal;
    porosity = countGlobal / volume;
    // for test only
    if (rank==0){
        printf("Final saturation=%.12f\n",sw_new);
        printf("Final porosity=%.12f\n",porosity);
    }
    
    
    
    
}



int main(int argc, char **argv)
{

	int rank,nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    {

        /*		bool MULTINPUT=false;

            int NWP,SOLID,rank_offset;
            SOLID=atoi(argv[1]);
            NWP=atoi(argv[2]);

            if (rank==0){
                printf("Solid Label: %i \n",SOLID);
                printf("NWP Label: %i \n",NWP);
            }
            if (argc > 3){
                rank_offset = atoi(argv[3]);
            }
            else{
                MULTINPUT=true;
                rank_offset=0;
            }
         */
        string filename;
        if (argc > 1)
            filename=argv[1];
        else{
            ERROR("lbpm_serial_decomp: no in put database provided \n");
        }
        int rank_offset=0;

        //.......................................................................
        // Reading the domain information file
        //.......................................................................
        int nprocs, nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
        double Lx, Ly, Lz;
        int64_t Nx,Ny,Nz;
        int64_t i,j,k,n;
        int BC=0;
        int64_t xStart,yStart,zStart;
        //  char fluidValue,solidValue;

        xStart=yStart=zStart=0;
        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );

        // Read domain parameters
        auto Filename = domain_db->getScalar<std::string>( "Filename" );
        auto L = domain_db->getVector<double>( "L" );
        auto size = domain_db->getVector<int>( "n" );
        auto SIZE = domain_db->getVector<int>( "N" );
        auto nproc = domain_db->getVector<int>( "nproc" );
        if (domain_db->keyExists( "offset" )){
            auto offset = domain_db->getVector<int>( "offset" );
            xStart = offset[0];
            yStart = offset[1];
            zStart = offset[2];
        }
        auto ReadValues = domain_db->getVector<char>( "ReadValues" );
        auto WriteValues = domain_db->getVector<char>( "WriteValues" );
        auto ReadType = domain_db->getScalar<std::string>( "ReadType" );
        if (ReadType == "8bit"){
        }
        else if (ReadType == "16bit"){
        }
        else{
            printf("INPUT ERROR: Valid ReadType are 8bit, 16bit \n");
            ReadType = "8bit";
        }

        nx = size[0];
        ny = size[1];
        nz = size[2];
        nprocx = nproc[0];
        nprocy = nproc[1];
        nprocz = nproc[2];
        Nx = SIZE[0];
        Ny = SIZE[1];
        Nz = SIZE[2];
        
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
               
               //printf("Set domain \n");
               Nx += 2;
               Nx = Ny = Nz;    // Cubic domain
               int N = Nx*Ny*Nz;
              
        
               // Define Dm.Communication sub-domain -- everywhere
               for (int k=0; k<Nz; k++){
                   for (int j=0; j<Ny; j++){
                       for (int i=0; i<Nx; i++){
                           n = k*Nx*Ny+j*Nx+i;
                           Dm->id[n] = 1;
                       }
                   }
               }
               Dm->CommInit();
        
        for (int k=0; k<Nz; k++){
            for (int j=0; j<Ny; j++){
                for (int i=0; i<Nx; i++){
                    n = k*Nx*Ny+j*Nx+i;
                    Dm->id[n] = 2;
                }
            }
        }
        
        Nx -=2;
        Ny -=2;
        Nz -=2;

        printf("Input media: %s\n",Filename.c_str());
        printf("Relabeling %lu values\n",ReadValues.size());
        for (int idx=0; idx<ReadValues.size(); idx++){
            char oldvalue=ReadValues[idx];
            char newvalue=WriteValues[idx];
            printf("oldvalue=%d, newvalue =%d \n",oldvalue,newvalue);
        }

        nprocs=nprocx*nprocy*nprocz;

        char *SegData = NULL;
        // Rank=0 reads the entire segmented data and distributes to worker processes
        if (rank==0){
            printf("Dimensions of segmented image: %ld x %ld x %ld \n",Nx,Ny,Nz);
            int64_t SIZE = Nx*Ny*Nz;
            SegData = new char[SIZE];
            if (ReadType == "8bit"){
                printf("Reading 8-bit input data \n");
                FILE *SEGDAT = fopen(Filename.c_str(),"rb");
                if (SEGDAT==NULL) ERROR("Error reading segmented data");
                size_t ReadSeg;
                ReadSeg=fread(SegData,1,SIZE,SEGDAT);
                if (ReadSeg != size_t(SIZE)) printf("lbpm_segmented_decomp: Error reading segmented data (rank=%i)\n",rank);
                fclose(SEGDAT);
            }
            else if (ReadType == "16bit"){
                printf("Reading 16-bit input data \n");
                short int *InputData;
                InputData = new short int[SIZE];
                FILE *SEGDAT = fopen(Filename.c_str(),"rb");
                if (SEGDAT==NULL) ERROR("Error reading segmented data");
                size_t ReadSeg;
                ReadSeg=fread(InputData,2,SIZE,SEGDAT);
                if (ReadSeg != size_t(SIZE)) printf("lbpm_segmented_decomp: Error reading segmented data (rank=%i)\n",rank);
                fclose(SEGDAT);
                for (int n=0; n<SIZE; n++){
                    SegData[n] = char(InputData[n]);
                }
            }
            printf("Read segmented data from %s \n",Filename.c_str());
        }

        // Get the rank info
        N = (nx+2)*(ny+2)*(nz+2);

        // number of sites to use for periodic boundary condition transition zone
        int64_t z_transition_size = (nprocz*nz - (Nz - zStart))/2;
        if (z_transition_size < 0) z_transition_size=0;

        char LocalRankFilename[40];
        char *loc_id;
        loc_id = new char [(nx+2)*(ny+2)*(nz+2)];

        std::vector<int> LabelCount(ReadValues.size(),0);
        // Set up the sub-domains
        if (rank==0){
            printf("Distributing subdomains across %i processors \n",nprocs);
            printf("Process grid: %i x %i x %i \n",nprocx,nprocy,nprocz);
            printf("Subdomain size: %i x %i x %i \n",nx,ny,nz);
            printf("Size of transition region: %ld \n", z_transition_size);

            for (int kp=0; kp<nprocz; kp++){
                for (int jp=0; jp<nprocy; jp++){
                    for (int ip=0; ip<nprocx; ip++){
                        // rank of the process that gets this subdomain
                        int rnk = kp*nprocx*nprocy + jp*nprocx + ip;
                        // Pack and send the subdomain for rnk
                        for (k=0;k<nz+2;k++){
                            for (j=0;j<ny+2;j++){
                                for (i=0;i<nx+2;i++){
                                    int64_t x = xStart + ip*nx + i-1;
                                    int64_t y = yStart + jp*ny + j-1;
                                    // int64_t z = zStart + kp*nz + k-1;
                                    int64_t z = zStart + kp*nz + k-1 - z_transition_size;
                                    if (x<xStart) 	x=xStart;
                                    if (!(x<Nx))	x=Nx-1;
                                    if (y<yStart) 	y=yStart;
                                    if (!(y<Ny))	y=Ny-1;
                                    if (z<zStart) 	z=zStart;
                                    if (!(z<Nz))	z=Nz-1;
                                    int64_t nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
                                    int64_t nglobal = z*Nx*Ny+y*Nx+x;
                                    loc_id[nlocal] = SegData[nglobal];
                                }
                            }
                        }
                        // relabel the data
                        for (k=0;k<nz+2;k++){
                            for (j=0;j<ny+2;j++){
                                for (i=0;i<nx+2;i++){
                                    n = k*(nx+2)*(ny+2) + j*(nx+2) + i;;
                                    char locval = loc_id[n];
                                    for (int idx=0; idx<ReadValues.size(); idx++){
                                        char oldvalue=ReadValues[idx];
                                        char newvalue=WriteValues[idx];
                                        if (locval == oldvalue){
                                            loc_id[n] = newvalue;
                                            LabelCount[idx]++;
                                            idx = ReadValues.size();
                                        }
                                    }
                                    //if (loc_id[n]==char(SOLID))     loc_id[n] = 0;
                                    //else if (loc_id[n]==char(NWP))  loc_id[n] = 1;
                                    //else                     loc_id[n] = 2;

                                }
                            }
                        }

                
                        
                        
                        /* Generally Saturate and the MPI bits violates the point of this pp but I dont care right now...*/
                        Saturate(loc_id, *Dm, (nx+2), (ny+2), (nz+2), rank);
                        
                        
                        // Write the data for this rank data
                        sprintf(LocalRankFilename,"ID.%05i",rnk+rank_offset);
                        FILE *ID = fopen(LocalRankFilename,"wb");
                        fwrite(loc_id,1,(nx+2)*(ny+2)*(nz+2),ID);
                        fclose(ID);
                    }
                }
            }
        }
        for (int idx=0; idx<ReadValues.size(); idx++){
            char label=ReadValues[idx];
            int count=LabelCount[idx];
            printf("Label=%d, Count=%d\n",label,count);
        }
        cout <<"Nx=" << Nx << endl;
        int val = 0;
        for (int k=1; k<Nz+2-1; k++){
            for (int j=1; j<Ny+2-1; j++){
                for (int i=1; i<Nx+2-1; i++){
                    n = k*Nx*Ny+j*Nx+i;
                    if (loc_id[n] == 0) val++;
                }
            }
        }
        
        cout << "val=" << val << endl;
        
        int twosphere = 6859 * 2;
        int diff = twosphere - val;
        printf("Solid Diff=%d  (neg means too many solid sites)\n",diff);
        
//        for (int idx=0; idx<ReadValues.size(); idx++){
//            char label=ReadValues[idx];
//            if (label == -1) {
//                int count=LabelCount[idx];
//                int twosphere = 6859 * 2;
//                int diff = twosphere - count;
//                printf("Solid Diff=%d  (neg means too many solid sites)\n",diff);
//            }
//        }
    }
    // ****************************************************
    MPI_Barrier(comm);
    MPI_Finalize();
    // ****************************************************
}
