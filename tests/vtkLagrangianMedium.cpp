#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>


#include <math.h>

#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/Database.h"
#include "common/SpherePack.h"

//#include "common/pmmc.h"
#include "common/Domain.h"
#include "common/SpherePack.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"


/* Exporting meshes for vcg */
#include <wrap/io_trimesh/export_off.h>
#include <wrap/io_trimesh/import_off.h>

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>     //class UpdateNormals
#include <vcg/complex/algorithms/update/curvature.h>   //class UpdateCurvature
#include <vcg/complex/algorithms/update/curvature_fitting.h>   //class UpdateCurvature
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>



using namespace std;


class solidMyEdge;
class solidMyFace;
class solidMyVertex;
struct solidMyUsedTypes : public vcg::UsedTypes< vcg::Use<solidMyVertex>::AsVertexType, vcg::Use<solidMyEdge>::AsEdgeType,vcg::Use<solidMyFace>::AsFaceType>{};
class solidMyVertex:public vcg::Vertex<solidMyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDirf, vcg::vertex::Mark>{};
class solidMyFace:public vcg::Face< solidMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3f, vcg::face::BitFlags, vcg::face::Mark > {};
class solidMyEdge:public vcg::Edge<solidMyUsedTypes>{};
class solidMyMesh:public vcg::tri::TriMesh< std::vector<solidMyVertex>, std::vector<solidMyFace> , std::vector<solidMyEdge>  > {};

// Inline function to read line without a return argument
static inline void fgetl( char * str, int num, FILE * stream )
{
  char* ptr = fgets( str, num, stream );
  if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}

void ReadLagrangianPoints(int nspheres, double *List_cx, double *List_cy, double *List_cz)
{
    // Read in the full sphere pack
    //...... READ IN THE SPHERES...................................
    cout << "Reading the lagrangian pack file..." << endl;
    FILE *fid = fopen("pack.off","rb");
    INSIST(fid!=NULL,"Error opening pack.off");
    //.........Trash the header lines..........
    char line[100];
    fgetl(line, 100, fid); // should say: OFF
    fgetl(line, 100, fid); // contains the number of vertices, faces, edges
    char* dat = line;
    int nvertices = strtod(dat,&dat);
    int nfaces = strtod(dat,&dat);
    cout << "nvertices=" << nvertices << " nfaces=" << nfaces << endl;
    
    //........read the spheres..................
    // We will read until a blank like or end-of-file is reached
    int count = 0;
    while ( !feof(fid) && fgets(line,100,fid)!=NULL && count < nvertices ) {
        char* line2 = line;
        List_cx[count] = strtod(line2,&line2);
        List_cy[count] = strtod(line2,&line2);
        List_cz[count] = strtod(line2,&line2);
       // List_rad[count] = strtod(line2,&line2);
        count++;
    }
    cout << "Number of lines extracted is: " << count << endl;
    //INSIST( count==nspheres, "Specified number of spheres is probably incorrect!" );
    // .............................................................
}


//           /* LAGRANGIAN POINTS */
//           double domain_length = 0.08682328464463308;
//           solidMyMesh lagrangianMesh;
//           vcg::tri::io::ImporterOFF<solidMyMesh>::Open(lagrangianMesh,"pack.off", 0);
//
//           cout << "num vertices=" << lagrangianMesh.VN() << " num faces=" << lagrangianMesh.FN() << endl;
//           double x,y,z;
//
//           for(auto vi = lagrangianMesh.vert.begin(); vi != lagrangianMesh.vert.end(); ++vi) {
//               cout << *vi[1] << " " << *vi[2] << " " << *vi[3] << " " << endl;
//
//
//
//           }

struct Point {
    size_t n; // local n
    double cx; // global cx
    double cy; // global cy
    double cz; // global cz
} ;

void LagrangianToEulerianMap(char *id, int nvertices, double *List_cx, double *List_cy, double *List_cz,
        double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, vector<Point> &P)
{
    // Use sphere lists to determine which nodes are in porespace
    // Write out binary file for nodes
    double cx,cy,cz;
    int i,j,k,n;
    //............................................
    int i2,j2,k2;
    
    for (n = 0; n < Nx*Ny*Nz; n++) id[n] = 0;

    // .........Loop over the spheres.............
    size_t n2 = 0;
    for (size_t p=0;p<nvertices;p++){
        // Get the sphere from the list, map to local min
        cx = List_cx[p];  cy = List_cy[p];  cz = List_cz[p];
        
        /* Get global indices */
        i2 = int( (cx/Lx + iproc)*double(Nx-2) ) + 1;
        j2 = int( (cy/Ly + jproc)*double(Ny-2) ) + 1;
        k2 = int( (cz/Lz + kproc)*double(Nz-2) ) + 1;
        // cout << "i2=" << i2 << " j2=" << j2 << " k2=" << k2 << endl;
        n2 = k2*(Nx-2)*(Ny-2)+j2*(Nx-2)+i2;
        id[n2] += 1;
        P.at(p).n = n2; P.at(p).cx = cx; P.at(p).cy = cy; P.at(p).cz = cz;
    }

    // end of function
}






int main(int argc, char **argv)
{
     int rank,nprocs;
       MPI_Init(&argc,&argv);
       MPI_Comm comm = MPI_COMM_WORLD;
       MPI_Comm_rank(comm,&rank);
       MPI_Comm_size(comm,&nprocs);
       {
           // parallel domain size (# of sub-domains)
           int nprocx,nprocy,nprocz;
           int iproc,jproc,kproc;
           //*****************************************
           
           if (rank == 0){
               printf("********************************************************\n");
               printf("Running Sphere Packing pre-processor for LBPM-WIA    \n");
               printf("********************************************************\n");
           }
           
           // Variables that specify the computational domain
           string filename;
           int Nx,Ny,Nz;        // local sub-domain size
           int nspheres;        // number of spheres in the packing
           double Lx,Ly,Lz;    // Domain length
           double D = 1.0;        // reference length for non-dimensionalization
           
           int i,j,k,n;
           filename = argv[1];
           auto db = std::make_shared<Database>( filename );
           auto domain_db = db->getDatabase( "Domain" );
           
           auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
           Nx = Dm->Nx;
           Ny = Dm->Ny;
           Nz = Dm->Nz;
           Lx = Dm->Lx;
           Ly = Dm->Ly;
           Lz = Dm->Lz;
           iproc = Dm->iproc();
           jproc = Dm->jproc();
           kproc = Dm->kproc();
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
           
           //.......................................................................
           // Filenames used
           char LocalRankString[8];
           char LocalRankFilename[40];
           char LocalRestartFile[40];
           sprintf(LocalRankString,"%05d",rank);
           sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
           sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
           
           //    printf("Local File Name =  %s \n",LocalRankFilename);
            // .......... READ THE INPUT FILE .......................................
            //    char value;
            char *id;
            id = new char[N];
            int sum = 0;
            double sum_local = 0;
            // double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
            // Consider reservoirs
            double iVol_global = (1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);  // got rid of the inverse to see the number
            printf("iVol_global=%f\n",iVol_global);
            double porosity, pore_vol;
            //...........................................................................
            DoubleArray SignDist(Nx,Ny,Nz);
            //.......................................................................
           cout << "Reading the lagrangian pack file..." << endl;
           FILE *fid = fopen("pack.off","rb");
           INSIST(fid!=NULL,"Error opening pack.off");
           //.........Trash the header lines..........
           char line[100];
           fgetl(line, 100, fid); // should say: OFF
           fgetl(line, 100, fid); // contains the number of vertices, faces, edges
           char* dat = line;
           int nvertices = strtod(dat,&dat);
           int nfaces = strtod(dat,&dat);
           fclose(fid);
           //........................................................................
           double *cx,*cy,*cz;
           cx = new double[nvertices];
           cy = new double[nvertices];
           cz = new double[nvertices];
           //........................................................................
            MPI_Barrier(comm);
            // Broadcast the sphere packing to all processes
            MPI_Bcast(cx,nvertices,MPI_DOUBLE,0,comm);
            MPI_Bcast(cy,nvertices,MPI_DOUBLE,0,comm);
            MPI_Bcast(cz,nvertices,MPI_DOUBLE,0,comm);
            //.......................................................................
            MPI_Barrier(comm);
            if (rank == 0)  ReadLagrangianPoints(nvertices,cx,cy,cz);
            // if (rank == 0) cout << "Domain set." << endl;

           vector<Point> points;
            points.resize(nvertices);
           
            //.......................................................................
            LagrangianToEulerianMap(id,nvertices,cx,cy,cz,Lx,Ly,Lz,Nx,Ny,Nz,
                           iproc,jproc,kproc,nprocx,nprocy,nprocz,points);
            //.......................................................................
            // Assign the phase ID field based on the signed distance
            //.......................................................................
           for (int i = 0; i < 10; i++)
               cout << i << " Px=" << points.at(i).cx << " Py=" << points.at(i).cy << " Pz=" << points.at(i).cz << endl;
           
           
           for (k=1;k<Nz-1;k++){ printf("k=%d\n",k);
               for (j=1;j<Ny-1;j++){
                   for (i=1;i<Nx-1;i++){
                       n = k*(Nx-2)*(Ny-2)+j*(Nx-2)+i;
                       if (id[n] > 9) printf("%d",9);
                       else printf("%d",id[n]);
                   }
                   printf("\n");
               }
           }
           
           
           /*
           for (k=1;k<Nz-1;k++){ printf("k=%d\n",k);
               for (j=1;j<Ny-1;j++){
                   for (i=1;i<Nx-1;i++){
                       n = k*(Nx-2)*(Ny-2)+j*(Nx-2)+i;
                       if (id[n] > 9) printf("%d",9);
                       else printf("%d",id[n]);
                   }
                   printf("\n");
               }
           }
           */
           
           
           double Jwn_global, Kwn_global;
     
           
           
           
           int npx,npy,npz;
           npx = npy = npz = 1;
           
           
           
        solidMyMesh m_wn;
           vcg::tri::io::ImporterOFF<solidMyMesh>::Open(m_wn,"pack.off", 0);
           
           /* First steps are to delete any duplicate faces and vertices */
        
           int wnGlobalEuler = 0;
           
           size_t beforeVN = m_wn.VN();

           size_t beforeFN = m_wn.FN();
           int wnEN = 0;
           size_t wnAfterVN, wnAfterFN;
           
           
           /* Get all edges */
              vcg::tri::UpdateTopology<solidMyMesh>::AllocateEdge(m_wn);
           
//           vcg::tri::Clean<solidMyMesh>::RemoveDuplicateFace(m_wn);
//           vcg::tri::Clean<solidMyMesh>::RemoveDuplicateVertex(m_wn,false);
//           vcg::tri::Clean<solidMyMesh>::RemoveDuplicateEdge(m_wn);

           
           /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
//           int degenface = vcg::tri::Clean<solidMyMesh>::RemoveDegenerateFace(m_wn);
//           int degenvertex = vcg::tri::Clean<solidMyMesh>::RemoveDegenerateVertex(m_wn);
//           int degenedge = vcg::tri::Clean<solidMyMesh>::RemoveDegenerateEdge(m_wn);
              
              /* Final step: Faces, vertices, and edges must be compact for curvature calculation */
              vcg::tri::Allocator<solidMyMesh>::CompactEveryVector(m_wn);
           
            wnAfterVN = m_wn.VN(); wnAfterFN = m_wn.FN();
                wnEN = m_wn.EN();
           

               /* Hope that the damn mesh is finally OK... */
               vcg::tri::UpdateCurvatureFitting<solidMyMesh>::computeCurvature(m_wn);
               
               /* Fill global values */
               size_t wnGlobalEdges = 0;
               size_t wnGlobalFaces = 0;
               size_t wnGlobalVertices = 0;
               MPI_Allreduce(&wnEN,&wnGlobalEdges,1,MPI_INT,MPI_SUM,Dm->Comm);
               MPI_Allreduce(&wnAfterVN,&wnGlobalVertices,1,MPI_INT,MPI_SUM,Dm->Comm);
               MPI_Allreduce(&wnAfterFN,&wnGlobalFaces,1,MPI_INT,MPI_SUM,Dm->Comm);
               
               
               
               /* Global Euler characteristic */
               wnGlobalEuler = wnGlobalVertices - wnGlobalEdges + wnGlobalFaces;
            if (0 == Dm->rank()) std::cout << "wn Euler characteristic:" << wnGlobalEuler << std::endl;
               
               int wnNonZero = 0;
               int wnGlobalNonZero = 0;
               if (wnAfterVN > 0) wnNonZero = 1;
               MPI_Allreduce(&wnNonZero,&wnGlobalNonZero,1,MPI_INT,MPI_SUM,Dm->Comm);
               
               
               double area = vcg::tri::Stat<solidMyMesh>::ComputeMeshArea(m_wn);
               double area_global = 0.0;
               MPI_Allreduce(&area,&area_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
           if (0 == Dm->rank()) std::cout << "wn area_global:" << area_global << std::endl;
               
               /* mean and Gaussian curvature calculations */
               double k1 = 0.0;
               double k2 = 0.0;
               double k1_global = 0;
               double k2_global = 0;
               
               int k1nan = 0;
               int k2nan = 0;
               double tmp = 0;
               for(auto vi = m_wn.vert.begin(); vi != m_wn.vert.end(); ++vi) {
                   tmp = (*vi).K1();
                   if (isnan(tmp)) {k1nan+=1;}
                   else k1 += tmp;
                   tmp = (*vi).K2();
                   if (isnan(tmp)) {k2nan+=1;}
                   else k2 += tmp;
               }
               std::cout << "(" << Dm->rank() << ") WN: k1nan=" << k1nan << " k2nan=" << k2nan << " k1,k2=" << k1 << " " << k2 <<  std::endl;
               
               double hval = 0;
               double gval = 0;
               
               hval = 0.5*(k1+k2);
               gval = k1*k2;
               
               // avg k1 and k2 on a process
               if (m_wn.VN() > 0) k1 = 2*hval/double(m_wn.VN()); else k1 = 0;
               if (m_wn.VN() > 0) k2 = gval/double(m_wn.VN())/double(m_wn.VN()); else k2 = 0;
                   
               MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
               MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
               
               k1_global/=double(npx*npy*npz);
               k2_global/=double(npx*npy*npz);
               
               Jwn_global = k1_global;
               Kwn_global = k2_global;
               
               if (0 == Dm->rank()) std::cout << "Jwn_global=" << Jwn_global << " Kwn_global=" << Kwn_global << std::endl;
           //    tmp = Kwn_global  / (2*M_PI); // * area_global
           //    if (0 == Dm->rank())
           
           
           
           /*
            for (k=0;k<Nz;k++){
                for (j=0;j<Ny;j++){
                    for (i=0;i<Nx;i++){
                        n = k*Nx*Ny+j*Nx+i;
                        id[n] = 0;
                    }
                }
            }
            sum=0;
            pore_vol = 0.0;
            for ( k=1;k<Nz-1;k++){
                for ( j=1;j<Ny-1;j++){
                    for ( i=1;i<Nx-1;i++){
                        n = k*Nx*Ny+j*Nx+i;
                        if (SignDist(n) > 0.0){
                            id[n] = 2;
                   
                        }
                       
                    }
                }
            }
            
            int aa = 0;
           // if (0 == Dm->rank()) {
                for ( k=1;k<Nz-1;k++){ //printf("k=%d\n",k);
                    for ( j=1;j<Ny-1;j++){
                        for ( i=1;i<Nx-1;i++){
                            n = k*Nx*Ny+j*Nx+i;
                            aa = id[n];
                           // printf("%d ",aa);
                            if (aa == 2) sum++;
                        }
                     //   printf("\n");
                    }
                }
                
            //}
            
            
            
            
           
            porosity = 0;
            sum_local = double(sum);
            //printf("sum_local=%f",sum_local);
            MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
            //printf("porosity=%f\n",porosity);
            porosity /= iVol_global;
            // if (0 == Dm->rank()) printf("(%d) Media porosity before walls and reservoirs=%f\n",Dm->rank(),porosity);
            
            // Run Morphological opening to initialize 50% saturation
            double SW=0.00;
            // if (rank==0) printf("MorphOpen: Initializing with saturation %f \n",SW);
            // if (rank==0) printf("Saturate: Initializing with saturation %f \n",SW);
            // Saturate(SignDist, id, *Dm, Nx, Ny, Nz, rank, SW);
            // MorphOpen(SignDist, id, *Dm, Nx, Ny, Nz, rank, SW);
            
            
           
            //.......................................................................
            sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
            FILE *IDFILE = fopen(LocalRankFilename,"wb");
            if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
            fwrite(id,1,N,IDFILE);
            fclose(IDFILE);
            //......................................................................
           
           */
    
    
    
    
  
    
   
       
            
    
        
    }
    // ****************************************************
    MPI_Barrier(comm);
    MPI_Finalize();
    // ****************************************************
}
