// Header file for two-phase averaging class
#ifndef TwoPhase_INC
#define TwoPhase_INC

#include <array>
#include <vector>

#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

#include "analysis/interpolator.h"



#include <wrap/io_trimesh/export_off.h>

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

// Volume functionality here 
#include<vcg/complex/algorithms/inertia.h>

#include <vcg/complex/algorithms/convex_hull.h>

#include <vcg/complex/append.h>




class visMyEdge;
class visMyFace;
class visMyVertex;
struct visMyUsedTypes : public vcg::UsedTypes< vcg::Use<visMyVertex>::AsVertexType, vcg::Use<visMyEdge>::AsEdgeType,vcg::Use<visMyFace>::AsFaceType>{};
class visMyVertex:public vcg::Vertex<visMyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDirf, vcg::vertex::Mark>{};
class visMyFace:public vcg::Face< visMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3f, vcg::face::BitFlags, vcg::face::Mark > {};
class visMyEdge:public vcg::Edge<visMyUsedTypes>{};
class visMyMesh:public vcg::tri::TriMesh< std::vector<visMyVertex>, std::vector<visMyFace> , std::vector<visMyEdge>  > {
public:
    visMyMesh() = default;
};




/*
 namespace vcg {
 namespace tri {
 template < class MeshType>         // Classe templatata su Tipo Mesh
 class MiQParametrizer
 {

     typedef typename MeshType::CoordType CoordType;
     typedef typename MeshType::VertexType VertexType;
     typedef typename MeshType::FaceType FaceType;
     typedef typename MeshType::ScalarType ScalarType;
     typedef typename MeshType::VertexType PolyVertexType;

 public:
 */

class nwsMyEdge;
class nwsMyFace;
class nwsMyVertex;
struct nwsMyUsedTypes : public vcg::UsedTypes< vcg::Use<nwsMyVertex>::AsVertexType, vcg::Use<nwsMyEdge>::AsEdgeType,vcg::Use<nwsMyFace>::AsFaceType>{};

class nwsMyVertex:public vcg::Vertex<nwsMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nwsMyFace:public vcg::Face< nwsMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nwsMyEdge:public vcg::Edge<nwsMyUsedTypes>{};

class nwsMyMesh: public vcg::tri::TriMesh< std::vector<nwsMyVertex>, std::vector<nwsMyFace> , std::vector<nwsMyEdge>  > {};



class nwMyEdge;
class nwMyFace;
class nwMyVertex;
struct nwMyUsedTypes : public vcg::UsedTypes< vcg::Use<nwMyVertex>::AsVertexType, vcg::Use<nwMyEdge>::AsEdgeType,vcg::Use<nwMyFace>::AsFaceType>{};

class nwMyVertex:public vcg::Vertex<nwMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nwMyFace:public vcg::Face< nwMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nwMyEdge:public vcg::Edge<nwMyUsedTypes>{};



class nwMyMesh: public vcg::tri::TriMesh< std::vector<nwMyVertex>, std::vector<nwMyFace> , std::vector<nwMyEdge>  > {};



class nwRootEdge;
class nwRootFace;
class nwRootVertex;
struct nwRootUsedTypes : public vcg::UsedTypes< vcg::Use<nwRootVertex>::AsVertexType, vcg::Use<nwRootEdge>::AsEdgeType,vcg::Use<nwRootFace>::AsFaceType>{};

class nwRootVertex:public vcg::Vertex<nwRootUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nwRootFace:public vcg::Face< nwRootUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nwRootEdge:public vcg::Edge<nwRootUsedTypes, vcg::edge::EEAdj>{};



class nwRootMesh: public vcg::tri::TriMesh< std::vector<nwRootVertex>, std::vector<nwRootFace> , std::vector<nwRootEdge>  > {};











class nsMyEdge;
class nsMyFace;
class nsMyVertex;
struct nsMyUsedTypes : public vcg::UsedTypes< vcg::Use<nsMyVertex>::AsVertexType, vcg::Use<nsMyEdge>::AsEdgeType,vcg::Use<nsMyFace>::AsFaceType>{};

class nsMyVertex:public vcg::Vertex<nsMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};

class nsMyFace:public vcg::Face< nsMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nsMyEdge:public vcg::Edge<nsMyUsedTypes>{};

class nsMyMesh:public vcg::tri::TriMesh< std::vector<nsMyVertex>, std::vector<nsMyFace> , std::vector<nsMyEdge>  > {};

class wsMyEdge;
class wsMyFace;
class wsMyVertex;
struct wsMyUsedTypes : public vcg::UsedTypes< vcg::Use<wsMyVertex>::AsVertexType, vcg::Use<wsMyEdge>::AsEdgeType,vcg::Use<wsMyFace>::AsFaceType>{};
class wsMyVertex:public vcg::Vertex<wsMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class wsMyFace:public vcg::Face< wsMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class wsMyEdge:public vcg::Edge<wsMyUsedTypes>{};
class wsMyMesh:public vcg::tri::TriMesh< std::vector<wsMyVertex>, std::vector<wsMyFace> , std::vector<wsMyEdge>  > {};


// ----------------------------

class nMyEdge;
class nMyFace;
class nMyVertex;
struct nMyUsedTypes : public vcg::UsedTypes< vcg::Use<nMyVertex>::AsVertexType, vcg::Use<nMyEdge>::AsEdgeType,vcg::Use<nMyFace>::AsFaceType>{};
class nMyVertex:public vcg::Vertex<nMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nMyFace:public vcg::Face< nMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nMyEdge:public vcg::Edge<nMyUsedTypes>{};
class nMyMesh:public vcg::tri::TriMesh< std::vector<nMyVertex>, std::vector<nMyFace> , std::vector<nMyEdge>  > {};


class wMyEdge;
class wMyFace;
class wMyVertex;
struct wMyUsedTypes : public vcg::UsedTypes< vcg::Use<wMyVertex>::AsVertexType, vcg::Use<wMyEdge>::AsEdgeType,vcg::Use<wMyFace>::AsFaceType>{};
class wMyVertex:public vcg::Vertex<wMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class wMyFace:public vcg::Face< wMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class wMyEdge:public vcg::Edge<wMyUsedTypes>{};
class wMyMesh:public vcg::tri::TriMesh< std::vector<wMyVertex>, std::vector<wMyFace> , std::vector<wMyEdge>  > {};


class sMyEdge;
class sMyFace;
class sMyVertex;
struct sMyUsedTypes : public vcg::UsedTypes< vcg::Use<sMyVertex>::AsVertexType, vcg::Use<sMyEdge>::AsEdgeType,vcg::Use<sMyFace>::AsFaceType>{};
class sMyVertex:public vcg::Vertex<sMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class sMyFace:public vcg::Face< sMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class sMyEdge:public vcg::Edge<sMyUsedTypes>{};
class sMyMesh:public vcg::tri::TriMesh< std::vector<sMyVertex>, std::vector<sMyFace> , std::vector<sMyEdge>  > {};






//class swMyEdge;
//class swMyFace;
//class swMyVertex;
//struct swMyUsedTypes : public vcg::UsedTypes< vcg::Use<swMyVertex>::AsVertexType, vcg::Use<swMyEdge>::AsEdgeType,vcg::Use<swMyFace>::AsFaceType>{};
//class swMyVertex:public vcg::Vertex<swMyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDirf, vcg::vertex::Mark>{};
//class swMyFace:public vcg::Face< swMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3f, vcg::face::BitFlags, vcg::face::Mark > {};
//class swMyEdge:public vcg::Edge<swMyUsedTypes>{};
//class swMyMesh:public vcg::tri::TriMesh< std::vector<swMyVertex>, std::vector<swMyFace> , std::vector<swMyEdge>  > {};


/* CONVEX HULL*/
class chMyEdge;
class chMyFace;
class chMyVertex;
struct chMyUsedTypes : public vcg::UsedTypes< vcg::Use<chMyVertex>::AsVertexType, vcg::Use<chMyEdge>::AsEdgeType,vcg::Use<chMyFace>::AsFaceType>{};
class chMyVertex:public vcg::Vertex<chMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class chMyFace:public vcg::Face< chMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class chMyEdge:public vcg::Edge<chMyUsedTypes>{};
class chMyMesh:public vcg::tri::TriMesh< std::vector<chMyVertex>, std::vector<chMyFace> , std::vector<chMyEdge>  > {};




class periodicMyEdge;
class periodicMyFace;
class periodicMyVertex;
struct periodicMyUsedTypes : public vcg::UsedTypes< vcg::Use<periodicMyVertex>::AsVertexType, vcg::Use<periodicMyEdge>::AsEdgeType,vcg::Use<periodicMyFace>::AsFaceType>{};
class periodicMyVertex:public vcg::Vertex<periodicMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class periodicMyFace:public vcg::Face< periodicMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class periodicMyEdge:public vcg::Edge<periodicMyUsedTypes>{};
class periodicMyMesh:public vcg::tri::TriMesh< std::vector<periodicMyVertex>, std::vector<periodicMyFace> , std::vector<periodicMyEdge>  > {};




class TwoPhase{

public:
    
    int wsGlobal_dS;
    int nsGlobal_dS;
    int nwGlobal_dS;
    
    int wsGlobal_CCs;
    int nsGlobal_CCs;
    int nwGlobal_CCs;
    
	//...........................................................................
	int n_nw_pts,n_ns_pts,n_ws_pts,n_nws_pts,n_local_sol_pts,n_local_nws_pts;
	int n_nw_tris,n_ns_tris,n_ws_tris,n_nws_seg,n_local_sol_tris;
	//...........................................................................
private:
    int nc;
	int kstart,kfinish;

	double fluid_isovalue, solid_isovalue;
	double Volume;
	double nodeVolume;
	double subNodeVolume;
	// initialize lists for vertices for surfaces, common line
	DTMutableList<Point> nw_pts;
	DTMutableList<Point> ns_pts;
	DTMutableList<Point> ws_pts;
	DTMutableList<Point> nws_pts;
	DTMutableList<Point> local_sol_pts;
	DTMutableList<Point> local_nws_pts;
	DTMutableList<Point> tmp;

public:
	// initialize triangle lists for surfaces
    IntArray w_tris;
	IntArray nw_tris;
	IntArray ns_tris;
	IntArray ws_tris;
	IntArray nws_seg;
	IntArray local_sol_tris;

private:
	// Temporary storage arrays
	DoubleArray CubeValues;
	DoubleArray Values;
	DoubleArray DistanceValues;
	DoubleArray KGwns_values;
	DoubleArray KNwns_values;
	DoubleArray InterfaceSpeed;
	DoubleArray NormalVector;

    DoubleArray RecvBuffer;

	char *TempID;

	// CSV / text file where time history of averages is saved
    Array<FILE*> TIMELOG_subdomain;
	FILE *TIMELOG_local;
	FILE *TIMELOG_global;
	FILE *NWPLOG;
	FILE *WPLOG;

public: // Internal data
	//...........................................................................
	std::shared_ptr <Domain> Dm;
	int NumberComponents_WP,NumberComponents_NWP;

public: // local averages (to each MPI process)
        std::shared_ptr<IO::TriList> wn_mesh;//( new IO::TriList() );
        std::shared_ptr<IO::TriList> ns_mesh;//( new IO::TriList() );
        std::shared_ptr<IO::TriList> ws_mesh;//( new IO::TriList() );

        std::shared_ptr<IO::TriList> w_mesh;//( new IO::TriList() );
    size_t npx,npy,npz;
    
	double trimdist;                    // pixel distance to trim surface for specified averages
	double porosity, poreVol;
	double awn, ans, aws, lwns;
	double wp_volume, nwp_volume;
	double As;
	double vol_w, vol_n;                // volumes the exclude the interfacial region
	double pn, pw;                    // local phase averaged pressure
	double efawns;                      // averaged contact angle
    double euler,Kn,Jn,An,eulerW,eulerS;
    double nsEuler,wsEuler,nwEuler;
    double nEuler,wEuler,sEuler;
    double nGlobalEuler, sGlobalEuler, wGlobalEuler;
    double Jwn,Jwn2;                       // average mean curavture - wn interface
    double Kwn,Kwn2;                         // average Gaussian curavture - wn interface
    double JwnW;                         // average mean curavture - wn interface
    double KwnW;                         // average Gaussian curavture - wn interface
	double KNwns;                       // wns common curve normal curavture
	double KGwns;                       // wns common curve geodesic curavture
	double trawn;                       // trimmed interfacial area
	double trJwn;                       // trimmed interfacial area
	double trRwn;                       // trimmed interfacial area
	double wwndnw;
	double wwnsdnwn;
	double Jwnwwndnw;
	double swn;
    
    double Jws;
    double Kws;
    double Jns;
    double Kns;
    
    double JwsS;
    double KwsS;
    double JnsS;
    double KnsS;
    
    DoubleArray qDistances;
    
    int NSpheresOnRank;
    
    bool ComputeAccurateNWPhaseVolume;
    
    
    
    double Jwn2_global;
    double Kwn2_global;
    
    double Ca_global, Re_global;
    double alpha,sphere_diameter,average_rho_input, average_mu,muA,muB;
    size_t nwBorderVert, nsBorderVert, wsBorderVert;
    size_t nwBorderFace, nsBorderFace, wsBorderFace;
    size_t nwBorderVert_global, nsBorderVert_global, wsBorderVert_global;
    size_t nwBorderFace_global, nsBorderFace_global, wsBorderFace_global;
  
	double aw, an, Xn, pwn;
    double awn2, awn2_global;
    
    double igl_mean_curvature, igl_Gaussian_curvature;

	DoubleArray vn;
	DoubleArray vw;
    
    DoubleArray time_avg_vw;
    DoubleArray time_avg_vn;
       
//    DoubleArray no_interface_time_avg_vw;
//    DoubleArray no_interface_time_avg_vn;
//
//    DoubleArray no_interface_vw;
//    DoubleArray no_interface_vn;
    
        DoubleArray nw;
	DoubleArray vnw;
	DoubleArray vws;
	DoubleArray vns;
	DoubleArray vnws;
	DoubleArray Gwn;
	DoubleArray Gns;
	DoubleArray Gws;
	DoubleArray Gwns;

public: // sub-domain averages
	DoubleArray trimdist_sub;               // pixel distance to trim surface for specified averages
	DoubleArray porosity_sub, poreVol_sub;
    DoubleArray awn_sub, ans_sub, aws_sub, lwns_sub, tmp_tmp;
	DoubleArray wp_volume_sub, nwp_volume_sub;
	DoubleArray As_sub;
	DoubleArray vol_w_sub, vol_n_sub;       // volumes the exclude the interfacial region
	DoubleArray pn_sub, pw_sub;           // local phase averaged pressure
	DoubleArray efawns_sub;                 // averaged contact angle
    DoubleArray euler_sub, Kn_sub, Jn_sub, An_sub, eulerW_sub,eulerS_sub;
    
    DoubleArray Jwn_sub;                    // average mean curavture - wn interface
	DoubleArray Kwn_sub;                    // average Gaussian curavture - wn interface
    DoubleArray JwnW_sub;                    // average mean curavture - wn interface
    DoubleArray KwnW_sub;                    // average Gaussian curavture - wn interface
    
    DoubleArray Jws_sub;
    DoubleArray Kws_sub;
    
    DoubleArray Jns_sub;
    DoubleArray Kns_sub;
    
    DoubleArray JwsS_sub;
    DoubleArray KwsS_sub;
    
    DoubleArray JnsS_sub;
    DoubleArray KnsS_sub;

    //IntArray ID;
    
    DoubleArray NxArray;
    DoubleArray NyArray;
    DoubleArray NzArray;
    
    DoubleArray rankSDs;
    
    
    DoubleArray TotalMass;
    DoubleArray DensityA, DensityB;

    DoubleArray Eta_sub; // Interfacial speed
    double Eta;
    double Eta_global;

    DoubleArray Eta_sub1; // Interfacial speed
    double Eta1;
    double Eta_global1;
    

    
	DoubleArray KNwns_sub;                  // wns common curve normal curavture
	DoubleArray KGwns_sub;                  // wns common curve geodesic curavture
	DoubleArray trawn_sub;                  // trimmed interfacial area
	DoubleArray trJwn_sub;                  // trimmed interfacial area
	DoubleArray trRwn_sub;                  // trimmed interfacial area
	DoubleArray wwndnw_sub;
	DoubleArray wwnsdnwn_sub;
	DoubleArray Jwnwwndnw_sub;
	DoubleArray swn_sub;
 	DoubleArray vn_sub;
	DoubleArray vws_sub;
	DoubleArray vns_sub;
	DoubleArray vw_sub;
	DoubleArray vnw_sub;
	DoubleArray vnws_sub;
	DoubleArray Gwn_sub;
	DoubleArray Gns_sub;
	DoubleArray Gws_sub;
	DoubleArray Gwns_sub;
    
    DoubleArray an_sub1;
    DoubleArray an_sub2;
    DoubleArray an_sub3;
    
    double vnwdnw;
    double vnwdnw_global;
    DoubleArray vnwdnw_sub;
    
    DoubleArray wwn;
    DoubleArray wwn_global;
    DoubleArray wwn_sub;
    
    DoubleArray nwOwnO;
    DoubleArray nwOwnO_global;
    DoubleArray nwOwnO_sub;
    
    DoubleArray nwOwnO2;
    DoubleArray nwOwnO_global2;
    DoubleArray nwOwnO_sub2;
    
    DoubleArray nwOwsO;
    DoubleArray nwOwsO_global;
    DoubleArray nwOwsO_sub;

    DoubleArray nwOwsO2;
    DoubleArray nwOwsO_global2;
    DoubleArray nwOwsO_sub2;
   
    DoubleArray nwOnsO2;
    DoubleArray nwOnsO_global2;
    DoubleArray nwOnsO_sub2;
    
    DoubleArray nw_sub;

 
    vector<double> MagVel;
    
    DoubleArray time_avg_vw_sub;
    DoubleArray time_avg_vn_sub;
       
    DoubleArray no_interface_time_avg_vw_sub;
    DoubleArray no_interface_time_avg_vn_sub;
       
    DoubleArray no_interface_vw_sub;
    DoubleArray no_interface_vn_sub;

	DoubleArray aw_sub, an_sub, Xn_sub, pwn_sub;

public: // Global averages (all processes)
    
    int efawns_count;

    double mass_local;
    double mass_global;
    double porosity_global;
	double pn_global,pw_global;               // local phase averaged pressure
	double vol_w_global, vol_n_global;          // volumes the exclude the interfacial region
	double awn_global,ans_global,aws_global;
	double lwns_global;
	double efawns_global;                       // averaged contact angle
    double euler_global,Kn_global,Jn_global,An_global,eulerW_global,eulerS_global;
	double Jwn_global;                          // average mean curavture - wn interface
	double Kwn_global;                          // average Gaussian curavture - wn interface
    double JwnW_global;                          // average mean curavture - wn interface
    double KwnW_global;                          // average Gaussian curavture - wn interface
    
    double affinity;
    
    double Jws_global;
    double Kws_global;                      
    double Jns_global;
    double Kns_global;
    
    double JwsS_global;
    double KwsS_global;
    double JnsS_global;
    double KnsS_global;
    
    double local_nw_phase_volume;
    double global_nw_phase_volume;
    
    double euler_global_NS;
    
    
	double KNwns_global;                        // wns common curve normal curavture
	double KGwns_global;                        // wns common curve geodesic curavture
	double trawn_global;                        // trimmed interfacial area
	double trJwn_global;                        // trimmed interfacial area
	double trRwn_global;                        // trimmed interfacial area
	double nwp_volume_global;                   // volume for the non-wetting phase
	double wp_volume_global;                    // volume for the wetting phase
	double As_global;
	double wwndnw_global;
	double wwnsdnwn_global;
	double Jwnwwndnw_global;
	double swn_global;
	double dEs,dAwn,dAns;
    
    double total_volume_global;
   

        double maxMagVel_global;

	double aw_global, an_global, Xn_global, pwn_global;
						// Global surface energy (calculated by rank=0)
	DoubleArray vn_global;
	DoubleArray vw_global;
    
    DoubleArray press_tmp;
    
    DoubleArray time_avg_vw_global;
    DoubleArray time_avg_vn_global;
    
    DoubleArray no_interface_time_avg_vw_global;
    DoubleArray no_interface_time_avg_vn_global;
    
    DoubleArray no_interface_vw_global;
    DoubleArray no_interface_vn_global;
    
        DoubleArray nw_global;

	DoubleArray vnw_global;
	DoubleArray vws_global;
	DoubleArray vns_global;
	DoubleArray vnws_global;
	DoubleArray Gwn_global;
	DoubleArray Gns_global;
	DoubleArray Gws_global;
	DoubleArray Gwns_global;

public: // Internal data
    std::array<size_t,3> subdivide;
	size_t Nx,Ny,Nz;
	IntArray PhaseID;	// Phase ID array (solid=0, non-wetting=1, wetting=2)
    IntArray ID;
	BlobIDArray Label_WP;   // Wetting phase label
	BlobIDArray Label_NWP;  // Non-wetting phase label index (0:nblobs-1)
	std::vector<BlobIDType> Label_NWP_map;  // Non-wetting phase label for each index
	DoubleArray SDn;
    DoubleArray TempPhase;
    DoubleArray TempPhase2;
    DoubleArray ParticularSDs;
    
    
	DoubleArray SDw;
	DoubleArray SDs;
	DoubleArray Phase;
	DoubleArray Press;
	DoubleArray dPdt;
    DoubleArray PhiField;
    DoubleArray BlobID;
    DoubleArray MeanCurvature;
	DoubleArray GaussCurvature;
    DoubleArray MeanCurvatureW;
    DoubleArray GaussCurvatureW;
    DoubleArray VFPhase;
    
	DoubleArray SDs_x;		// Gradient of the signed distance
	DoubleArray SDs_y;
	DoubleArray SDs_z;
	DoubleArray SDn_x;		// Gradient of the signed distance
	DoubleArray SDn_y;
	DoubleArray SDn_z;
    DoubleArray SDw_x;              // Gradient of the signed distance
    DoubleArray SDw_y;
    DoubleArray SDw_z;

    int method = 0;

   DoubleArray normalw_x;
   DoubleArray normalw_y;
   DoubleArray normalw_z;
    
    IntArray BoolField;

	DoubleArray DelPhi;
    DoubleArray zero;
    
    DoubleArray VFmask, VFboundary, vfmask;
    
    DoubleArray wnCube;
    DoubleArray wsCube;
    DoubleArray nsCube;
    
    DoubleArray InactiveCellFluidCenterOfMass;
    
    IntArray InactiveCellsActiveNeighborInterpolationAbscissa;
    
    DoubleArray ParticularSDs_x;
    DoubleArray ParticularSDs_y;
    DoubleArray ParticularSDs_z;
    
    DoubleArray BCarray;
    DoubleArray Aarray;
    DoubleArray Darray;
    DoubleArray Qarray;
    
    DoubleArray BCarray2;
    DoubleArray Aarray2;
    DoubleArray Darray2;
    DoubleArray Qarray2;
    
    DoubleArray GradPhiX;
    DoubleArray GradPhiY;
    DoubleArray GradPhiZ;
	
    DoubleArray temp_tplus;
    DoubleArray SDw_tplus;
    DoubleArray SDw_tplusplus;
    DoubleArray SDn_tplus;
    DoubleArray SDn_tplusplus;
	DoubleArray Phase_tplus;
	DoubleArray Phase_tminus;
	DoubleArray Vel_x;		// Velocity
	DoubleArray Vel_y;
	DoubleArray Vel_z;
    /* Time-averaged velocities */
    DoubleArray Time_Avg_Vel_x;
    DoubleArray Time_Avg_Vel_y;
    DoubleArray Time_Avg_Vel_z;
	//	Container for averages;
	DoubleArray ComponentAverages_WP;
	DoubleArray ComponentAverages_NWP;

    DoubleArray IntegrationRegionOne;
    
    int numToAvg;

public: // Public interfaces
	TwoPhase(std::shared_ptr <Domain> Dm, std::array<size_t,3> subdivide = { 1, 1, 1 } );
	~TwoPhase();
	void Initialize();
//	void SetupCubes(Domain &Dm);
    
	void UpdateMeshValues();
	void UpdateSolid();
	void ComputeDelPhi();
	void ComputeSDw();
	void ColorToSignedDistance(double Beta, DoubleArray &ColorData, DoubleArray &DistData);
    void SmoothSameSign(double Beta, DoubleArray &ColorData, DoubleArray &DistData);
    void MedianFilter(double Beta, DoubleArray &ColorData, DoubleArray &DistData);
	void ComputeLocal();
    
    void PosAdjustSDs();
    void NegAdjustSDs();
    
    void CuttingPlaneZ(int zval);
    
    void MeshComputations();
    void ZeroSolidPhaseLatticeSites(DoubleArray &Field);
    
    void MeanFilter(DoubleArray &Field);
    void InactiveCellMeanFilter(DoubleArray &Field);
    void ComputePhiDistance(DoubleArray &Field);
    
    void FindIsolatedLatticeSites();

    void ToBinary();
    
    void InitializeRankSignDistance();
    void RankSignDistance();
    void PhaseInit(double affinity);
    void FluidSolidVolumeFraction(double affinity);
    void FluidSolidVolumeFraction_Analysis();
    
    void FluidCentroidInit(double affinity);
    
	void AssignComponentLabels();
	void ComponentAverages();
	inline void Reduce() { ReduceGlobal(); }
	void WriteSurfaces(int logcount);
	void NonDimensionalize(double D, double viscosity, double IFT);
    double SignedVolumeOfTriangle(Point p1, Point p2, Point p3);
	void PrintAll(int timestep);
	int GetCubeLabel(int i, int j, int k, IntArray &BlobLabel);
    
    void nwComputeQuantities(nwMyMesh & /* some tri mesh */);
    void nsComputeQuantities(nsMyMesh & /* some tri mesh */);
    void wsComputeQuantities(wsMyMesh & /* some tri mesh */);

    void nwsComputeQuantities(nwsMyMesh & /* some tri mesh */);
    
    void nComputeQuantities(nMyMesh & /* some tri mesh */);
    void sComputeQuantities(sMyMesh & /* some tri mesh */);
    void wComputeQuantities(wMyMesh & /* some tri mesh */);
    
    void ComputeAccurateNWInterfaceCurvatures(nwMyMesh & /* some tri mesh */);
    
	void SortBlobs();
	void PrintComponents(int timestep);
    void SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha);
    double Volume_w(){
        return wp_volume_global;
    }
    double Volume_n(){
        return nwp_volume_global;
    }

    void CreateIntegrationBoundaryRegionOne();
    void CreateIntegrationBoundaryRegionTwo(DoubleArray & Field);
    void AggregateLabels( const std::string& filename );

private: // Private interfaces
	void ReduceGlobal();
};


#endif

