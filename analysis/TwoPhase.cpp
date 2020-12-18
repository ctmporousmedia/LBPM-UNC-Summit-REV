#include "analysis/TwoPhase.h"

#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "what is happening here?"
#endif



//#include <Eigen/Dense>
//#include <igl/principal_curvature.h>
//#include <igl/cotmatrix.h>
//#include <igl/massmatrix.h>
//#include <igl/invert_diag.h>
//
//#include <igl/avg_edge_length.h>
//#include <igl/barycenter.h>
//#include <igl/grad.h>

// #include </Users/cpf/Repos/libigl/include/igl/cotmatrix.h>
//#include <igl/massmatrix.h>
//#include <igl/invert_diag.h>
//
//#include <igl/gaussian_curvature.h>
//
//
//#include <igl/readOFF.h>
////#include <igl/opengl/glfw/Viewer.h>
//#include <igl/jet.h>
//#include <igl/writeOFF.h>


#include <math.h>

#include <iostream>
#include <fstream>
using namespace std;

#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/Array.hpp"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

#include <iterator>
#include <vector>

#include "analysis/distance.h"



#define BLOB_AVG_COUNT 40

// Array access for averages defined by the following
#define VOL 0
#define TRIMVOL 1
#define PRS 2
#define AWN 3
#define AWNS 4
#define LWNS 5
#define JWN 6
#define KWN 7
#define CWNS 8
#define KNWNS 9
#define KGWNS 10
#define VX 11
#define VY 12
#define VZ 13
#define VSQ 14
#define VNWX 15
#define VNWY 16
#define VNWZ 17
#define VNWSX 18
#define VNWSY 19
#define VNWSZ 20
#define GWNXX 21
#define GWNYY 22
#define GWNZZ 23
#define GWNXY 24
#define GWNXZ 25
#define GWNYZ 26
#define TRAWN 27
#define TRJWN 28
#define CMX 29
#define CMY 30
#define CMZ 31
#define EULER 32
#define INTCURV 33
#define GXX 34
#define GYY 35
#define GZZ 36
#define GXY 37
#define GXZ 38
#define GYZ 39

#define PI 3.14159265359


// Constructor
TwoPhase::TwoPhase(std::shared_ptr <Domain> dm, std::array<size_t,3> sub ):
n_nw_pts(0), n_ns_pts(0), n_ws_pts(0), n_nws_pts(0), n_local_sol_pts(0), n_local_nws_pts(0),
n_nw_tris(0), n_ns_tris(0), n_ws_tris(0), n_nws_seg(0), n_local_sol_tris(0),
nc(0), kstart(0), kfinish(0), fluid_isovalue(0), solid_isovalue(0),    Volume(0),nodeVolume(0),subNodeVolume(0),
TIMELOG_local(nullptr), TIMELOG_global(nullptr), NWPLOG(nullptr), WPLOG(nullptr),
Dm(dm), NumberComponents_WP(0), NumberComponents_NWP(0), trimdist(0),
porosity(0), poreVol(0), awn(0), ans(0), aws(0), lwns(0), aw(0), an(0), Xn(0), pwn(0), wp_volume(0), nwp_volume(0),
As(0), vol_w(0), vol_n(0), vnwdnw(0), vnwdnw_global(0), Jws(0), Kws(0), euler(0), euler_global(0), eulerW(0), eulerW_global(0),
pn(0), pw(0), pn_global(0), pw_global(0), swn(0), swn_global(0), vol_w_global(0), vol_n_global(0),
awn_global(0), ans_global(0),aws_global(0), lwns_global(0), efawns(0), efawns_global(0), JnsS(0), JnsS_global(0), KnsS(0),KnsS_global(0),
JwsS(0),JwsS_global(0), KwsS(0),KwsS_global(0),eulerS(0),eulerS_global(0),
Jwn(0), Jwn_global(0),  Kwn(0), Kwn_global(0), KNwns(0), KNwns_global(0), JwnW_global(0), KwnW_global(0), JwnW(0), KwnW(0),
KGwns(0), KGwns_global(0), trawn(0), trawn_global(0), trJwn(0), trJwn_global(0),
trRwn(0), trRwn_global(0), nwp_volume_global(0), wp_volume_global(0), Kws_global(0), Jws_global(0), Kns_global(0), Jns_global(0),
As_global(0), wwndnw_global(0), wwnsdnwn_global(0), Jwnwwndnw_global(0), dEs(0), dAwn(0), dAns(0),
subdivide( sub )
{
    Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;

    npx = Dm->nprocx();
    npy = Dm->nprocy();
    npz = Dm->nprocz();
    Volume=1.0*double((Nx-2)*(Ny-2)*(Nz-2)*npx*npy*npz);

    nodeVolume = 1.0*double((Nx-2)*(Ny-2)*(Nz-2));
    
    temp_tplus.resize(Nx,Ny,Nz);   temp_tplus.fill(0);
    
    zero.resize(Nx,Ny,Nz);         zero.fill(0);
    VFmask.resize(Nx,Ny,Nz);       VFmask.fill(0);
    
    DensityA.resize(Nx,Ny,Nz);     DensityA.fill(0);
    DensityB.resize(Nx,Ny,Nz);     DensityB.fill(0);
    
    PhiField.resize(Nx,Ny,Nz);     PhiField.fill(0);
    ID.resize(Nx,Ny,Nz);           ID.fill(0);
    Label_WP.resize(Nx,Ny,Nz);      Label_WP.fill(0);
    Label_NWP.resize(Nx,Ny,Nz);     Label_NWP.fill(0);
    PhaseID.resize(Nx,Ny,Nz);       PhaseID.fill(0);
    NxArray.resize(Nx,Ny,Nz);      NxArray.fill(0);
    NyArray.resize(Nx,Ny,Nz);      NyArray.fill(0);
    NzArray.resize(Nx,Ny,Nz);      NzArray.fill(0);
    
    SDn.resize(Nx,Ny,Nz);           SDn.fill(0);
    SDw.resize(Nx,Ny,Nz);           SDw.fill(0);
    
    SDn_tplus.resize(Nx,Ny,Nz);           SDn_tplus.fill(0);
    SDw_tplus.resize(Nx,Ny,Nz);           SDw_tplus.fill(0);

    SDw_tplusplus.resize(Nx,Ny,Nz);           SDw_tplusplus.fill(0);
    SDn_tplusplus.resize(Nx,Ny,Nz);           SDn_tplusplus.fill(0);
    
    BCarray.resize(Nx,Ny,Nz); BCarray.fill(0);
    Aarray.resize(Nx,Ny,Nz);  Aarray.fill(0);
    Darray.resize(Nx,Ny,Nz);  Darray.fill(0);
    Qarray.resize(Nx,Ny,Nz);  Qarray.fill(0);
    
    BCarray2.resize(Nx,Ny,Nz); BCarray2.fill(0);
    Aarray2.resize(Nx,Ny,Nz);  Aarray2.fill(0);
    Darray2.resize(Nx,Ny,Nz);  Darray2.fill(0);
    Qarray2.resize(Nx,Ny,Nz);  Qarray2.fill(0);
    
    VFPhase.resize(Nx,Ny,Nz);   VFPhase.fill(0);
    
    TempPhase.resize(Nx,Ny,Nz);     TempPhase.fill(0);
    
    SDs.resize(Nx,Ny,Nz);           SDs.fill(0);
    Phase.resize(Nx,Ny,Nz);         Phase.fill(0);
    Press.resize(Nx,Ny,Nz);         Press.fill(0);
    dPdt.resize(Nx,Ny,Nz);          dPdt.fill(0);

    MeanCurvature.resize(Nx,Ny,Nz); MeanCurvature.fill(0);
    GaussCurvature.resize(Nx,Ny,Nz); GaussCurvature.fill(0);
    
    MeanCurvatureW.resize(Nx,Ny,Nz); MeanCurvatureW.fill(0);
    GaussCurvatureW.resize(Nx,Ny,Nz); GaussCurvatureW.fill(0);
    
    SDs_x.resize(Nx,Ny,Nz);         SDs_x.fill(0);      // Gradient of the signed distance
    SDs_y.resize(Nx,Ny,Nz);         SDs_y.fill(0);
    SDs_z.resize(Nx,Ny,Nz);         SDs_z.fill(0);
    SDn_x.resize(Nx,Ny,Nz);         SDn_x.fill(0);      // Gradient of the signed distance
    SDn_y.resize(Nx,Ny,Nz);         SDn_y.fill(0);
    SDn_z.resize(Nx,Ny,Nz);         SDn_z.fill(0);
    
    SDw_x.resize(Nx,Ny,Nz);         SDw_x.fill(0);      // Gradient of the signed distance
    SDw_y.resize(Nx,Ny,Nz);         SDw_y.fill(0);
    SDw_z.resize(Nx,Ny,Nz);         SDw_z.fill(0);
    
    
    
    normalw_x.resize(Nx,Ny,Nz);         normalw_x.fill(0);      // Gradient of the signed distance
    normalw_y.resize(Nx,Ny,Nz);         normalw_y.fill(0);      // Gradient of the signed distance
    normalw_z.resize(Nx,Ny,Nz);         normalw_z.fill(0);      // Gradient of the signed distance
    
    
    
    MagVel.resize(Nx*Ny*Nz);   std::fill(MagVel.begin(), MagVel.end(), 0);
    
    
    DelPhi.resize(Nx,Ny,Nz);        DelPhi.fill(0);
    Phase_tplus.resize(Nx,Ny,Nz);   Phase_tplus.fill(0);
    Phase_tminus.resize(Nx,Ny,Nz);  Phase_tminus.fill(0);
    Vel_x.resize(Nx,Ny,Nz);         Vel_x.fill(0);        // Gradient of the phase indicator field
    Vel_y.resize(Nx,Ny,Nz);         Vel_y.fill(0);
    Vel_z.resize(Nx,Ny,Nz);         Vel_z.fill(0);
    BlobID.resize(Nx,Ny,Nz);        BlobID.fill(0);
    
    GradPhiX.resize(Nx,Ny,Nz);  GradPhiX.fill(0);
    GradPhiY.resize(Nx,Ny,Nz);  GradPhiY.fill(0);
    GradPhiZ.resize(Nx,Ny,Nz);  GradPhiZ.fill(0);
    
    //.........................................
    // Allocate cube storage space
    CubeValues.resize(2,2,2);
    nw_tris.resize(3,20);
    ns_tris.resize(3,20);
    ws_tris.resize(3,20);
    nws_seg.resize(2,20);
    local_sol_tris.resize(3,18);
    nw_pts=DTMutableList<Point>(20);
    ns_pts=DTMutableList<Point>(20);
    ws_pts=DTMutableList<Point>(20);
    nws_pts=DTMutableList<Point>(20);
    local_nws_pts=DTMutableList<Point>(20);
    local_sol_pts=DTMutableList<Point>(20);
    tmp=DTMutableList<Point>(20);
    //.........................................
    Values.resize(20);
    DistanceValues.resize(20);
    KGwns_values.resize(20);
    KNwns_values.resize(20);
    InterfaceSpeed.resize(20);
    NormalVector.resize(60);
    //.........................................
    vn.resize(3);
    vw.resize(3);
    
    nw.resize(3);
    vnw.resize(3);
    vws.resize(3);
    vns.resize(3);
    vnws.resize(3);
    Gwn.resize(6);
    Gns.resize(6);
    Gws.resize(6);
    Gwns.resize(6);
    vn_global.resize(3);
    vw_global.resize(3);
    
    
    nw_global.resize(3);
    vnw_global.resize(3);
    vws_global.resize(3);
    vns_global.resize(3);
    vnws_global.resize(3);
    Gwn_global.resize(6);
    Gns_global.resize(6);
    Gws_global.resize(6);
    Gwns_global.resize(6);
    // Allocate subdomain variables
    
    Eta_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Eta_sub1.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trimdist_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    swn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    porosity_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    poreVol_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    awn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    tmp_tmp.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    ans_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    aws_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    lwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    wp_volume_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    nwp_volume_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    As_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    vol_w_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    vol_n_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    pn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    pw_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    efawns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    euler_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    eulerW_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    An_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    Jwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    JwnW_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KwnW_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jws_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kws_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    JwsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KwsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    JnsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KnsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    eulerS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    KNwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KGwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trawn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trJwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trRwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    wwndnw_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    wwnsdnwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jwnwwndnw_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    
    nw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    vn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    vw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    
    vws_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    vns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    vnw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    vnws_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    Gwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Gns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Gws_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Gwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    pwn_sub.resize(  { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Xn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    aw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    an_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );

    an_sub1.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    an_sub2.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    an_sub3.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    
    qDistances.resize(Nx*Ny*Nz*18); qDistances.fill(-9);
    
    // w^wn
    wwn.resize(3);
    wwn_global.resize(3);
    wwn_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOwnO.resize(3);
    nwOwnO_global.resize(3);
    nwOwnO_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOwnO2.resize(3);
    nwOwnO_global2.resize(3);
    nwOwnO_sub2.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOwsO.resize(3);
    nwOwsO_global.resize(3);
    nwOwsO_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );

    nwOwsO2.resize(3);
    nwOwsO_global2.resize(3);
    nwOwsO_sub2.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOnsO2.resize(3);
    nwOnsO_global2.resize(3);
    nwOnsO_sub2.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    // Create global files
    if (Dm->rank()==0){
        TIMELOG_global = fopen("timelog.tcat","a+");
        if (fseek(TIMELOG_global,0,SEEK_SET) == fseek(TIMELOG_global,0,SEEK_CUR)) {
            // If timelog is empty, write a short header to list the averages
            const char header[] = "time porosity "
            "s^w p^w* p^n* eps^{wn}* eps^{ns}* eps^{ws}* J_n^{wn}* K_n^{wn}* l^{wns}* c^{wns} kappa_N^{wns}* kappa_G^{wns}* "   // Scalar averages
            "v_x^w* v_y^w* v_z^w* v_x^n* v_y^n* v_z^n* v_x^{wn}* v_y^{wn}* v_z^{wn}* v_x^{wns}* v_y^{wns}* v_z^{wns}* "         // Velocity averages
            "G_{xx}^{wn}* G_{yy}^{wn}* G_{zz}^{wn}* G_{xy}^{wn}* G_{xz}^{wn}* G_{yz}^{wn}* "
            "G_{xx}^{ws}* G_{yy}^{ws}* G_{zz}^{ws}* G_{xy}^{ws}* G_{xz}^{ws}* G_{yz}^{ws}* "
            "G_{xx}^{ns}* G_{yy}^{ns}* G_{zz}^{ns}* G_{xy}^{ns}* G_{xz}^{ns}* G_{yz}^{ns}* "
            "p_n^{wn}* X^n* eps^w eps^n "
            "G_{xx}^{wns}* G_{yy}^{wns}* G_{zz}^{wns}* G_{xy}^{wns}* G_{xz}^{wns}* G_{yz}^{wns}* "
            "v_x^{ws}* v_y^{ws}* v_z^{ws}* v_x^{ns}* v_y^{ns}* v_z^{ns}* "
            "nw_x- nw_y- nw_z- Eta- wwnx- wwny- wwnz- Ca Re";
            fprintf(TIMELOG_global,"%s\n",header);
        }
        
        /*NWPLOG = fopen("components.NWP.tcat","a+");
        fprintf(NWPLOG,"time label vol p^n awn ans J_n^{wn} K_n^{wn} lwns c^{wns} kappa_N^{wns} kappa_G^{wns} ");
        fprintf(NWPLOG,"v_x^n v_y^n v_z^n v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} ");
        fprintf(NWPLOG,"G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} ");
        fprintf(NWPLOG,"G_{xx}^{ns} G_{yy}^{ns} G_{zz}^{ns} G_{xy}^{ns} G_{xz}^{ns} G_{yz}^{ns} ");
        fprintf(NWPLOG,"Euler Cx Cy Cz vsq \n");
        
        WPLOG = fopen("components.WP.tcat","a+");
        fprintf(WPLOG,"time label vol p^w awn aws J_n^{wn} K_n^{wn} lwns c^{wns} kappa_N^{wns} kappa_G^{wns} ");
        fprintf(WPLOG,"v_x^w v_y^w v_z^w v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} ");
        fprintf(WPLOG,"G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} ");
        fprintf(WPLOG,"G_{xx}^{ws} G_{yy}^{ws} G_{zz}^{ws} G_{xy}^{ws} G_{xz}^{ws} G_{yz}^{ws} ");
        fprintf(WPLOG,"Cx Cy Cz vsq \n");*/

        NWPLOG = fopen("components.NWP.tcat","a+");
        fprintf(NWPLOG,"time label vol \n");
        
        WPLOG = fopen("components.WP.tcat","a+");
        fprintf(WPLOG,"time label vol \n");        
    }
    
    // Create the local files
    
    if ( subdivide[0] * subdivide[1] * subdivide[2] > 1 ) {
        char tmp[64];
        sprintf(tmp,"%05d",Dm->rank());
        std::string rankStr( tmp );
        std::string filename = "timelog.tcat." + rankStr;
        TIMELOG_local = fopen(filename.c_str(),"a+");
        
        const char header[] = "time porosity "                  // Timestep, porosity, Change in Surface Energy
        "s^w p^w p^n eps^{wn} eps^{ns} eps^{ws} J_n^{wn} K_n^{wn} l^{wns} c^{wns} kappa_N^{wns} kappa_G^{wns} " // Scalar averages
        "v_x^w v_y^w v_z^w v_x^n v_y^n v_z^n v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} " // Velocity averages
        "G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} "
        "G_{xx}^{ws} G_{yy}^{ws} G_{zz}^{ws} G_{xy}^{ws} G_{xz}^{ws} G_{yz}^{ws} "
        "G_{xx}^{ns} G_{yy}^{ns} G_{zz}^{ns} G_{xy}^{ns} G_{xz}^{ns} G_{yz}^{ns} "
        "p_n^{wn} X^n eps^w eps^n "
        "G_{xx}^{wns} G_{yy}^{wns} G_{zz}^{wns} G_{xy}^{wns} G_{xz}^{wns} G_{yz}^{wns} "
        "v_x^{ws} v_y^{ws} v_z^{ws} v_x^{ns} v_y^{ns} v_z^{ns} "
        "nw_x nw_y nw_z Eta wwnx wwny wwnz Ca Re";
        
        fprintf(TIMELOG_local,"%s\n",header);
        //   SUB-DIVIDE SECTION
        // for ( auto file : TIMELOG_subdomain ) fprintf(file,"%s\n",header);
    }
    
    // Initialize variables
    Initialize();
}


// Destructor
TwoPhase::~TwoPhase()
{
    //delete [] TempID;
    if ( NWPLOG!=nullptr ) { fclose(NWPLOG); }
    if ( WPLOG!=nullptr ) { fclose(WPLOG); }
    if ( TIMELOG_local!=nullptr ) { fclose(TIMELOG_local); }

}


double TwoPhase::SignedVolumeOfTriangle(Point p1, Point p2, Point p3) {
    double v321 = p3.x * p2.y * p1.z;
    double v231 = p2.x * p3.y * p1.z;
    double v312 = p3.x * p1.y * p2.z;
    double v132 = p1.x * p3.y * p2.z;
    double v213 = p2.x * p1.y * p3.z;
    double v123 = p1.x * p2.y * p3.z;
    return (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);
}


void TwoPhase::SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha)
{
    double Fx = force_x;
    double Fy = force_y;
    double Fz = force_z;
    double rho_n = rhoA;
    double rho_w = rhoB;
    double nu_n = (tauA-0.5)/3.f;
    double nu_w = (tauB-0.5)/3.f;
    double gamma_wn = 5.796*alpha;
    
}


inline void MeanFilter2(DoubleArray &Temp, DoubleArray &Mesh, Array<int> & domain){}
inline void MeanFilter3(DoubleArray &Temp, DoubleArray &Mesh, Array<int> & domain){}


void TwoPhase::ColorToSignedDistance(double Beta, DoubleArray &ColorData, DoubleArray &DistData)
{
    Dm->CommunicateMeshHalo(ColorData);
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                DistData(i,j,k) = ColorData(i,j,k);
            }
        }
    }
}

void TwoPhase::SmoothSameSign(double Beta, DoubleArray &ColorData, DoubleArray &DistData){}
void TwoPhase::MedianFilter(double Beta, DoubleArray &ColorData, DoubleArray &DistData){}
void TwoPhase::ComputeDelPhi(){}

void TwoPhase::ComputeSDw() {
    
    int i,j,k,n;
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                n = k*Nx*Ny+j*Nx+i;
                if (!(Dm->id[n] > 0)){
                    // Solid phase
                    SDw(i,j,k) = -1;
                }
                else {
                    SDw(i,j,k) = -SDn(i,j,k);
                }
            }
        }
    }
}


void TwoPhase::PosAdjustSDs() {
    
    int i,j,k,n;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    for (int k=0; k<Nz; k++) {
        for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx; i++) {
                size_t ijk = i + j * strideY + k * strideZ;
                SDs(ijk) += 1.1;
            }
        }
    }
}


void TwoPhase::NegAdjustSDs() {
    
    int i,j,k,n;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    for (int k=0; k<Nz; k++) {
        for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx; i++) {
                size_t ijk = i + j * strideY + k * strideZ;
                SDs(ijk) -= 1.1;
            }
        }
    }
}


void TwoPhase::Initialize()
{
    
    Ca_global = 0;
    Re_global = 0;
    
    awn = 0;
    awn_global = 0;

    ans = 0;
    ans_global = 0;
    
    aws = 0;
    aws_global = 0;

    Jwn = Kwn = Jwn_global = Kwn_global = 0;
    local_nw_phase_volume = 0.0;
    global_nw_phase_volume = 0.0;
    mass_local = 0.0;
    mass_global = 0.0;
    porosity = 0.0;
    porosity_sub.fill(0);
    porosity_global = 0;

    aw_sub.fill(0);
    aw = 0;
    aw_global = 0;

    an_sub.fill(0);
    an = 0;
    an_global = 0;
    
    efawns_count = 0;

    vn.fill(0);
    vn_global.fill(0);
    
    vol_w_sub.fill( 0 );
    vol_n_sub.fill( 0 );
    vol_w = vol_n =0.0;
    vol_w_global = 0;
    vol_n_global = 0;
    
    pn_sub.fill( 0.0 );
    pw_sub.fill( 0.0 );
    pn = pw = pwn = 0.0;
    pn_global = pw_global = 0.0;
    
    Eta_sub.fill(0);
    Eta = Eta_global = 0;
    total_volume_global = 0;

    Eta_sub1.fill(0);
    Eta1 = Eta_global1 = 0;
    
    /*wsGlobal_dS = 0;
    nsGlobal_dS =0;
    nwGlobal_dS =0;
    
    wsGlobal_CCs =0;
    nsGlobal_CCs =0;
    nwGlobal_CCs =0;
    
    nEuler = wEuler = sEuler = 0.0;
    nwEuler = nsEuler = wsEuler = 0.0;
    nGlobalEuler = sGlobalEuler = wGlobalEuler = 0.0;*/
    
    trimdist=0.0;
    fluid_isovalue = 0.0;
    solid_isovalue= 0.0;
    // Initialize the averaged quantities
    porosity = 0.0;
    awn = aws = ans = lwns = 0.0;
    nwp_volume = wp_volume = 0.0;
    As = 0.0;
    pn = pw = pwn = 0.0;
    vw.fill( 0 );
    vn.fill( 0 );
    vnw.fill( 0 );
    vws.fill( 0 );
    vns.fill( 0 );
    vnws.fill( 0 );
    Gwn.fill( 0 );
    Gws.fill( 0 );
    Gns.fill( 0 );
    Gwns.fill( 0 );
    vol_w = vol_n =0.0;
    KGwns = KNwns = 0.0;
    Jwn = Kwn = efawns = 0.0;
    Jws = Kws = Jns = Kns = 0;
    trJwn = trawn = trRwn = 0.0;
    euler = Jn = An = Kn = eulerW = 0.0;
    wwndnw = 0.0; wwnsdnwn = 0.0; Jwnwwndnw=0.0;
    vnwdnw = 0.0;
    // Initialize subdomain variables
    trimdist_sub.fill( 0 );
    porosity_sub.fill( 0 );
    porosity_global = 0.0;
    poreVol_sub.fill( 0 );
    awn_sub.fill( 0 );
    tmp_tmp.fill( 0 );
    ans_sub.fill( 0 );
    aws_sub.fill( 0 );
    lwns_sub.fill( 0 );
    wp_volume_sub.fill( 0 );
    nwp_volume_sub.fill( 0 );
    As_sub.fill( 0 );
    vol_w_sub.fill( 0 );
    vol_n_sub.fill( 0 );
    pn_sub.fill( 0 );
    pw_sub.fill( 0 );
    pwn_sub.fill ( 0 );
    efawns_sub.fill( 0 );
    euler_sub.fill( 0 );
    euler = 0;
    euler_global = 0;
    eulerW_sub.fill( 0 );
    eulerW = 0;
    eulerW_global = 0;
    Kn_sub.fill( 0 );
    Jn_sub.fill( 0 );
    An_sub.fill( 0 );

    GaussCurvature.fill( 0 );
    GaussCurvatureW.fill( 0 );
    MeanCurvature.fill( 0 );
    MeanCurvatureW.fill( 0 );
    
    Jwn_sub.fill( 0 );
    Kwn_sub.fill( 0 );
    Jws_sub.fill( 0 );
    Kws_sub.fill( 0 );
    Jns_sub.fill( 0 );
    Kns_sub.fill( 0 );
    
    KwnW_sub.fill( 0 );
    JwnW_sub.fill( 0 );
    JwnW = 0;
    KwnW = 0;
    JwnW_global = 0;
    KwnW_global = 0;
    
    JwsS_sub.fill( 0 );
    KwsS_sub.fill( 0 );
    JnsS_sub.fill( 0 );
    KnsS_sub.fill( 0 );
    eulerS_sub.fill( 0 );
    eulerS = 0;
    eulerS_global = 0;
    
    KNwns_sub.fill( 0 );
    KGwns_sub.fill( 0 );
    trawn_sub.fill( 0 );
    trJwn_sub.fill( 0 );
    trRwn_sub.fill( 0 );
    wwndnw_sub.fill( 0 );
    wwnsdnwn_sub.fill( 0 );
    Jwnwwndnw_sub.fill( 0 );
    vn_sub.fill( 0 );
    vw_sub.fill( 0 );
    vws_sub.fill( 0 );
    vns_sub.fill( 0 );
    vnw_sub.fill( 0 );
    vnws_sub.fill( 0 );
    Gwn_sub.fill( 0 );
    Gns_sub.fill( 0 );
    Gws_sub.fill( 0 );
    Gwns_sub.fill( 0 );
    
    //    vnwdnw_sub.fill( 0 );
    //    vnwdnw = 0;
    //    vnwdnw_global = 0;
    
    wwn.fill( 0 );
    wwn_global.fill( 0 );
    wwn_sub.fill( 0 );
    
    nwOwnO_sub.fill( 0 );
    nwOwnO.fill( 0 );
    nwOwnO_global.fill( 0 );
    
    nwOwnO_sub2.fill( 0 );
    nwOwnO2.fill( 0 );
    nwOwnO_global2.fill( 0 );
    
    nwOwsO_sub.fill( 0 );
    nwOwsO.fill( 0 );
    nwOwsO_global.fill( 0 );

    nwOnsO_sub2.fill( 0 );
    nwOnsO2.fill( 0 );
    nwOnsO_global2.fill( 0 );
}


void TwoPhase::UpdateSolid()
{
    Dm->CommunicateMeshHalo(SDs);
    //...........................................................................
    // Gradient of the Signed Distance function
    //...........................................................................
    pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_z);
    //...........................................................................
}


/*double TwoPhase::FiveCases(double & phi1, double & phi2, double & phi3, double & phi4) {
    double temp;
    double numerator, denominator;
    
    if (phi4 <= 0.0) {
        return 1.0;
    }
    
    if (phi3 <= 0 && phi4 > 0) {
        numerator = phi4*phi4*phi4;
        denominator = (phi4-phi1)*(phi4-phi2)*(phi4-phi3);
        if (denominator == 0) return 1.0;
        else { return 1.0-numerator/denominator;}
    }
    
    if (phi2 <=0 && phi3 > 0) {
        numerator = phi1*phi2*(phi3*phi3 + phi3*phi4 + phi4*phi4)
        + phi3*phi4*( phi3*phi4 - (phi1 + phi2)*(phi3 + phi4) );
        denominator = (phi1-phi3)*(phi2-phi3)*(phi1-phi4)*(phi2-phi4);
        if (denominator == 0) return 0;
        else { return numerator/denominator;}
    }
    
    if (phi1 <=0 && phi2 > 0) {
        numerator = -1.0*phi1*phi1*phi1;
        denominator = (phi2-phi1)*(phi3-phi1)*(phi4-phi1);
        if (denominator == 0) return 0;
        else { return numerator/denominator;}
    }
    
    if (phi1 > 0.0) {
        return 0.0;
    }
    
}*/

/*oid TwoPhase::ComputeVolumeFraction(DoubleArray & A, DoubleArray & VFPhase) {
    
    // Cardinality
    double Cardinality = 24.0;
    double phi1,phi2,phi3,phi4;
    double fraction = 0;
    double x, X, y, Y, z, Z;
    double sum;
    
    // Get centroid value
    
    double C0, C1, C2, C3, C4, C5, C6, C7, C8;
    
    vector<double> vec;
    vec.resize(4);
    
    
    for (int k=1; k<Nz-1; k++)
        for (int j=1; j<Ny-1; j++)
            for (int i=1; i<Nx-1; i++) {
                sum = 0.0;
                
                // checked
                C0 = A(i,j,k);
                C1 = A(i+1,j,k);
                C2 = A(i,j+1,k);
                C3 = A(i,j,k+1);
                
                C4 = A(i+1,j+1,k);
                C5 = A(i+1,j,k+1);
                C6 = A(i,j+1,k+1);
                C7 = A(i+1,j+1,k+1);
                C8 = 0.125*( C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7 );
                
                // checked
                X = 0.25*( C1 + C4 + C7 + C5 );
                x = 0.25*( C0 + C2 + C6 + C3 );
                Y = 0.25*( C2 + C4 + C7 + C6 );
                y = 0.25*( C0 + C1 + C5 + C3 );
                Z = 0.25*( C3 + C5 + C7 + C6 );
                z = 0.25*( C0 + C1 + C4 + C2 );
                
                // BOTTOM
                // 1
                vec.at(0) = C0; vec.at(1) = C1; vec.at(2) = z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 2
                vec.at(0) = C1; vec.at(1) = C4; vec.at(2) = z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 3
                vec.at(0) = C4; vec.at(1) = C2; vec.at(2) = z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 4
                vec.at(0) = C2; vec.at(1) = C0; vec.at(2) = z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                
                
                // TOP
                // 5
                vec.at(0) = C3; vec.at(1) = C5; vec.at(2) = Z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 6
                vec.at(0) = C5; vec.at(1) = C7; vec.at(2) = Z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 7
                vec.at(0) = C7; vec.at(1) = C6; vec.at(2) = Z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 8
                vec.at(0) = C6; vec.at(1) = C3; vec.at(2) = Z; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                
                
                // LEFT
                // 9
                vec.at(0) = C0; vec.at(1) = C1; vec.at(2) = y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 10
                vec.at(0) = C1; vec.at(1) = C5; vec.at(2) = y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 11
                vec.at(0) = C5; vec.at(1) = C3; vec.at(2) = y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 12
                vec.at(0) = C3; vec.at(1) = C0; vec.at(2) = y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                
                
                // RIGHT
                // 13
                vec.at(0) = C2; vec.at(1) = C4; vec.at(2) = Y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 14
                vec.at(0) = C4; vec.at(1) = C7; vec.at(2) = Y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 15
                vec.at(0) = C7; vec.at(1) = C6; vec.at(2) = Y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 16
                vec.at(0) = C6; vec.at(1) = C2; vec.at(2) = Y; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                
                
                // BACK
                // 17
                vec.at(0) = C0; vec.at(1) = C2; vec.at(2) = x; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 18
                vec.at(0) = C2; vec.at(1) = C6; vec.at(2) = x; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 19
                vec.at(0) = C6; vec.at(1) = C3; vec.at(2) = x; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 20
                vec.at(0) = C3; vec.at(1) = C0; vec.at(2) = x; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                
                
                // FRONT
                // 21
                vec.at(0) = C1; vec.at(1) = C4; vec.at(2) = X; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 22
                vec.at(0) = C4; vec.at(1) = C7; vec.at(2) = X; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 23
                vec.at(0) = C7; vec.at(1) = C5; vec.at(2) = X; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                // 24
                vec.at(0) = C5; vec.at(1) = C1; vec.at(2) = X; vec.at(3) = C8;
                sort(vec.begin(), vec.end()); phi1 = vec.at(0); phi2 = vec.at(1); phi3 = vec.at(2); phi4 = vec.at(3);
                fraction = FiveCases(phi1, phi2, phi3, phi4);  sum += fraction;
                
                
                
                // RETURN
                VFPhase(i,j,k) = sum / Cardinality;
            }
    
    Dm->CommunicateMeshHalo(VFPhase);
}*/


void TwoPhase::UpdateMeshValues()
{
    int i,j,k,n;
    double fx,fy,fz;

    // Make sure that SDn is communicated before use
    for (k=1; k<Nz-1; k++)
    for (j=1; j<Ny-1; j++)
    for (i=1; i<Nx-1; i++){
        fx = 0.5*(SDn(i+1,j,k) - SDn(i-1,j,k));
        fy = 0.5*(SDn(i,j+1,k) - SDn(i,j-1,k));
        fz = 0.5*(SDn(i,j,k+1) - SDn(i,j,k-1));
        DelPhi(i,j,k) = sqrt(fx*fx+fy*fy+fz*fz);
    }
    Dm->CommunicateMeshHalo(DelPhi);
    
    
    Dm->CommunicateMeshHalo(SDn_tplus);
    for (k=1; k<Nz-1; k++){
        for (j=1; j<Ny-1; j++){
            for (i=1; i<Nx-1; i++){
                n = k*Nx*Ny+j*Nx+i;
                dPdt(n) = SDn(n)-SDn_tplus(n);
            }
        }
    }
    Dm->CommunicateMeshHalo(dPdt);
    //...........................................................................
    pmmc_MeshGradient(SDn,SDn_x,SDn_y,SDn_z,Nx,Ny,Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDn_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDn_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDn_z);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_z);
    //...........................................................................
    Dm->CommunicateMeshHalo(Press);


    
    Dm->CommunicateMeshHalo(SDn);
    
    Dm->CommunicateMeshHalo(SDn_tplus);
    
    pmmc_MeshCurvature(SDn, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(MeanCurvature);
    //...........................................................................
    Dm->CommunicateMeshHalo(GaussCurvature);

    for (k=1; k<Nz-1; k++){
        for (j=1; j<Ny-1; j++){
            for (i=1; i<Nx-1; i++){
                n = k*Nx*Ny+j*Nx+i;
                SDw(i,j,k)=-SDn(i,j,k);
                SDw_tplus(i,j,k) =-SDn_tplus(i,j,k);
            }
        }
    }

    Dm->CommunicateMeshHalo(SDw_tplus);
    
    Dm->CommunicateMeshHalo(SDw);
    //...........................................................................
    pmmc_MeshGradient(SDw,SDw_x,SDw_y,SDw_z,Nx,Ny,Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDw_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDw_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDw_z);
    //...........................................................................
    
    for (int n=0; n<Nx*Ny*Nz; n++) { double magSDw = std::sqrt(SDw_x(n)*SDw_x(n) + SDw_y(n)*SDw_y(n) + SDw_z(n)*SDw_z(n));
        
        if (magSDw > 0.0) {normalw_x(n) = SDw_x(n)/magSDw; }
        else { normalw_x(n) = 0.0; }
        if (magSDw > 0.0) {normalw_y(n) = SDw_y(n)/magSDw; }
        else { normalw_y(n) = 0.0; }
        if (magSDw > 0.0) { normalw_z(n) = SDw_z(n)/magSDw; }
        else { normalw_z(n) = 0.0; }

    }

    //...........................................................................

    Dm->CommunicateMeshHalo(SDs);
    //...........................................................................
    // Gradient of the Signed Distance function
    //...........................................................................
    pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_z);

    //Initializing the blob ID
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                n = k*Nx*Ny+j*Nx+i;
                if (!(Dm->id[n] > 0)){
                    // Solid phase
                    PhaseID(i,j,k) = 0;
                }
                else if (SDn(i,j,k) < 0.0){
                    // wetting phase
                    PhaseID(i,j,k) = 2;
                }
                else {
                    // non-wetting phase
                    PhaseID(i,j,k) = 1;
                }
            }
        }
    }
    
    
}


//inline double normalEqnsExtrapolation(double y0, double y1, double dist) {}


void TwoPhase::FindIsolatedLatticeSites() {}


void TwoPhase::ComputePhiDistance(DoubleArray &Field) {
    double SqrtMagGradPhase = 0;
    int ijk;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    
    Array<char> id_solid(Nx,Ny,Nz);

    double factor,temp,value;
    factor=0.5/0.9;
    // Initialize to -1,1 (segmentation)
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                int n = k*Nx*Ny+j*Nx+i;
                double value = Field(n);
                if (value > 0)    id_solid(n) = 1;
                else              id_solid(n) = 0;
                temp = factor*log((1.0+value)/(1.0-value));
                if (value > 0.8) Field(n) = 2.94*factor;
                else if (value < -0.8) Field(n) = -2.94*factor;
                else Field(n) = temp;
            }
        }
    }
   
    CalcDist(Field,id_solid,*Dm);
    
}


void TwoPhase::InactiveCellMeanFilter(DoubleArray &Field) {

    DoubleArray SavedField; SavedField.resize(Nx,Ny,Nz); SavedField.fill(0);

    double sf;
    int ijk;
    int strideY = Nx;
    int strideZ = Nx*Ny;

    for (int k=0; k < Nz; k++) {
        for (int j=0; j < Ny; j++) {
            for (int i=0; i < Nx; i++) {
                ijk = i + j * strideY + k * strideZ;
                SavedField(ijk) = Field(ijk);
            }
        }
    }

    int kmin=1; int kmax=Nz-1;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=1;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=Nz-1; // for BCC test - should normally be -4

    int nn;
    for (int k=kmin; k < kmax; k++) {
        for (int j=1; j < Ny-1; j++) {
            for (int i=1; i < Nx-1; i++) {

                ijk = i + j * strideY + k * strideZ;

                //Field(ijk) = Field(ijk) * (1./3.);

                if (VFmask(ijk) > 0.5 && VFmask(ijk) < 1.0) {

                    nn = ijk-1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideY;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideY;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideZ;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideZ;
                    Field(ijk) += (SavedField(nn) );

                    nn = ijk-strideY-1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideY+1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideY-1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideY+1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideZ-1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideZ+1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideZ-1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideZ+1;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideZ-strideY;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideZ+strideY;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk+strideZ-strideY;
                    Field(ijk) += (SavedField(nn) );
                    
                    nn = ijk-strideZ+strideY;
                    Field(ijk) += (SavedField(nn) );
                    
                    Field(ijk) /= 18.0;
                }
            }
        }
    }
    Dm->CommunicateMeshHalo(Field);
}

void TwoPhase::MeanFilter(DoubleArray &Field) {
    
    DoubleArray SavedField; SavedField.resize(Nx,Ny,Nz); SavedField.fill(0);

    double sf;
    int ijk;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    for (int k=0; k < Nz; k++) {
        for (int j=0; j < Ny; j++) {
            for (int i=0; i < Nx; i++) {
                ijk = i + j * strideY + k * strideZ;
                SavedField(ijk) = Field(ijk);
            }
        }
    }

    int kmin=1; int kmax=Nz-1;

    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=1;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=Nz-1; // for BCC test - should normally be -4

    int nn;

    for (int k=kmin; k < kmax; k++) {
        for (int j=1; j < Ny-1; j++) {
            for (int i=1; i < Nx-1; i++) {

                ijk = i + j * strideY + k * strideZ;

                Field(ijk) = Field(ijk) * (1./3.);

                nn = ijk-1;
                Field(ijk) += (SavedField(nn) * (1./18.));

                nn = ijk+1;
                Field(ijk) += (SavedField(nn) * (1./18.));

                nn = ijk-strideY;
                Field(ijk) += (SavedField(nn) * (1./18.));

                nn = ijk+strideY;
                Field(ijk) += (SavedField(nn) * (1./18.));

                nn = ijk-strideZ;
                Field(ijk) += (SavedField(nn) * (1./18.));

                nn = ijk+strideZ;
                Field(ijk) += (SavedField(nn) * (1./18.));

                nn = ijk-strideY-1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk+strideY+1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk+strideY-1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk-strideY+1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk-strideZ-1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk+strideZ+1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk+strideZ-1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk-strideZ+1;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk-strideZ-strideY;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk+strideZ+strideY;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk+strideZ-strideY;
                Field(ijk) += (SavedField(nn) * (1./36.));

                nn = ijk-strideZ+strideY;
                Field(ijk) += (SavedField(nn) * (1./36.));
            }
        }
    }
    Dm->CommunicateMeshHalo(Field);
}


/*void TwoPhase::ToBinary(DoubleArray & Field) {
    int ijk;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    for (int k=0; k < Nz; k++) {
        for (int j=0; j < Ny; j++) {
            for (int i=0; i < Nx; i++) {
                ijk = i + j * strideY + k * strideZ;
                if ( Field(ijk) == affinity) Field(ijk) = -1;
                if ( Field(ijk) > 0) Field(ijk) = 1;
                if (Field(ijk) <= 0) Field(ijk) = -1;
            }
        }
    }
}*/


void TwoPhase::ZeroSolidPhaseLatticeSites(DoubleArray & Field) {

    int count;
    double vfs;

    int strideY = Nx;
    int strideZ = Nx*Ny;

    int n,nn,ijk,m,mm;
    double denA,denB;
    double m0;
    int id;

    size_t N = Nx*Ny*Nz;

    for (int k=1; k < Nz-1; k++) {
        for (int j=1; j < Ny-1; j++) {
            for (int i=1; i < Nx-1; i++) {
                ijk = i + j * strideY + k * strideZ;
                if (VFmask(ijk) == 1) {
                    if (Field(ijk) == affinity) {
                        Field(ijk) = affinity;
                    }
                }
            }
        }
    }
    Dm->CommunicateMeshHalo(Field);
}

/*void TwoPhase::SavePhaseIntoOldPhase(DoubleArray & NewField, DoubleArray & Field) {
    
    int ijk;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    for (int k=1; k < Nz-1; k++) {
        for (int j=1; j < Ny-1; j++) {
            for (int i=1; i < Nx-1; i++) {
                ijk = i + j * strideY + k * strideZ;
                NewField(ijk) = Field(ijk);
            }
        }
    }
    
}*/


// Done before analysis is performed to compute VFs 
void TwoPhase::PhaseInit(double affinity)
{
    int i,j,k;
    size_t n;
    
    int N = Nx*Ny*Nz;
    
    visMyMesh vis_ws;
    visMyMesh vis_ns;
    visMyMesh vis_nw;
    
    int ncubes=(Nx-1)*(Ny-1)*(Nz-1);
    int npx = Dm->nprocx();
    int npy = Dm->nprocy();
    int npz = Dm->nprocz();
    
    Point Q,D,E,F;
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    
    // If external boundary conditions are set, do not average over the inlet
    int kmin=1; int kmax=Nz-1;
        if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=4;
        if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=Nz-4; // for BCC test - should normally be -4
    
    int proc_dat = (Dm->nprocz()-1);
    
    int x,y,z;
    double ws_vol, ns_vol, nw_vol;
    int CHerror = 0;
    
    double area_x, area_X, area_y, area_Y, area_z, area_Z;
    double nx,ny,nz;
    
    DoubleArray nwOwnO_tmp(3);
    DoubleArray nwOwsO_tmp(3);
    DoubleArray nw_tmp(3);
    DoubleArray wwn_tmp(3);
    DoubleArray vnw_tmp(3), vnws_tmp(3), vws_tmp(3), vns_tmp(3), Gwn_tmp(6), Gws_tmp(6), Gns_tmp(6), Gwns_tmp(6);
    //size_t N = Nx*Ny*Nz;
    int ii,jj,kk;
    for (int ll = 0; ll < NSpheresOnRank; ll++) {
        for (int h = 0; h < N; h++) {
            kk = h/(Nx*Ny); jj = (h-Nx*Ny*kk)/Nx; ii = h-Nx*Ny*kk-Nx*jj;
            ParticularSDs(ii,jj,kk) = rankSDs(h + ll*N);
        }
        for (int k = kmin; k < kmax; k++) {
            for (int j=1; j<Ny-1; j++) {
                for (int i=1; i<Nx-1; i++) {
                    n = i + j*Nx + k*Nx*Ny;
                    /* Volume portion of the analysis - replaces coarse Marching Cubes: Interial calculation is slow */
                    n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
                    n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
                    //...........................................................................
                    pmmc_ConstructLocalCube(ParticularSDs, zero, solid_isovalue, fluid_isovalue,
                                            nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                            local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                            n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                            n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg, i, j, k, Nx, Ny, Nz);
                    if (n_ws_tris > 0) {
                        /* Allocate new mesh and vertex iterator */
                        wsMyMesh m_ws;
                        wsMyMesh::VertexIterator vi;
                        /* Create convex hull mesh object */
                        chMyMesh CH_ws;
                        
                        /* Add points to the mesh */
                        for (int r=0;r<n_ws_tris;r++) {
                            D = ws_pts(ws_tris(0,r)); E = ws_pts(ws_tris(1,r)); F = ws_pts(ws_tris(2,r));
                            
                            D.x += 1.0*Dm->iproc()*(Nx-2); E.x += 1.0*Dm->iproc()*(Nx-2); F.x += 1.0*Dm->iproc()*(Nx-2);
                            D.y += 1.0*Dm->jproc()*(Nx-2); E.y += 1.0*Dm->jproc()*(Nx-2); F.y += 1.0*Dm->jproc()*(Nx-2);
                            D.z += 1.0*Dm->kproc()*(Nx-2); E.z += 1.0*Dm->kproc()*(Nx-2); F.z += 1.0*Dm->kproc()*(Nx-2);
                            
                            vi = vcg::tri::Allocator<wsMyMesh>::AddVertices(m_ws,3);
                            vi->P()=wsMyMesh::CoordType ( D.x, D.y, D.z); ++vi;
                            vi->P()=wsMyMesh::CoordType ( E.x, E.y, E.z); ++vi;
                            vi->P()=wsMyMesh::CoordType ( F.x, F.y, F.z); ++vi;
                        }
                        
                        /* Add corner points to surface for volume calculations */
                        for (int p = 0; p < 8; p++) {
                            int m = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                            if ( ParticularSDs(m) < 0 ) {
                                z = m/(Nx*Ny); y = (m-Nx*Ny*z)/Nx; x = m-Nx*Ny*z-Nx*y;
                                Q.z = z; Q.z += 1.0*Dm->kproc()*(Nx-2);
                                Q.y = y; Q.y += 1.0*Dm->jproc()*(Nx-2);
                                Q.x = x; Q.x += 1.0*Dm->iproc()*(Nx-2);
                                
                                vi = vcg::tri::Allocator<wsMyMesh>::AddVertices(m_ws,1);
                                vi->P()=wsMyMesh::CoordType ( Q.x, Q.y, Q.z); vi++;
                            }
                        }
                        
                        int beforeVN, beforeFN, wsAfterVN, wsAfterFN;
                        beforeVN = wsAfterVN = beforeFN = wsAfterFN = 0;
                        beforeVN = m_ws.VN(); beforeFN = m_ws.FN();
                        vcg::tri::Clean<wsMyMesh>::RemoveDuplicateFace(m_ws);
                        vcg::tri::Clean<wsMyMesh>::RemoveDuplicateVertex(m_ws,false);
                        
                        /* Generate convex hull mesh */
                        vcg::tri::ConvexHull<wsMyMesh,chMyMesh>::ComputeConvexHull(m_ws,CH_ws);
                        
                        /* For visualization purposes */
                        vcg::tri::Append<visMyMesh, chMyMesh>::Mesh(vis_ws, CH_ws, false /* selected */, false /* adjFlag */);
                        
                        /* Compute volume fraction of phase */
                        vcg::tri::Inertia<chMyMesh> CH_inertia(CH_ws);
                        double CH_vol = CH_inertia.Mass();
                        
                        vcg::tri::UpdateNormal<chMyMesh>::PerFaceNormalized(CH_ws);
                        
                        /* Overwrite VFmask with volume fraction of solid phase */
                        VFmask(i,j,k) = CH_vol;
                        
                        vcg::tri::UpdateNormal<chMyMesh>::PerFaceNormalized(CH_ws);
                        
                        if (CH_vol >= 0.5 && CH_vol < 1) {
                            /* For visualization purposes */
                            vcg::tri::Append<visMyMesh, chMyMesh>::Mesh(vis_ws, CH_ws, false /* selected */, false /* adjFlag */);
                        }
                    }
                } // i
            } // j
        } // k
    }
    
    
   for (int k=1; k<Nz-1; k++) { // kmin kmax
          for (int j=1; j<Ny-1; j++) {
              for (int i=1; i<Nx-1; i++) {
                  n = i + j*Nx + k*Nx*Ny;
                  if (Dm->id[n] == 0 && VFmask(n) == 0)
                  {
                      VFmask(n) = 1.0;
                  }
              }
          }
      }
      
    Dm->CommunicateMeshHalo(VFmask);

    char LocalRankString[8];
    char LocalExtension[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s%s","s_inactive_vis_data.",LocalRankString,".off");
    vcg::tri::io::ExporterOFF<visMyMesh>::Save(vis_ws,LocalRankFilename, 1);
    //        if (Dm->rank() == 0 ) { vcg::tri::io::ExporterOFF<visMyMesh>::Save(vis_ws,"s_vis_ws_data0.off", 1); }
    //  if (Dm->rank() == 0 ) { vcg::tri::io::ExporterOFF<visMyMesh>::Save(vis_ns,"s_vis_ns_data0.off", 1); }
}


/*void TwoPhase::PrintPhaseI(DoubleArray & field) {
    
    ofstream file;
    file.open ("phase_field_sample.txt");
    for (int k = 1; k < Nz-1; k++)
        for (int j = 1; j < Ny-1; j++)
            for (int i = 1; i < Nx-1; i++) {
                file << field(i,j,k) << " ";
            }
    file.close();
    
}*/


/*void TwoPhase::PrintPhaseII(DoubleArray & field) {
    
    ofstream file;
    file.open ("VF_field.txt");
    for (int k = 1; k < Nz-1; k++)
        for (int j = 1; j < Ny-1; j++)
            for (int i = 1; i < Nx-1; i++) {
                file << field(i,j,k) << " ";
            }
    file.close();
    
}*/


//void TwoPhase::NonWettingPhaseVolume() {}


void TwoPhase::FluidCentroidInit(double affinity){}


void TwoPhase::InitializeRankSignDistance()
{
    size_t N = Nx*Ny*Nz;
    
    // ! parallel plates
    rankSDs.resize(N*NSpheresOnRank); rankSDs.fill(0);
    ParticularSDs.resize(Nx,Ny,Nz); ParticularSDs.fill(0);
    
    ParticularSDs_x.resize(Nx,Ny,Nz); ParticularSDs_x.fill(0);
    ParticularSDs_y.resize(Nx,Ny,Nz); ParticularSDs_y.fill(0);
    ParticularSDs_z.resize(Nx,Ny,Nz); ParticularSDs_z.fill(0);
    
    SDn_tplus.resize(Nx,Ny,Nz);           SDn_tplus.fill(0);
    SDn_tplusplus.resize(Nx,Ny,Nz);       SDn_tplusplus.fill(0);

}


void TwoPhase::RankSignDistance() {}
void TwoPhase::CuttingPlaneZ(int zval) {}

double MySqr(double a_fVal) {  return a_fVal*a_fVal; }

void TwoPhase::ComputeLocal()
{
    
    int i,j,k;
    
    nwMyMesh m_nw;
    nsMyMesh m_ns;
    wsMyMesh m_ws;
    
    nMyMesh m_n;
    sMyMesh m_s;
    wMyMesh m_w;
    
    int strideY = Nx;
    int strideZ = Nx*Ny;
    
    //size_t fullncubes=npx*npy*npz*(Nx-1)*(Ny-1)*(Nz-1);
    int ncubes=(Nx-1)*(Ny-1)*(Nz-1);
    int npx = Dm->nprocx();
    int npy = Dm->nprocy();
    int npz = Dm->nprocz();

    Point P,A,B,C;
    Point Q,D,E,F;
    Point R,G,H,I;
    // Point S,J,K,L;
    
    int total_ranks = Dm->nprocz() * Dm->nprocy() * Dm->nprocx();
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    
    // If external boundary conditions are set, do not average over the inlet
    int kmin=1; int kmax=Nz-1;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=7;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-7;
    
    double dummy = 0;
    const int N2x = ( Nx-2 ) / subdivide[0];
    const int N2y = ( Ny-2 ) / subdivide[1];
    const int N2z = ( kmax - kmin ) / subdivide[2];
    
    subNodeVolume = 1.0*(N2x*N2y*N2z);
    
    size_t wsFacesCount = 0;
    size_t nfaces,nvertices;
    nfaces = nvertices = 0;
    
    size_t nsFacesCount = 0;
    size_t nsnfaces,nsnvertices;
    nsnfaces = nsnvertices = 0;
    
    double total_mass = 0;

    double WPHASEWALL = 1E-16;
    
    double minJwn = 1E16;
    double maxJwn = -1E16;
    DoubleArray nwOwnO_tmp(3);
    DoubleArray nwOwsO_tmp(3);
    DoubleArray nw_tmp(3);
    DoubleArray wwn_tmp(3);
    DoubleArray vnw_tmp(3), vnws_tmp(3), vws_tmp(3), vns_tmp(3), Gwn_tmp(6), Gws_tmp(6), Gns_tmp(6), Gwns_tmp(6);
    
    size_t N = Nx*Ny*Nz;
    int ii,jj,kk;

    for (int ll = 0; ll < NSpheresOnRank; ll++) {
        for (int h = 0; h < N; h++) {
            kk = h/(Nx*Ny); jj = (h-Nx*Ny*kk)/Nx; ii = h-Nx*Ny*kk-Nx*jj;
            ParticularSDs(ii,jj,kk) = rankSDs(h + ll*N);
        }
        
        for (k=1; k<Nz-1; k++){
            for (j=1; j<Ny-1; j++){
                for (i=1; i<Nx-1; i++){
                    // Compute all of the derivatives using finite differences
                    ParticularSDs_x(i,j,k) = 0.5*(ParticularSDs(i+1,j,k) - ParticularSDs(i-1,j,k));
                    ParticularSDs_y(i,j,k) = 0.5*(ParticularSDs(i,j+1,k) - ParticularSDs(i,j-1,k));
                    ParticularSDs_z(i,j,k) = 0.5*(ParticularSDs(i,j,k+1) - ParticularSDs(i,j,k-1));
                }
            }
        }
            
        for (int k=kmin; k<kmax; k++) { 
            const int k2 = std::min<int>( ( k - kmin) / N2z, subdivide[2] - 1 );  // kmin
            for (int j=1; j<Ny-1; j++){
                const int j2 = std::min<int>( ( j - 1 ) / N2y, subdivide[1] - 1 );
                for (int i=1; i<Nx-1; i++) {
                    const int i2 = std::min<int>( ( i - 1 ) / N2x, subdivide[0] - 1 );
                    
                    vnw_tmp.fill( 0 );  vws_tmp.fill( 0 );  vns_tmp.fill( 0 );  vnws_tmp.fill( 0 );
                    Gwn_tmp.fill( 0 );  Gws_tmp.fill( 0 );  Gns_tmp.fill( 0 );  Gwns_tmp.fill( 0 );
                    wwn_tmp.fill( 0 ); nwOwnO_tmp.fill( 0 ); nwOwsO_tmp.fill( 0 );
                    
                    // Compute volume averages
                    int n = i + j*Nx + k*Nx*Ny;
                    
                    //...........................................................................
                    n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
                    n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
                    //...........................................................................
                    pmmc_ConstructLocalCube(ParticularSDs, SDw, solid_isovalue, fluid_isovalue,
                                            nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                            local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                            n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                            n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                                            i, j, k, Nx, Ny, Nz);
                    //.......................................................................................
                    
                    if (n_ws_tris > 0) {
                        for (int r=0;r<n_ws_tris;r++) {
                            A = ws_pts(ws_tris(0,r));
                            B = ws_pts(ws_tris(1,r));
                            C = ws_pts(ws_tris(2,r));

                            double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                            double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                            double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                            double normal_x = ParticularSDs_x(i,j,k);
                            double normal_y = ParticularSDs_y(i,j,k);
                            double normal_z = ParticularSDs_z(i,j,k);

                            if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                P = A;
                                A = C;
                                C = P;
                            }

                            // Remap the points
                            A.x += 1.0*Dm->iproc()*(Nx-2);
                            A.y += 1.0*Dm->jproc()*(Nx-2);
                            A.z += 1.0*Dm->kproc()*(Nx-2);
                            B.x += 1.0*Dm->iproc()*(Nx-2);
                            B.y += 1.0*Dm->jproc()*(Nx-2);
                            B.z += 1.0*Dm->kproc()*(Nx-2);
                            C.x += 1.0*Dm->iproc()*(Nx-2);
                            C.y += 1.0*Dm->jproc()*(Nx-2);
                            C.z += 1.0*Dm->kproc()*(Nx-2);
                            wsMyMesh::FaceIterator fi = vcg::tri::Allocator<wsMyMesh>::AddFace(m_ws,wsMyMesh::CoordType ( A.x, A.y, A.z),wsMyMesh::CoordType ( B.x, B.y, B.z),wsMyMesh::CoordType ( C.x, C.y, C.z));
                        }
                    }
                    if (n_ns_tris > 0 ) {
                        for (int r=0;r<n_ns_tris;r++) {
                            A = ns_pts(ns_tris(0,r));
                            B = ns_pts(ns_tris(1,r));
                            C = ns_pts(ns_tris(2,r));

                            double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                            double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                            double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                            double normal_x = ParticularSDs_x(i,j,k);
                            double normal_y = ParticularSDs_y(i,j,k);
                            double normal_z = ParticularSDs_z(i,j,k);

                            if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                P = A;
                                A = C;
                                C = P;
                            }

                            // Remap the points
                            A.x += 1.0*Dm->iproc()*(Nx-2);
                            A.y += 1.0*Dm->jproc()*(Nx-2);
                            A.z += 1.0*Dm->kproc()*(Nx-2);
                            B.x += 1.0*Dm->iproc()*(Nx-2);
                            B.y += 1.0*Dm->jproc()*(Nx-2);
                            B.z += 1.0*Dm->kproc()*(Nx-2);
                            C.x += 1.0*Dm->iproc()*(Nx-2);
                            C.y += 1.0*Dm->jproc()*(Nx-2);
                            C.z += 1.0*Dm->kproc()*(Nx-2);
                            nsMyMesh::FaceIterator fi = vcg::tri::Allocator<nsMyMesh>::AddFace(m_ns,nsMyMesh::CoordType ( A.x, A.y, A.z),nsMyMesh::CoordType ( B.x, B.y, B.z),nsMyMesh::CoordType ( C.x, C.y, C.z));

                        }
                    }
                    if (n_nw_tris > 0) {
                        for (int r=0;r<n_nw_tris;r++) {
                            A = nw_pts(nw_tris(0,r));
                            B = nw_pts(nw_tris(1,r));
                            C = nw_pts(nw_tris(2,r));
                            
                            double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                            double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                            double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                            
                            double normal_x = -SDn_x(i,j,k);
                            double normal_y = -SDn_y(i,j,k);
                            double normal_z = -SDn_z(i,j,k);
                            
                            if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                P = A;
                                A = C;
                                C = P;
                            }
                            
                            // Remap the points
                            A.x += 1.0*Dm->iproc()*(Nx-2);
                            A.y += 1.0*Dm->jproc()*(Nx-2);
                            A.z += 1.0*Dm->kproc()*(Nx-2);
                            B.x += 1.0*Dm->iproc()*(Nx-2);
                            B.y += 1.0*Dm->jproc()*(Nx-2);
                            B.z += 1.0*Dm->kproc()*(Nx-2);
                            C.x += 1.0*Dm->iproc()*(Nx-2);
                            C.y += 1.0*Dm->jproc()*(Nx-2);
                            C.z += 1.0*Dm->kproc()*(Nx-2);
                            nwMyMesh::FaceIterator fi = vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                        }
                    }                    
                    
                    // wns common curve averages
                    if ( n_local_nws_pts > 0 ) {
                        wwnsdnwn_sub(i2,j2,k2) += pmmc_CommonCurveSpeed(CubeValues, dPdt, vnws_tmp,
                                                                        SDn_x, SDn_y, SDn_z,
                                                                        ParticularSDs_x,ParticularSDs_y,ParticularSDs_z,
                                                                        local_nws_pts,i,j,k,n_local_nws_pts);
                        vnws_sub(i2,j2,k2,0) += vnws_tmp(0);
                        vnws_sub(i2,j2,k2,1) += vnws_tmp(1);
                        vnws_sub(i2,j2,k2,2) += vnws_tmp(2);

                        pmmc_CurveCurvature(SDn, ParticularSDs,
                                            SDn_x, SDn_y, SDn_z,
                                            ParticularSDs_x,ParticularSDs_y,ParticularSDs_z,
                                            KNwns_values, KGwns_values, KNwns, KGwns,
                                            nws_pts, n_nws_pts, i, j, k);
                        KNwns_sub(i2,j2,k2) += KNwns;
                        KGwns_sub(i2,j2,k2) += KGwns;
                        lwns_sub(i2,j2,k2) +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
                        pmmc_CurveOrientation(Gwns_tmp,local_nws_pts,n_local_nws_pts,i,j,k);
                        for (int m=0; m<6; m++)  Gwns_sub(i2,j2,k2,m) += Gwns_tmp(m);
                    }

                    if (n_ws_tris > 0) {
                        aws_sub(i2,j2,k2) += pmmc_CubeSurfaceOrientation(Gws_tmp,ws_pts,ws_tris,n_ws_tris);
                        for (int m=0; m<6; m++) { Gws_sub(i2,j2,k2,m) += Gws_tmp(m); }
                        //aws_sub(i2,j2,k2) += pmmc_CubeSurfaceNormalVector(nwOwsO_tmp,ws_pts,ws_tris,n_ws_tris);
                        //for (int m=0; m<3; m++)  nwOwsO_sub(i2,j2,k2,m) += nwOwsO_tmp(m);
                        pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, ws_pts, ws_tris, NormalVector, InterfaceSpeed, vws_tmp, i, j, k, n_ws_pts, n_ws_tris);
                        
                        vws_sub(i2,j2,k2,0) += vws_tmp(0);
                        vws_sub(i2,j2,k2,1) += vws_tmp(1);
                        vws_sub(i2,j2,k2,2) += vws_tmp(2);

                    }
                    if (n_ns_tris > 0) {
                        ans_sub(i2,j2,k2) += pmmc_CubeSurfaceOrientation(Gns_tmp,ns_pts,ns_tris,n_ns_tris);
                        for (int m=0; m<6; m++) {  Gns_sub(i2,j2,k2,m) += Gns_tmp(m); }
                        pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z,  CubeValues, ns_pts,ns_tris, NormalVector, InterfaceSpeed, vns_tmp, i, j, k, n_ns_pts, n_ns_tris);
                        vns_sub(i2,j2,k2,0) += vns_tmp(0);
                        vns_sub(i2,j2,k2,1) += vns_tmp(1);
                        vns_sub(i2,j2,k2,2) += vns_tmp(2);

                        euler_sub(i2,j2,k2) += geomavg_EulerCharacteristic(ns_pts,ns_tris,n_ns_pts,n_ns_tris,i,j,k);
                    }
                    if (n_nw_tris > 0) {
                        
                        euler_sub(i2,j2,k2) += geomavg_EulerCharacteristic(nw_pts,nw_tris,n_nw_pts,n_nw_tris,i,j,k);
                        
                        pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris, NormalVector, InterfaceSpeed, vnw_tmp, i, j, k, n_nw_pts, n_nw_tris);
                        
                        vnw_sub(i2,j2,k2,0) += vnw_tmp(0);
                        vnw_sub(i2,j2,k2,1) += vnw_tmp(1);
                        vnw_sub(i2,j2,k2,2) += vnw_tmp(2);
                    }
                    if (n_nws_seg > 0) {
                        for (int r = 0; r < n_nws_seg; r++) {
                            nwsMyMesh m_nws;
                            Q = nws_pts(nws_seg(0,r));
                            R = nws_pts(nws_seg(1,r));
                            if (n_ns_tris > 0 ) {
                                for (int r=0;r<n_ns_tris;r++) {
                                    A = ns_pts(ns_tris(0,r));
                                    B = ns_pts(ns_tris(1,r));
                                    C = ns_pts(ns_tris(2,r));
                                    
                                    double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                                    double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                                    double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                                    
                                    double normal_x = ParticularSDs_x(i,j,k);
                                    double normal_y = ParticularSDs_y(i,j,k);
                                    double normal_z = ParticularSDs_z(i,j,k);
                                    
                                    if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                        P = A;
                                        A = C;
                                        C = P;
                                    }
                                    
                                    if (Q == A) {
                                        if (R == B) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                        if (R == C) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                    }
                                    if (Q == B) {
                                        if (R == A) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                        if (R == C) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                    }
                                    if (Q == C) {
                                        if (R == B) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                        if (R == A) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                    }
                                }
                            }
                            if (n_nw_tris > 0 ) {
                                for (int r=0;r<n_nw_tris;r++) {
                                    A = nw_pts(nw_tris(0,r));
                                    B = nw_pts(nw_tris(1,r));
                                    C = nw_pts(nw_tris(2,r));
                                    
                                    double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                                    double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                                    double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                                    
                                    double normal_x = -SDn_x(i,j,k);
                                    double normal_y = -SDn_y(i,j,k);
                                    double normal_z = -SDn_z(i,j,k);
                                    
                                    if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                        P = A;
                                        A = C;
                                        C = P;
                                    }
                                    if (Q == A) {
                                        if (R == B) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                        if (R == C) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                    }
                                    if (Q == B) {
                                        if (R == A) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                        if (R == C) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                    }
                                    if (Q == C) {
                                        if (R == B) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                        if (R == A) {
                                            A.x += 1.0*Dm->iproc()*(Nx-2);
                                            A.y += 1.0*Dm->jproc()*(Nx-2);
                                            A.z += 1.0*Dm->kproc()*(Nx-2);
                                            B.x += 1.0*Dm->iproc()*(Nx-2);
                                            B.y += 1.0*Dm->jproc()*(Nx-2);
                                            B.z += 1.0*Dm->kproc()*(Nx-2);
                                            C.x += 1.0*Dm->iproc()*(Nx-2);
                                            C.y += 1.0*Dm->jproc()*(Nx-2);
                                            C.z += 1.0*Dm->kproc()*(Nx-2);
                                            nwsMyMesh::FaceIterator fi = vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        }
                                    }
                                }
                            }
                            nwsComputeQuantities(m_nws);
                        }
                    }
                    // wn interface averages
                }
            }
        }
    }
    
    aws_sub.fill(0);
    ans_sub.fill(0);
    wsComputeQuantities(m_ws);
    nsComputeQuantities(m_ns);
    awn_sub.fill(0);
    nwComputeQuantities(m_nw);
    
    
    if (efawns_count == 0) efawns_count = 1;
    if (efawns > 0) { efawns /= double(efawns_count); efawns = 180.0 - efawns; } // if phase affinity < 0

    
    for (int k=kmin; k<kmax; k++) {
        const int k2 = std::min<int>( ( k - kmin) / N2z, subdivide[2] - 1 );
        for (int j=1; j<Ny-1; j++){
            const int j2 = std::min<int>( ( j - 1 ) / N2y, subdivide[1] - 1 );
            for (int i=1; i<Nx-1; i++) {
                const int i2 = std::min<int>( ( i - 1 ) / N2x, subdivide[0] - 1 );
                vnw_tmp.fill( 0 );  vws_tmp.fill( 0 );  vns_tmp.fill( 0 );  vnws_tmp.fill( 0 );
                Gwn_tmp.fill( 0 );  Gws_tmp.fill( 0 );  Gns_tmp.fill( 0 );  Gwns_tmp.fill( 0 );
                wwn_tmp.fill( 0 ); nwOwnO_tmp.fill( 0 ); nwOwsO_tmp.fill( 0 );
                
                // Compute volume averages
                int n = i + j*Nx + k*Nx*Ny;
                double vs = VFmask(n);  
                double vf = 1.0 - vs;
                double phase_val = SDn(n);  
                double delPhi_val = DelPhi(n);
                char id_val = Dm->id[n];

                if (vf == 1.0) {
                    if ( phase_val > 0 ) {
                        nwp_volume_sub(i2,j2,k2) += 1.0;
                        an_sub(i2,j2,k2) += 1.0;
                        porosity_sub(i2,j2,k2) += 1.0;
                        // velocity
                        vn_sub(i2,j2,k2,0) += Vel_x(n);
                        vn_sub(i2,j2,k2,1) += Vel_y(n);
                        vn_sub(i2,j2,k2,2) += Vel_z(n);
                        if (delPhi_val < 1e-4){
                            vol_n_sub(i2,j2,k2) += 1.0;
                            pn_sub(i2,j2,k2) += Press(n);
                        }
                    }
                    else {
                        wp_volume_sub(i2,j2,k2) += 1.0;
                        aw_sub(i2,j2,k2) += 1.0;
                        porosity_sub(i2,j2,k2) += 1.0;
                        // velocity
                        vw_sub(i2,j2,k2,0) += Vel_x(n);
                        vw_sub(i2,j2,k2,1) += Vel_y(n);
                        vw_sub(i2,j2,k2,2) += Vel_z(n);
                        if (delPhi_val < 1e-4){
                            vol_w_sub(i2,j2,k2) += 1.0;
                            pw_sub(i2,j2,k2) += Press(n);
                        }
                    } 
                }
                if (vf < 1.0 && vf > 0.0) {
                    /* Add to the overall volume */
                    if (id_val == 0) { // inactive lattice site
                        if ( phase_val > 0 ) {
                            nwp_volume_sub(i2,j2,k2) += vf;
                            an_sub(i2,j2,k2) += vf;
                            porosity_sub(i2,j2,k2) += vf;
                            if (delPhi_val < 1e-4) {
                                vol_n_sub(i2,j2,k2) += vf;
                            }
                        }
                        else {
                            wp_volume_sub(i2,j2,k2) += vf;
                            aw_sub(i2,j2,k2) += vf;
                            porosity_sub(i2,j2,k2) += vf;
                            if (delPhi_val < 1e-4) {
                                vol_w_sub(i2,j2,k2) += vf;
                            }
                        }
                    }
                    else {
                        if ( phase_val > 0 ) {
                            nwp_volume_sub(i2,j2,k2) += vf;
                            an_sub(i2,j2,k2) += vf;
                            porosity_sub(i2,j2,k2) += vf;
                            // velocity
                            vn_sub(i2,j2,k2,0) += Vel_x(n);
                            vn_sub(i2,j2,k2,1) += Vel_y(n);
                            vn_sub(i2,j2,k2,2) += Vel_z(n);
                            if (delPhi_val < 1e-4){
                                vol_n_sub(i2,j2,k2) += vf;
                                pn_sub(i2,j2,k2) += Press(n);
                            }
                        }
                        else {
                            wp_volume_sub(i2,j2,k2) += vf;
                            aw_sub(i2,j2,k2) += vf;
                            porosity_sub(i2,j2,k2) += vf;
                            // velocity
                            vw_sub(i2,j2,k2,0) += Vel_x(n);
                            vw_sub(i2,j2,k2,1) += Vel_y(n);
                            vw_sub(i2,j2,k2,2) += Vel_z(n);
                            if (delPhi_val < 1e-4){
                                vol_w_sub(i2,j2,k2) += vf;
                                pw_sub(i2,j2,k2) += Press(n);
                            }
                        }
                    }
                } 

                //if (vs > 0.0)  porosity_sub(i2,j2,k2) += vs;
                
                //...........................................................................
                n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
                n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
                //...........................................................................
                pmmc_ConstructLocalCube(SDs, SDn, solid_isovalue, fluid_isovalue,
                                        nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                        local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                        n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                        n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                                        i, j, k, Nx, Ny, Nz);
                //.......................................................................................
                
                if (n_nw_tris > 0) {
                    // Write the triangle lists to text file
                    for (int r=0;r<n_nw_tris;r++) {

                        A = nw_pts(nw_tris(0,r));
                        B = nw_pts(nw_tris(1,r));
                        C = nw_pts(nw_tris(2,r));


                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                        double normal_x = -SDn_x(i,j,k);
                        double normal_y = -SDn_y(i,j,k);
                        double normal_z = -SDn_z(i,j,k);

                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }

                        // Remap the points
                        A.x += 1.0*Dm->iproc()*(Nx-2);
                        A.y += 1.0*Dm->jproc()*(Nx-2);
                        A.z += 1.0*Dm->kproc()*(Nx-2);
                        B.x += 1.0*Dm->iproc()*(Nx-2);
                        B.y += 1.0*Dm->jproc()*(Nx-2);
                        B.z += 1.0*Dm->kproc()*(Nx-2);
                        C.x += 1.0*Dm->iproc()*(Nx-2);
                        C.y += 1.0*Dm->jproc()*(Nx-2);
                        C.z += 1.0*Dm->kproc()*(Nx-2);
                        nwMyMesh::FaceIterator fi = vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                    }
                }
                if ( n_nw_tris >0 ) {

                    double value = pmmc_CubeSurfaceOrientation(Gwn_tmp,nw_pts,nw_tris,n_nw_tris);    // computes normal vector, triangle area, and orientation tensor
                    for (int m=0; m<6; m++)  Gwn_sub(i2,j2,k2,m) += Gwn_tmp(m);

                    pmmc_na(SDn_x, SDn_y, SDn_z,  nwOwnO_tmp, nw_pts, nw_tris,  i,  j,  k, n_nw_pts, n_nw_tris);
                    for (int m=0; m<3; m++)   nwOwnO_sub(i2,j2,k2,m) += nwOwnO_tmp(m);

                    pwn_sub(i2,j2,k2) += pmmc_CubeSurfaceInterpValue(CubeValues,Press,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);

                    // Compute the normal speed of the interface
                    Eta_sub(i2,j2,k2) += pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris, NormalVector, InterfaceSpeed, wwn_tmp, i, j, k, n_nw_pts, n_nw_tris);

                    wwn_sub(i2,j2,k2,0) += wwn_tmp(0);
                    wwn_sub(i2,j2,k2,1) += wwn_tmp(1);
                    wwn_sub(i2,j2,k2,2) += wwn_tmp(2);


                    Eta_sub1(i2,j2,k2) += pmmc_interface_velocity_dot_nn(Vel_x, Vel_y, Vel_z, SDn_x, SDn_y, SDn_z, nw_pts, nw_tris,  i,  j,  k, n_nw_pts, n_nw_tris);
                }

            }
        }
    }
    
    vol_n = vol_n_sub.sum();
    vol_w = vol_w_sub.sum();
    wp_volume = wp_volume_sub.sum();
    nwp_volume = nwp_volume_sub.sum();

    porosity = porosity_sub.sum();
    an = an_sub.sum();
    aw = aw_sub.sum();
    pn = pn_sub.sum();
    pw = pw_sub.sum();
    pwn = pwn_sub.sum();
    lwns = lwns_sub.sum();
    euler = euler_sub.sum();
    Eta = Eta_sub.sum();
    Eta1 = Eta_sub1.sum();

    // Computed in XXComputeMeshQuantities function
    /*aws_sub.fill(0);
    ans_sub.fill(0);
    awn_sub.fill(0);
    
    wsComputeQuantities(m_ws);
    nsComputeQuantities(m_ns);
    nwComputeQuantities(m_nw);
    
    // Perform the local reductions

    Jwnwwndnw = Jwnwwndnw_sub.sum();
    efawns = efawns_sub.sum();
    wwnsdnwn = wwnsdnwn_sub.sum();*/
    
    
    for ( size_t i=0; i<subdivide[0]; i++) {
        for ( size_t j=0; j<subdivide[1]; j++) {
            for ( size_t k=0; k<subdivide[2]; k++) {
                for ( int m=0; m<3; m++) {
                    
                    // Needs to be updated
                    vws(m) += vws_sub(i,j,k,m);
                    vns(m) += vns_sub(i,j,k,m);
                    vnw(m) += vnw_sub(i,j,k,m);
                    wwn(m) += wwn_sub(i,j,k,m);
                    vnws(m) += vnws_sub(i,j,k,m);
                    vw(m) += vw_sub(i,j,k,m);
                    vn(m) += vn_sub(i,j,k,m);
                    
                }
                for ( int m=0; m<6; m++) {
                    Gwn(m) += Gwn_sub(i,j,k,m);
                    Gws(m) += Gws_sub(i,j,k,m);
                    Gns(m) += Gns_sub(i,j,k,m);
                    Gwns(m) += Gwns_sub(i,j,k,m);
                }
            }
        }
    }
}


void TwoPhase::MeshComputations(){}


void TwoPhase::AssignComponentLabels()
{
    const int LabelNWP=1;
    const int LabelWP=2;
    NumberComponents_NWP = ComputeGlobalBlobIDs(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->rank_info,SDs,SDn,solid_isovalue,fluid_isovalue,Label_NWP,Dm->Comm);
}

void TwoPhase::ComponentAverages()
{

    int i,j,k,n;
    int kmin,kmax;
	int LabelWP,LabelNWP;
	double TempLocal;

    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    
    ComponentAverages_WP.resize(BLOB_AVG_COUNT,NumberComponents_WP);
    ComponentAverages_NWP.resize(BLOB_AVG_COUNT,NumberComponents_NWP);
    
    ComponentAverages_WP.fill(0.0);
    ComponentAverages_NWP.fill(0.0);
    
    //if (Dm->rank()==0){
    //    printf("Number of wetting phase components is %i \n",NumberComponents_WP);
    //    printf("Number of non-wetting phase components is %i \n",NumberComponents_NWP);
    //}
    
    /* Component analysis code */
    kmin=1; kmax=Nz-1;
	if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=4;
	if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-4;
	for (k=kmin; k<kmax; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				LabelWP=GetCubeLabel(i,j,k,Label_WP);
				LabelNWP=GetCubeLabel(i,j,k,Label_NWP);
				n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
				n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
				// Initialize the averaged quantities
				awn = aws = ans = lwns = 0.0;
				vnw(0) = vnw(1) = vnw(2) = 0.0;
				vnws(0) = vnws(1) = vnws(2) = 0.0;
				Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
				Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
				Gws(0) = Gws(1) = Gws(2) = 0.0;
				Gws(3) = Gws(4) = Gws(5) = 0.0;
				Gns(0) = Gns(1) = Gns(2) = 0.0;
				Gns(3) = Gns(4) = Gns(5) = 0.0;
				Jwn = Kwn = efawns = 0.0;
				//trawn=trJwn=0.0;
				euler=0.0;
				//...........................................................................
				//...........................................................................
				// Compute volume averages
				for (int p=0;p<8;p++){
					n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
					if ( Dm->id[n] > 0 ){
						// 1-D index for this cube corner
						// compute the norm of the gradient of the phase indicator field
						// Compute the non-wetting phase volume contribution

						if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0.0 && !(LabelNWP < 0) ){
							// volume
							ComponentAverages_NWP(VOL,LabelNWP) += 0.125;
							// velocity
							ComponentAverages_NWP(VX,LabelNWP) += 0.125*Vel_x(n);
							ComponentAverages_NWP(VY,LabelNWP) += 0.125*Vel_y(n);
							ComponentAverages_NWP(VZ,LabelNWP) += 0.125*Vel_z(n);
							// center of mass
							ComponentAverages_NWP(CMX,LabelNWP) += 0.125*(i+cube[p][0]+Dm->iproc()*Nx);
							ComponentAverages_NWP(CMY,LabelNWP) += 0.125*(j+cube[p][1]+Dm->jproc()*Ny);
							ComponentAverages_NWP(CMZ,LabelNWP) += 0.125*(k+cube[p][2]+Dm->kproc()*Nz);

							// twice the kinetic energy
							ComponentAverages_NWP(VSQ,LabelNWP) += 0.125*(Vel_x(n)*Vel_x(n)+Vel_y(n)*Vel_y(n)+Vel_z(n)*Vel_z(n));

							// volume the for pressure averaging excludes the interfacial region
							if (DelPhi(n) < 1e-4 ){
								ComponentAverages_NWP(TRIMVOL,LabelNWP) += 0.125;
								ComponentAverages_NWP(PRS,LabelNWP) += 0.125*Press(n);
							}
						}
						else if (!(LabelWP < 0)){
							ComponentAverages_WP(VOL,LabelWP) += 0.125;
							// velocity
							ComponentAverages_WP(VX,LabelWP) += 0.125*Vel_x(n);
							ComponentAverages_WP(VY,LabelWP)+= 0.125*Vel_y(n);
							ComponentAverages_WP(VZ,LabelWP) += 0.125*Vel_z(n);
							// Center of mass
							ComponentAverages_WP(CMX,LabelWP) += 0.125*(i+cube[p][0]+Dm->iproc()*Nx);
							ComponentAverages_WP(CMY,LabelWP) += 0.125*(j+cube[p][1]+Dm->jproc()*Ny);
							ComponentAverages_WP(CMZ,LabelWP) += 0.125*(k+cube[p][2]+Dm->kproc()*Nz);
							// twice the kinetic energy
							ComponentAverages_WP(VSQ,LabelWP) += 0.125*(Vel_x(n)*Vel_x(n)+Vel_y(n)*Vel_y(n)+Vel_z(n)*Vel_z(n));

							// volume the for pressure averaging excludes the interfacial region
							if (DelPhi(n) < 1e-4){
								ComponentAverages_WP(TRIMVOL,LabelWP) += 0.125;
								ComponentAverages_WP(PRS,LabelWP) += 0.125*Press(n);
							}
						}
					}
				}
				//...........................................................................
				// Construct the interfaces and common curve
				pmmc_ConstructLocalCube(SDs, SDn, solid_isovalue, fluid_isovalue,
						nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
						local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
						n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
						n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
						i, j, k, Nx, Ny, Nz);

				//...........................................................................
				// wn interface averages
				if (n_nw_pts>0  && LabelNWP >=0 && LabelWP >=0 ){
					// Mean curvature
                    //printf("Pre GWN size: %i \n",Gwn.size());
                    //TempLocal = pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
                    //printf("Post GWN size: %i \n",Gwn.size());
					//ComponentAverages_WP(AWN,LabelWP) += TempLocal;
					//ComponentAverages_NWP(AWN,LabelNWP) += TempLocal;

					TempLocal = pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
					ComponentAverages_WP(JWN,LabelWP) += TempLocal;
					ComponentAverages_NWP(JWN,LabelNWP) += TempLocal;

					// Trimmed Mean curvature
					//pmmc_CubeTrimSurfaceInterpValues(CubeValues,MeanCurvature,SDs,nw_pts,nw_tris,Values,DistanceValues,
					//		i,j,k,n_nw_pts,n_nw_tris,trimdist,trawn,trJwn);
					//ComponentAverages_WP(TRAWN,LabelWP) += trawn;
					//ComponentAverages_WP(TRJWN,LabelWP) += trJwn;
					//ComponentAverages_NWP(TRAWN,LabelNWP) += trawn;
					//ComponentAverages_NWP(TRJWN,LabelNWP) += trJwn;

					// Gaussian curvature
					TempLocal = pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
					ComponentAverages_WP(KWN,LabelWP) += TempLocal;
					ComponentAverages_NWP(KWN,LabelNWP) += TempLocal;

					// Compute the normal speed of the interface
					pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris,
							NormalVector, InterfaceSpeed, vnw, i, j, k, n_nw_pts, n_nw_tris);
					ComponentAverages_WP(VNWX,LabelWP) += vnw(0);
					ComponentAverages_WP(VNWY,LabelWP) += vnw(1);
					ComponentAverages_WP(VNWZ,LabelWP) += vnw(2);
					ComponentAverages_NWP(VNWX,LabelNWP) += vnw(0);
					ComponentAverages_NWP(VNWY,LabelNWP) += vnw(1);
					ComponentAverages_NWP(VNWZ,LabelNWP) += vnw(2);

					// Interfacial Area
					TempLocal = pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
					ComponentAverages_WP(AWN,LabelWP) += TempLocal;
					ComponentAverages_NWP(AWN,LabelNWP) += TempLocal;

					ComponentAverages_WP(GWNXX,LabelWP) += Gwn(0);
					ComponentAverages_WP(GWNYY,LabelWP) += Gwn(1);
					ComponentAverages_WP(GWNZZ,LabelWP) += Gwn(2);
					ComponentAverages_WP(GWNXY,LabelWP) += Gwn(3);
					ComponentAverages_WP(GWNXZ,LabelWP) += Gwn(4);
					ComponentAverages_WP(GWNYZ,LabelWP) += Gwn(5);

					ComponentAverages_NWP(GWNXX,LabelNWP) += Gwn(0);
					ComponentAverages_NWP(GWNYY,LabelNWP) += Gwn(1);
					ComponentAverages_NWP(GWNZZ,LabelNWP) += Gwn(2);
					ComponentAverages_NWP(GWNXY,LabelNWP) += Gwn(3);
					ComponentAverages_NWP(GWNXZ,LabelNWP) += Gwn(4);
					ComponentAverages_NWP(GWNYZ,LabelNWP) += Gwn(5);

				}
				//...........................................................................
				// Common curve averages
				if (n_local_nws_pts > 0  && LabelNWP >=0 && LabelWP >=0){
					// Contact angle
					TempLocal = pmmc_CubeContactAngle(CubeValues,Values,SDn_x,SDn_y,SDn_z,SDs_x,SDs_y,SDs_z,
							local_nws_pts,i,j,k,n_local_nws_pts);
					ComponentAverages_WP(CWNS,LabelWP) += TempLocal;
					ComponentAverages_NWP(CWNS,LabelNWP) += TempLocal;

					// Kinematic velocity of the common curve
					pmmc_CommonCurveSpeed(CubeValues, dPdt, vnws, SDn_x, SDn_y, SDn_z,SDs_x,SDs_y,SDs_z,
							local_nws_pts,i,j,k,n_local_nws_pts);
					ComponentAverages_WP(VNWSX,LabelWP) += vnws(0);
					ComponentAverages_WP(VNWSY,LabelWP) += vnws(1);
					ComponentAverages_WP(VNWSZ,LabelWP) += vnws(2);
					ComponentAverages_NWP(VNWSX,LabelNWP) += vnws(0);
					ComponentAverages_NWP(VNWSY,LabelNWP) += vnws(1);
					ComponentAverages_NWP(VNWSZ,LabelNWP) += vnws(2);

					// Curvature of the common curve
                    pmmc_CurveCurvature(SDn, SDs, SDn_x, SDn_y, SDn_z, SDs_x, SDs_y,
                                        SDs_z, KNwns_values, KGwns_values, KNwns, KGwns,
                                        nws_pts, n_nws_pts, i, j, k);
					ComponentAverages_WP(KNWNS,LabelWP) += KNwns;
					ComponentAverages_WP(KGWNS,LabelWP) += KGwns;
					ComponentAverages_NWP(KNWNS,LabelNWP) += KNwns;
					ComponentAverages_NWP(KGWNS,LabelNWP) += KGwns;

					// Length of the common curve
					TempLocal = pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
					ComponentAverages_NWP(LWNS,LabelNWP) += TempLocal;
                    ComponentAverages_WP(LWNS,LabelWP) += TempLocal;
				}
				//...........................................................................
				// Solid interface averages
				//if (n_local_sol_pts > 0  && LabelWP >=0){
					//As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);
					// Compute the surface orientation and the interfacial area

				if (n_ns_pts > 0  && LabelNWP >=0 ){

                    //As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);

					TempLocal = pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
					ComponentAverages_NWP(AWNS,LabelNWP) += TempLocal;
                    ComponentAverages_NWP(GXX,LabelNWP) += Gns(0);
					ComponentAverages_NWP(GYY,LabelNWP) += Gns(1);
					ComponentAverages_NWP(GZZ,LabelNWP) += Gns(2);
					ComponentAverages_NWP(GXY,LabelNWP) += Gns(3);
					ComponentAverages_NWP(GXZ,LabelNWP) += Gns(4);
					ComponentAverages_NWP(GYZ,LabelNWP) += Gns(5);
				}

                if (n_ws_pts > 0 && LabelWP >= 0){

				    TempLocal = pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
					ComponentAverages_WP(AWNS,LabelWP) += TempLocal;

                    ComponentAverages_WP(GXX,LabelWP) += Gws(0);
					ComponentAverages_WP(GYY,LabelWP) += Gws(1);
				    ComponentAverages_WP(GZZ,LabelWP) += Gws(2);
					ComponentAverages_WP(GXY,LabelWP) += Gws(3);
					ComponentAverages_WP(GXZ,LabelWP) += Gws(4);
					ComponentAverages_WP(GYZ,LabelWP) += Gws(5);
				}

				//...........................................................................


				/* Compute the Euler characteristic
				 * count all vertices, edges and faces (triangles)
				 * Euler Number = vertices - edges + faces
				 * double geomavg_EulerCharacteristic(PointList, PointCount, TriList, TriCount);
				 */
				if (LabelNWP >= 0 ){
				  // find non-wetting phase interface
//					double euler;
					n_nw_pts=n_nw_tris=0;
					geomavg_MarchingCubes(SDn,fluid_isovalue,i,j,k,nw_pts,n_nw_pts,nw_tris,n_nw_tris);
                    euler = geomavg_EulerCharacteristic(nw_pts,nw_tris,n_nw_pts,n_nw_tris,i,j,k);
					ComponentAverages_NWP(EULER,LabelNWP) += euler;
					// Compute Euler characteristic from integral of gaussian curvature
					//euler = pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,
					//		i,j,k,n_nw_pts,n_nw_tris);
					//ComponentAverages_NWP(INTCURV,LabelNWP) += euler;

					// Compute the Euler characteristic from vertices - faces + edges
                    n_ns_pts=n_ns_tris=0;
					geomavg_MarchingCubes(SDn,fluid_isovalue,i,j,k,ns_pts,n_ns_pts,ns_tris,n_ns_tris);
					euler = geomavg_EulerCharacteristic(ns_pts,ns_tris,n_ns_pts,n_ns_tris,i,j,k);
					ComponentAverages_NWP(EULER,LabelNWP) += euler;                    

				}
			}
		}
	}

	MPI_Barrier(Dm->Comm);
	/*if (Dm->rank()==0){
		printf("Component averages computed locally -- reducing result... \n");
	}*/
	// Globally reduce the non-wetting phase averages
	RecvBuffer.resize(BLOB_AVG_COUNT,NumberComponents_NWP);

/*	for (int b=0; b<NumberComponents_NWP; b++){
		MPI_Barrier(Dm->Comm);
		MPI_Allreduce(&ComponentAverages_NWP(0,b),&RecvBuffer(0),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,Dm->Comm);
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(idx);
	}
	*/
	MPI_Barrier(Dm->Comm);
	MPI_Allreduce(ComponentAverages_NWP.data(),RecvBuffer.data(),BLOB_AVG_COUNT*NumberComponents_NWP,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	//	MPI_Reduce(ComponentAverages_NWP.data(),RecvBuffer.data(),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,0,Dm->Comm);
  
	/*if (Dm->rank()==0){
		printf("rescaling... \n");
	}*/

	for (int b=0; b<NumberComponents_NWP; b++){
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(idx,b);
	}

	for (int b=0; b<NumberComponents_NWP; b++){
        
		if (ComponentAverages_NWP(VOL,b) > 0.0){
			double Vn,pn,awn,ans,Jwn,Kwn,lwns,cwns,vsq;

			Vn = ComponentAverages_NWP(VOL,b);
			awn = ComponentAverages_NWP(AWN,b);
			ans = ComponentAverages_NWP(AWNS,b);
			vn(0) = ComponentAverages_NWP(VX,b)/Vn;
			vn(1) = ComponentAverages_NWP(VY,b)/Vn;
			vn(2) = ComponentAverages_NWP(VZ,b)/Vn;
			vsq = ComponentAverages_NWP(VSQ,b)/Vn;
            NULL_USE(ans);

			if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
				pn = ComponentAverages_NWP(PRS,b)/ComponentAverages_NWP(TRIMVOL,b);
			}
			else pn = 0.0;

			if (awn != 0.0){
				Jwn = ComponentAverages_NWP(JWN,b)/awn;
				Kwn = ComponentAverages_NWP(KWN,b)/awn;
				vnw(0) = ComponentAverages_NWP(VNWX,b)/awn;
				vnw(1) = ComponentAverages_NWP(VNWY,b)/awn;
				vnw(2) = ComponentAverages_NWP(VNWZ,b)/awn;
				Gwn(0) = ComponentAverages_NWP(GWNXX,b)/awn;
				Gwn(1) = ComponentAverages_NWP(GWNYY,b)/awn;
				Gwn(2) = ComponentAverages_NWP(GWNZZ,b)/awn;
				Gwn(3) = ComponentAverages_NWP(GWNXY,b)/awn;
				Gwn(4) = ComponentAverages_NWP(GWNXZ,b)/awn;
				Gwn(5) = ComponentAverages_NWP(GWNYZ,b)/awn;
			}
			else Jwn=Kwn=0.0;

            if (ans != 0.0){
                Gns(0) = ComponentAverages_NWP(GXX,b)/ans;
				Gns(1) = ComponentAverages_NWP(GYY,b)/ans;
				Gns(2) = ComponentAverages_NWP(GZZ,b)/ans;
				Gns(3) = ComponentAverages_NWP(GXY,b)/ans;
				Gns(4) = ComponentAverages_NWP(GXZ,b)/ans;
				Gns(5) = ComponentAverages_NWP(GYZ,b)/ans;
            }
            else Gns(0) = Gns(1) = Gns(2) = Gns(3) = Gns(4) = Gns(5) = 0.0;

			lwns = ComponentAverages_NWP(LWNS,b);
			if (lwns != 0.0){
				cwns = ComponentAverages_NWP(CWNS,b)/lwns;
				vnws(0) = ComponentAverages_NWP(VNWSX,b)/lwns;
				vnws(1) = ComponentAverages_NWP(VNWSY,b)/lwns;
				vnws(2) = ComponentAverages_NWP(VNWSZ,b)/lwns;
                KNwns = ComponentAverages_NWP(KNWNS,b)/lwns;
                KGwns = ComponentAverages_NWP(KGWNS,b)/lwns;
			}
			else  {
                cwns=0.0; 
                KNwns = 0.0; 
                KGwns = 0.0;
            }

			ComponentAverages_NWP(PRS,b) = pn;
			ComponentAverages_NWP(VX,b) = vn(0);
			ComponentAverages_NWP(VY,b) = vn(1);
			ComponentAverages_NWP(VZ,b) = vn(2);
			ComponentAverages_NWP(VSQ,b) = vsq;

			ComponentAverages_NWP(JWN,b) = Jwn;
			ComponentAverages_NWP(KWN,b) = Kwn;
			ComponentAverages_NWP(VNWX,b) = vnw(0);
			ComponentAverages_NWP(VNWY,b) = vnw(1);
			ComponentAverages_NWP(VNWZ,b) = vnw(2);
			
			ComponentAverages_NWP(GWNXX,b) = Gwn(0);
			ComponentAverages_NWP(GWNYY,b) = Gwn(1);
			ComponentAverages_NWP(GWNZZ,b) = Gwn(2);
			ComponentAverages_NWP(GWNXY,b) = Gwn(3);
			ComponentAverages_NWP(GWNXZ,b) = Gwn(4);
			ComponentAverages_NWP(GWNYZ,b) = Gwn(5);

            ComponentAverages_NWP(GXX,b) = Gns(0);
			ComponentAverages_NWP(GYY,b) = Gns(1);
			ComponentAverages_NWP(GZZ,b) = Gns(2);
			ComponentAverages_NWP(GXY,b) = Gns(3);
			ComponentAverages_NWP(GXZ,b) = Gns(4);
			ComponentAverages_NWP(GYZ,b) = Gns(5);

			ComponentAverages_NWP(CWNS,b) = cwns;
			ComponentAverages_NWP(VNWSX,b) = vnws(0);
			ComponentAverages_NWP(VNWSY,b) = vnws(1);
			ComponentAverages_NWP(VNWSZ,b) = vnws(2);
            ComponentAverages_NWP(KNWNS,b) = KNwns;
            ComponentAverages_NWP(KGWNS,b) = KGwns;

			ComponentAverages_NWP(CMX,b) /= Vn;
			ComponentAverages_NWP(CMY,b) /= Vn;
			ComponentAverages_NWP(CMZ,b) /= Vn;

			ComponentAverages_NWP(EULER,b) /= (2*PI);

		}

	}

	/*if (Dm->rank()==0){
		printf("reduce WP averages... \n");
	}*/

	// reduce the wetting phase averages
	for (int b=0; b<NumberComponents_WP; b++){
		MPI_Barrier(Dm->Comm);
//		MPI_Allreduce(&ComponentAverages_WP(0,b),RecvBuffer.data(),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,Dm->Comm);
		MPI_Reduce(&ComponentAverages_WP(0,b),RecvBuffer.data(),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,0,Dm->Comm);
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_WP(idx,b)=RecvBuffer(idx);
	}
    
	for (int b=0; b<NumberComponents_WP; b++){
		if (ComponentAverages_WP(VOL,b) > 0.0){
			double Vw,pw,awn,aws,Jwn,Kwn,lwns,cwns,vsq;
			Vw = ComponentAverages_WP(VOL,b);
			awn = ComponentAverages_WP(AWN,b);
			aws = ComponentAverages_WP(AWNS,b);
			vw(0) = ComponentAverages_WP(VX,b)/Vw;
			vw(1) = ComponentAverages_WP(VY,b)/Vw;
			vw(2) = ComponentAverages_WP(VZ,b)/Vw;
			vsq = ComponentAverages_WP(VSQ,b)/Vw;
            //NULL_USE(ans);

			if (ComponentAverages_WP(TRIMVOL,b) > 0.0){
				pw = ComponentAverages_WP(PRS,b)/ComponentAverages_WP(TRIMVOL,b);
			}
			else pw = 0.0;

			if (awn != 0.0){
				Jwn = ComponentAverages_WP(JWN,b)/awn;
				Kwn = ComponentAverages_WP(KWN,b)/awn;
				vnw(0) = ComponentAverages_WP(VNWX,b)/awn;
				vnw(1) = ComponentAverages_WP(VNWY,b)/awn;
				vnw(2) = ComponentAverages_WP(VNWZ,b)/awn;
				Gwn(0) = ComponentAverages_WP(GWNXX,b)/awn;
				Gwn(1) = ComponentAverages_WP(GWNYY,b)/awn;
				Gwn(2) = ComponentAverages_WP(GWNZZ,b)/awn;
				Gwn(3) = ComponentAverages_WP(GWNXY,b)/awn;
				Gwn(4) = ComponentAverages_WP(GWNXZ,b)/awn;
				Gwn(5) = ComponentAverages_WP(GWNYZ,b)/awn;
			}
			else Jwn=Kwn=0.0;

            if (aws != 0.0){
                Gws(0) = ComponentAverages_WP(GXX,b)/aws;
				Gws(1) = ComponentAverages_WP(GYY,b)/aws;
				Gws(2) = ComponentAverages_WP(GZZ,b)/aws;
				Gws(3) = ComponentAverages_WP(GXY,b)/aws;
				Gws(4) = ComponentAverages_WP(GXZ,b)/aws;
				Gws(5) = ComponentAverages_WP(GYZ,b)/aws;
            }
            else Gws(0) = Gws(1) = Gws(2) = Gws(3) = Gws(4) = Gws(5) = 0.0;

			lwns = ComponentAverages_WP(LWNS,b);
			if (lwns != 0.0){
				cwns = ComponentAverages_WP(CWNS,b)/lwns;
				vnws(0) = ComponentAverages_WP(VNWSX,b)/lwns;
				vnws(1) = ComponentAverages_WP(VNWSY,b)/lwns;
				vnws(2) = ComponentAverages_WP(VNWSZ,b)/lwns;
                KNwns = ComponentAverages_WP(KNWNS,b)/lwns;
                KGwns = ComponentAverages_WP(KGWNS,b)/lwns;
			}
			else  {
                cwns=0.0; 
                KNwns = 0.0; 
                KGwns = 0.0;
            }

			ComponentAverages_WP(PRS,b) = pw;
			ComponentAverages_WP(VX,b) = vw(0);
			ComponentAverages_WP(VY,b) = vw(1);
			ComponentAverages_WP(VZ,b) = vw(2);
			ComponentAverages_WP(VSQ,b) = vsq;

			ComponentAverages_WP(JWN,b) = Jwn;
			ComponentAverages_WP(KWN,b) = Kwn;
			ComponentAverages_WP(VNWX,b) = vnw(0);
			ComponentAverages_WP(VNWY,b) = vnw(1);
			ComponentAverages_WP(VNWZ,b) = vnw(2);
			
			ComponentAverages_WP(GWNXX,b) = Gwn(0);
			ComponentAverages_WP(GWNYY,b) = Gwn(1);
			ComponentAverages_WP(GWNZZ,b) = Gwn(2);
			ComponentAverages_WP(GWNXY,b) = Gwn(3);
			ComponentAverages_WP(GWNXZ,b) = Gwn(4);
			ComponentAverages_WP(GWNYZ,b) = Gwn(5);

            ComponentAverages_WP(GXX,b) = Gws(0);
			ComponentAverages_WP(GYY,b) = Gws(1);
			ComponentAverages_WP(GZZ,b) = Gws(2);
			ComponentAverages_WP(GXY,b) = Gws(3);
			ComponentAverages_WP(GXZ,b) = Gws(4);
			ComponentAverages_WP(GYZ,b) = Gws(5);

			ComponentAverages_WP(CWNS,b) = cwns;
			ComponentAverages_WP(VNWSX,b) = vnws(0);
			ComponentAverages_WP(VNWSY,b) = vnws(1);
			ComponentAverages_WP(VNWSZ,b) = vnws(2);
            ComponentAverages_WP(KNWNS,b) = KNwns;
            ComponentAverages_WP(KGWNS,b) = KGwns;


			//ComponentAverages_WP(TRJWN,b) = trJwn;
		}
	}        
}
void TwoPhase::WriteSurfaces(int logcount){}

void TwoPhase::ReduceGlobal()
{
    
    double iVol_global = 1.0/double(Volume);
    
    // Sum local averages to global
    MPI_Barrier(Dm->Comm);
    MPI_Allreduce(&porosity,&porosity_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    // MPI_Allreduce(&dPdt,&dPdt_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    //  dPdt_global /= double(npx * npy * npz);
    
    
    /* Volume Reductions */
    MPI_Allreduce(&nwp_volume,&nwp_volume_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wp_volume,&wp_volume_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&an,&an_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&aw,&aw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&Eta,&Eta_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Eta1,&Eta_global1,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    //MPI_Allreduce(&awn,&awn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&euler,&euler_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    //   std::cout << "euler_global=" << euler_global << std::endl;
    
    MPI_Allreduce(&vol_w,&vol_w_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vol_n,&vol_n_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&local_nw_phase_volume,&global_nw_phase_volume,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&mass_local,&mass_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    if (awn_global == 0.0) {
        Jwn_global = 0.0;
        Kwn_global = 0.0;
    }
    if (aws_global == 0.0){
        Jws_global = 0.0;
        Kws_global = 0.0;
    }
    if (ans_global == 0.0) {
        Jns_global = 0.0;
        Kns_global = 0.0;
    }
    
    /* Common curve curvatures  */
    MPI_Allreduce(&KGwns,&KGwns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&KNwns,&KNwns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    /* Common Curve Reduction */
    MPI_Allreduce(&efawns,&efawns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    efawns_global /= double(npx * npy * npz);

    /* Pressure */
    MPI_Allreduce(&pw,&pw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&pwn,&pwn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&pn,&pn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);


    if (vol_n_global > 0) {
        pn_global /= vol_n_global;
    }
    else {
        pn_global = 0.0;
    }

    if (vol_w_global > 0) {
        pw_global /= vol_w_global;
    }
    else {
        pw_global = 0.0;
    }
    
    if (awn_global > 0 ) {
        pwn_global /= awn_global;
    }
    else {
        pwn_global= 0.0;
    }
    
    /* Velocity-Related Reductions */
    MPI_Allreduce(&vn(0),&vn_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vw(0),&vw_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwn(0),&wwn_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&vws(0),&vws_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vns(0),&vns_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vnw(0),&vnw_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vnws(0),&vnws_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vnwdnw,&vnwdnw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwndnw,&wwndnw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwnsdnwn,&wwnsdnwn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&nwOwnO(0),&nwOwnO_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    if (wp_volume_global > 0.0){
        vw_global(0) /= wp_volume_global;
        vw_global(1) /= wp_volume_global;
        vw_global(2) /= wp_volume_global;
    }
    else {
        vw_global(0) = 0.0;
        vw_global(1) = 0.0;
        vw_global(2) = 0.0;
    }
    if (nwp_volume_global > 0.0){
        vn_global(0) /= nwp_volume_global;
        vn_global(1) /= nwp_volume_global;
        vn_global(2) /= nwp_volume_global;
    }
    else {
        vn_global(0) = 0.0;
        vn_global(1) = 0.0;
        vn_global(2) = 0.0;
    }
    if (awn_global > 0.0) {
        for (int i=0; i<3; i++) vnw_global(i)/= awn_global;
    }
    else {
        for (int i=0; i<3; i++) vnw_global(i) = 0;
    }
    if (awn_global > 0) {
        wwndnw_global /= awn_global;
    }
    else {
        wwndnw_global = 0;
    }
    
    /* Orientation Reductions */
    MPI_Allreduce(&Gwn(0),&Gwn_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gns(0),&Gns_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gws(0),&Gws_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gwns(0),&Gwns_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    if (awn_global > 0.0) {
        for (int i=0; i<6; i++) Gwn_global(i)/= awn_global;
    }
    else {
        for (int i=0; i<6; i++) Gwn_global(i) = 0.0;
    }
    
    if (ans_global > 0.0) {
        for (int i=0; i<3; i++) vns_global(i) /= ans_global;
        for (int i=0; i<6; i++) Gns_global(i) /= ans_global;
    }
    else {
        for (int i=0; i<3; i++) vns_global(i) = 0.0;
        for (int i=0; i<6; i++) Gns_global(i) = 0.0;
    }
    if (aws_global > 0.0) {
        for (int i=0; i<3; i++) vws_global(i) /= aws_global;
        for (int i=0; i<6; i++) Gws_global(i) /= aws_global;
    }
    else {
        for (int i=0; i<3; i++) vws_global(i) = 0.0;
        for (int i=0; i<6; i++) Gws_global(i) = 0.0;
    }
    
    /* Curvature Reductions */
    MPI_Allreduce(&lwns,&lwns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    if (lwns_global > 0.0){
        //efawns_global /= lwns_global;
        KNwns_global /= lwns_global;
        KGwns_global /= lwns_global;
        for (int i=0; i<3; i++) vnws_global(i) /= lwns_global;
        for (int i=0; i<6; i++) Gwns_global(i) /= lwns_global;
    }
    else {
        efawns_global = 1.0;
        KNwns_global = 0.0;
        KGwns_global = 0.0;
        for (int i=0; i<3; i++) vnws_global(i) = 0.0;
        for (int i=0; i<6; i++) Gwns_global(i) = 0.0;
    }
    
    MPI_Allreduce(&nwOwsO(0),&nwOwsO_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Barrier(Dm->Comm);
    // <v_wn . n_w >_OwnO
    vnwdnw_global = vnwdnw_global*iVol_global;
    Eta_global = Eta_global*iVol_global;
    Eta_global1 = Eta_global1*iVol_global;

    // w^wn
    wwn_global(0) = wwn_global(0)*iVol_global;
    wwn_global(1) = wwn_global(1)*iVol_global;
    wwn_global(2) = wwn_global(2)*iVol_global;
    porosity_global *= iVol_global;
    ans_global = ans_global*iVol_global;
    aws_global = aws_global*iVol_global;
    lwns_global = lwns_global*iVol_global;
    an_global = an_global*iVol_global;
    aw_global = aw_global*iVol_global;

}

void TwoPhase::NonDimensionalize(double D, double viscosity, double IFT)
{
    NULL_USE( viscosity );
    NULL_USE( IFT );
    awn_global *= D;
    ans_global *= D;
    ans_global *= D;
    lwns_global *= D*D;
}

void TwoPhase::PrintAll(int timestep)
{

    if ( TIMELOG_global ){
        
        // double alpha,sphere_diameter,average_rho_input, average_mu;
        //Ca_global = abs(vn_global(2) *average_mu)/(5.796*alpha*1.0);

        Ca_global = abs(vn_global(2)*muA + vw_global(2)*muB)/(5.796*alpha*1.0);
        Re_global = abs(vn_global(2) *1.0 * sphere_diameter)/(average_mu);
        
        if (Ca_global < 1E-15) Ca_global = 0.0;
        if (Re_global < 1E-15) Re_global = 0.0;
        
        double sat_w = 1.0 - nwp_volume_global/(nwp_volume_global+wp_volume_global);
        // double poros = porosity_global;
        // efawns_global = std::acos(efawns_global)*180/PI;

        nw_global(0) = nw_global(1) = nw_global(2) = 0;
        // "time porosity s^w p^w* p^n* eps^{wn}* eps^{ns}* eps^{ws}* J_n^{wn}* K_n^{wn}* l^{wns}* c^{wns} kappa_N^{wns}* kappa_G^{wns}* "
        fprintf(TIMELOG_global,"%i %-15.15E ",timestep,porosity_global);  // porosity and change in surface energy
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",sat_w,pw_global,pn_global);                  // saturation and pressure
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",awn_global,ans_global,aws_global);             // interfacial areas
        fprintf(TIMELOG_global,"%-15.15E %-15.15E ",Jwn_global, Kwn_global);                            // curvature of wn interface
        fprintf(TIMELOG_global,"%-15.15E ",lwns_global);                                            // common curve length
        fprintf(TIMELOG_global,"%.10f ",efawns_global);                                          // average contact angle
        fprintf(TIMELOG_global,"%-15.15E %-15.15E ",KNwns_global, KGwns_global);                        // curvature of wn interface
        
        // ""v_x^w* v_y^w* v_z^w* v_x^n* v_y^n* v_z^n* v_x^{wn}* v_y^{wn}* v_z^{wn}* v_x^{wns}* v_y^{wns}* v_z^{wns}* "
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vw_global(0),vw_global(1),vw_global(2));    // average velocity of w phase
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vn_global(0),vn_global(1),vn_global(2));    // average velocity of n phase
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vnw_global(0),vnw_global(1),vnw_global(2)); // velocity of wn interface
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vnws_global(0),vnws_global(1),vnws_global(2)); // velocity of wn interface
        
        // "G_{xx}^{wn}* G_{yy}^{wn}* G_{zz}^{wn}* G_{xy}^{wn}* G_{xz}^{wn}* G_{yz}^{wn}* "
        // "G_{xx}^{ws}* G_{yy}^{ws}* G_{zz}^{ws}* G_{xy}^{ws}* G_{xz}^{ws}* G_{yz}^{ws}* "
        // "G_{xx}^{ns}* G_{yy}^{ns}* G_{zz}^{ns}* G_{xy}^{ns}* G_{xz}^{ns}* G_{yz}^{ns}* "
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwn_global(0),Gwn_global(1),Gwn_global(2),Gwn_global(3),Gwn_global(4),Gwn_global(5)); // orientation of wn interface
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gws_global(0),Gws_global(1),Gws_global(2),Gws_global(3),Gws_global(4),Gws_global(5)); // orientation of ns interface
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gns_global(0),Gns_global(1),Gns_global(2),Gns_global(3),Gns_global(4),Gns_global(5)); // orientation of ws interface
        
        // "p_n^{wn}* X^n* eps^w eps^n "
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E ", pwn_global, euler_global, aw_global, an_global);
        
        // "G_{xx}^{wns}* G_{yy}^{wns}* G_{zz}^{wns}* G_{xy}^{wns}* G_{xz}^{wns}* G_{yz}^{wns}* "
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwns_global(0),Gwns_global(1),Gwns_global(2),Gwns_global(3),Gwns_global(4),Gwns_global(5)); // orientation of wns interface
        
        // "v_x^{ws}* v_y^{ws}* v_z^{ws}* v_x^{ns}* v_y^{ns}* v_z^{ns}* "
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                vws_global(0),vws_global(1),vws_global(2),vns_global(0),vns_global(1),vns_global(2)); // velocity ws and ns
        
        // "nw_x- nw_y- nw_z- vnwdnw- wwnx- wwny- wwnz- nwOwnOx- nwOwnOy- nwOwnOz- nwOwsOx- nwOwsOy- nwOwsOz- "
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",
                nwOwnO_global2(0),nwOwnO_global2(1),nwOwnO_global2(2));
        fprintf(TIMELOG_global,"%-15.15E ",
                Eta_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",
                wwn_global(0),wwn_global(1),wwn_global(2));
        // "Jwn2 Kwn2"
//        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E\n", Jwn2_global, Kwn2_global, Ca_global, Re_global, an_global1, an_global2, an_global3, Eta_global1, dendt_global,dendt_global3);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E\n", Ca_global, Re_global);
        
        
        
        
        
        
        
        
        fflush(TIMELOG_global);
        
        
    }
    if ( TIMELOG_local ){
        double iVol_nodal = 1.0/nodeVolume;
        
        
        
        porosity *= iVol_nodal;
        euler *= iVol_nodal;
        eulerW *= iVol_nodal;
        eulerS *= iVol_nodal;
        
        
        if (vol_n > 0) {
            pn /= vol_n;
        }
        else {
            pn = 0.0;
        }
        
        if (vol_w > 0) {
            pw /= vol_w;
        }
        else {
            pw = 0.0;
        }
        
        
        if (awn > 0.0) {
            //            Jwn /= awn;
            //            Kwn /= awn;
            //            pwn /= awn;
        }
        else {
            JwnW = 0;
            KwnW = 0;
            Jwn = 0;
            Kwn = 0;
            pwn = 0;
        }
        
        if (aws > 0.0) {
            //            Jws /= aws;
            //            Kws /= aws;
        }
        else {
            Jws = 0;
            Kws = 0;
            JwsS = 0;
            KwsS = 0;
        }
        
        if (ans > 0.0) {
            //            Jns /= ans;
            //            Kns /= ans;
        }
        else {
            Jns = 0;
            Kns = 0;
            JnsS = 0;
            KnsS = 0;
        }
        
        if (lwns > 0.0) {
            KNwns /= lwns;
            KGwns /= lwns;
            // efawns /= lwns;
        }
        else {
            KNwns = 0;
            KGwns = 0;
            // efawns = 1.0;
        }
        
        
        for (int m = 0; m < 3; m++) {
            if (wp_volume > 0.0) {vw(m) /= wp_volume;   } else { vw(m) = 0; }
            if (nwp_volume > 0.0) {vn(m) /= nwp_volume; } else { vn(m) = 0; }
            if (awn > 0.0) {vnw(m) /= awn;} else { vnw(m) = 0; }
            if (aws > 0.0) {vws(m) /= aws;} else { vws(m) = 0; }
            if (ans > 0.0) {vns(m) /= ans;} else { vns(m) = 0; }
            if (lwns > 0.0) {vnws(m) /= lwns;} else { vnws(m) = 0; }
        }
        
        for (int m = 0; m < 6; m++ ) {
            if (awn > 0.0) {Gwn(m) /= awn; } else { Gwn(m) = 0; }
            if (ans > 0.0) {Gns(m) /= ans; } else { Gns(m) = 0; }
            if (aws > 0.0) {Gws(m) /= aws; } else { Gws(m) = 0; }
            if (lwns > 0.0){Gwns(m) /= lwns; } else { Gwns(m) = 0; }
        }
        
        //        wwn(0) *= iVol_nodal;
        if (awn > 0) { wwn(0) /= awn; }  else{ wwn(0) = 0;}
        //        wwn(1) *= iVol_nodal;
        if (awn > 0) { wwn(1) /= awn; }  else{ wwn(1) = 0;}
        //        wwn(2) *= iVol_nodal;
        if (awn > 0) { wwn(2) /= awn; }  else{ wwn(2) = 0;}
        
        //        nwOwnO(0) *= iVol_nodal;
        //        nwOwnO(1) *= iVol_nodal;
        //        nwOwnO(2) *= iVol_nodal;
        //
        //        nwOwsO(0) *= iVol_nodal;
        //        nwOwsO(1) *= iVol_nodal;
        //        nwOwsO(2) *= iVol_nodal;
        
        // vnwdnw *= iVol_nodal;
        //        Eta *= iVol_nodal;
        
        // 1^{wn}
        //        awn *= iVol_nodal;
        //        // 1^{ws}
        //        aws *= iVol_nodal;
        //        // 1^{ns}
        //        ans *= iVol_nodal;
        //        // 1^{nws}
        //        lwns *= iVol_nodal;
        //        // 1^w
        //        aw *= iVol_nodal;
        // 1^n
        // an *= iVol_nodal;
        // printf("vnwdnw=%f\n",vnwdnw);
        
        nw(0) = nw(1) = nw(2) = 0;
        
        
        //double Ca = fabs(vn(2) *average_mu)/(5.796*alpha*1.0);
        double Ca = fabs(vn(2)*muA + vw(2)*muB)/(5.796*alpha*1.0);
        double Re = fabs(vn(2) *1.0 * sphere_diameter)/(average_mu);
        
        double sat_w = 1.0 - nwp_volume/(nwp_volume+wp_volume);
        //        an = nwp_volume_global/(nwp_volume_global+wp_volume_global);
        //        aw = wp_volume_global/(nwp_volume_global+wp_volume_global);
        
        // double poros = 1.0 - porosity;
        // efawns = std::acos(efawns)*180/PI;
        //               efawns = 180 - efawns;
        //               efawns = (efawns >= 180) ? efawns -180 : efawns;
        fprintf(TIMELOG_local,"%i %-15.15E ",timestep, porosity);                        // porosity
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",sat_w,pw,pn);                     // saturation and pressure
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",awn,ans,aws);                       // interfacial areas
        fprintf(TIMELOG_local,"%-15.15E %-15.15E ",Jwn, Kwn);                               // curvature of wn interface
        fprintf(TIMELOG_local,"%-15.15E ",lwns);                                        // common curve length
        fprintf(TIMELOG_local,"%-15.15E ",efawns);                                      // average contact angle
        fprintf(TIMELOG_local,"%-15.15E %-15.15E ",KNwns, KGwns);                           // curvature of wn interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vw(0),vw(1),vw(2));              // average velocity of w phase
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vn(0),vn(1),vn(2));              // average velocity of n phase
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vnw(0),vnw(1),vnw(2));           // velocity of wn interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vnws(0),vnws(1),vnws(2));        // velocity of wn interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwn(0),Gwn(1),Gwn(2),Gwn(3),Gwn(4),Gwn(5));    // orientation of wn interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gws(0),Gws(1),Gws(2),Gws(3),Gws(4),Gws(5));    // orientation of ns interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gns(0),Gns(1),Gns(2),Gns(3),Gns(4),Gns(5));    // orientation of ws interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E ", pwn, euler, aw, an);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwns(0),Gwns(1),Gwns(2),Gwns(3),Gwns(4),Gwns(5)); // orientation of wns interface
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                vws(0),vws(1),vws(2),vns(0),vns(1),vns(2)); // velocity ws and ns
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",
                nwOwnO2(0),nwOwnO2(1),nwOwnO2(2)); // velocity ws and ns
        fprintf(TIMELOG_local,"%-15.15E ",
                Eta);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",
                wwn(0),wwn(1),wwn(2));
//        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E\n",Jwn2, Kwn2, Ca, Re, an1, an2, an3, Eta1, dendt, dendt3);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E\n",Ca, Re);
        //        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",trawn, trJwn, trRwn);               // Trimmed curvature
        //        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",wwndnw, wwnsdnwn, Jwnwwndnw);       // kinematic quantities
        //        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E ",euler, Kn, Jn, An);           // minkowski measures
        //        fprintf(TIMELOG_local,"%-15.15E\n",swn);
        
        
        //        const char header[] = "time porosity "                  // Timestep, porosity, Change in Surface Energy
        //                  "s^w p^w p^n eps^{wn} eps^{ns} eps^{ws} J_n^{wn} K_n^{wn} l^{wns} c^{wns} kappa_N^{wns} kappa_G^{wns} " // Scalar averages
        //                  "v_x^w v_y^w v_z^w v_x^n v_y^n v_z^n v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} " // Velocity averages
        //                  "G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} "
        //                  "G_{xx}^{ws} G_{yy}^{ws} G_{zz}^{ws} G_{xy}^{ws} G_{xz}^{ws} G_{yz}^{ws} "
        //                  "G_{xx}^{ns} G_{yy}^{ns} G_{zz}^{ns} G_{xy}^{ns} G_{xz}^{ns} G_{yz}^{ns} "
        //                  "p_n^{wn} X^n eps^w eps^n "
        //                  "G_{xx}^{wns} G_{yy}^{wns} G_{zz}^{wns} G_{xy}^{wns} G_{xz}^{wns} G_{yz}^{wns} "
        //                  "v_x^{ws} v_y^{ws} v_z^{ws} v_x^{ns} v_y^{ns} v_z^{ns} nw_x nw_y nw_z vnwdnw wwnx wwny wwnz Jwn2* Kwn2*";
        fflush(TIMELOG_local);
    }
}


void TwoPhase::PrintComponents(int timestep)
{
    if (Dm->rank()==0){
        printf("PRINT %i COMPONENT AVERAGES: time = %i \n",(int)ComponentAverages_NWP.size(1),timestep);
        for (int b=0; b<NumberComponents_NWP; b++){
            //if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
            fprintf(NWPLOG,"%i ",timestep);
            if ( Label_NWP_map.empty() ) {
                // print index
                fprintf(NWPLOG,"%i ",b);
            } else {
                // print final global id
                fprintf(NWPLOG,"%i ",Label_NWP_map[b]);
            }
            fprintf(NWPLOG,"%-15.15E\n",ComponentAverages_NWP(VOL,b));
            //			fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(TRIMVOL,b));
            /*fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(PRS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(AWN,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(AWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(JWN,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(KWN,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(LWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(KNWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(KGWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNWX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNWY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNWZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNWSX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNWSY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNWSZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNXX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNYY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNZZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNXY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNXZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNYZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GXX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GYY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GZZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GXY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GXZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GYZ,b));
            //fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(TRAWN,b));
            //fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(TRJWN,b));
            //fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(INTCURV,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(EULER,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CMX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CMY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CMZ,b));
            fprintf(NWPLOG,"%-15.15E\n",ComponentAverages_NWP(VSQ,b));*/
            //				fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(NVERT,b));
            //			fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(NSIDE,b));
            //		fprintf(NWPLOG,"%-15.15E\n",ComponentAverages_NWP(NFACE,b));
            //}
        }
        fflush(NWPLOG);
        
        for (int b=0; b<NumberComponents_WP; b++){
            if (ComponentAverages_WP(TRIMVOL,b) > 0.0 ){
                fprintf(WPLOG,"%i ",timestep);
                fprintf(WPLOG,"%i ",b);
                fprintf(WPLOG,"%-15.15E\n",ComponentAverages_WP(VOL,b));
                //			fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(TRIMVOL,b));
                /*fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(PRS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(AWN,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(AWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(JWN,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(KWN,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(LWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(KNWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(KGWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VNWX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VNWY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VNWZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VNWSX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VNWSY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VNWSZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNXX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNYY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNZZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNXY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNXZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNYZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GXX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GYY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GZZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GXY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GXZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GYZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CMX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CMY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CMZ,b));
                fprintf(WPLOG,"%-15.15E\n",ComponentAverages_WP(VSQ,b));*/
                //fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(TRAWN,b));
                //fprintf(WPLOG,"%-15.15E\n",ComponentAverages_WP(TRJWN,b));
                
            }
        }
        fflush(WPLOG);
    }
}


inline int TwoPhase::GetCubeLabel(int i, int j, int k, IntArray &BlobLabel)
{
    int label;
    label=BlobLabel(i,j,k);
    label=max(label,BlobLabel(i+1,j,k));
    label=max(label,BlobLabel(i,j+1,k));
    label=max(label,BlobLabel(i+1,j+1,k));
    label=max(label,BlobLabel(i,j,k+1));
    label=max(label,BlobLabel(i+1,j,k+1));
    label=max(label,BlobLabel(i,j+1,k+1));
    label=max(label,BlobLabel(i+1,j+1,k+1));
    
    return label;
}

void TwoPhase::SortBlobs()
{
    //printf("Sorting the blobs based on volume \n");
    //printf("-----------------------------------------------\n");
    int TempLabel,a,aa,bb,i,j,k,idx;
    double TempValue;
    //.......................................................................
    // Sort NWP components by volume
    //.......................................................................
    IntArray OldLabel(NumberComponents_NWP);
    for (a=0; a<NumberComponents_NWP; a++)	OldLabel(a) = a;
    // Sort the blob averages based on volume
    for (aa=0; aa<NumberComponents_NWP-1; aa++){
        for ( bb=aa+1; bb<NumberComponents_NWP; bb++){
            if (ComponentAverages_NWP(VOL,aa) < ComponentAverages_NWP(VOL,bb)){
                // Exchange location of blobs aa and bb
                //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
                // switch the label
                TempLabel = OldLabel(bb);
                OldLabel(bb) = OldLabel(aa);
                OldLabel(aa) = TempLabel;
                // switch the averages
                for (idx=0; idx<BLOB_AVG_COUNT; idx++){
                    TempValue = ComponentAverages_NWP(idx,bb);
                    ComponentAverages_NWP(idx,bb) = ComponentAverages_NWP(idx,aa);
                    ComponentAverages_NWP(idx,aa) = TempValue;
                }
            }
        }
    }
    IntArray NewLabel(NumberComponents_NWP);
    for (aa=0; aa<NumberComponents_NWP; aa++){
        // Match the new label for original blob aa
        bb=0;
        while (OldLabel(bb) != aa)	bb++;
        NewLabel(aa) = bb;
    }
    // Re-label the blob ID
    //	printf("Re-labeling the blobs, now indexed by volume \n");
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                if (Label_NWP(i,j,k) > -1){
                    TempLabel = NewLabel(Label_NWP(i,j,k));
                    Label_NWP(i,j,k) = TempLabel;
                }
            }
        }
    }
    //.......................................................................
    // Sort WP components by volume
    //.......................................................................
    OldLabel.resize(NumberComponents_WP);
    for (a=0; a<NumberComponents_WP; a++)	OldLabel(a) = a;
    // Sort the blob averages based on volume
    for (aa=0; aa<NumberComponents_WP-1; aa++){
        for ( bb=aa+1; bb<NumberComponents_WP; bb++){
            if (ComponentAverages_WP(VOL,aa) < ComponentAverages_WP(VOL,bb)){
                // Exchange location of blobs aa and bb
                //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
                // switch the label
                TempLabel = OldLabel(bb);
                OldLabel(bb) = OldLabel(aa);
                OldLabel(aa) = TempLabel;
                // switch the averages
                for (idx=0; idx<BLOB_AVG_COUNT; idx++){
                    TempValue = ComponentAverages_WP(idx,bb);
                    ComponentAverages_WP(idx,bb) = ComponentAverages_WP(idx,aa);
                    ComponentAverages_WP(idx,aa) = TempValue;
                }
            }
        }
    }
    NewLabel.resize(NumberComponents_WP);
    for (aa=0; aa<NumberComponents_WP; aa++){
        // Match the new label for original blob aa
        bb=0;
        while (OldLabel(bb) != aa)	bb++;
        NewLabel(aa) = bb;
    }
    // Re-label the blob ID
    //	printf("Re-labeling the blobs, now indexed by volume \n");
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                if (Label_WP(i,j,k) > -1){
                    TempLabel = NewLabel(Label_WP(i,j,k));
                    Label_WP(i,j,k) = TempLabel;
                }
            }
        }
    }  
}

void TwoPhase::AggregateLabels( const std::string& filename )
{
	
	int nx = Dm->Nx;
	int ny = Dm->Ny;
	int nz = Dm->Nz;
		
	//printf("aggregate labels: local size=%i, global size = %i",local_size, full_size);
	// assign the ID from the phase indicator field
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				int n = k*nx*ny+j*nx+i;
				signed char local_id_val = Dm->id[n]; 
				if (local_id_val > 0){
					double value = TempPhase(i,j,k);
					if (value > 0.0)	local_id_val = 1;
					else 				local_id_val = 2;
				}
				 Dm->id[n] = local_id_val;
			}
		}
	}
	MPI_Barrier(Dm->Comm);

	Dm->AggregateLabels( filename );

}


