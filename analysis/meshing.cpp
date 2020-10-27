#include "analysis/TwoPhase.h"

#include <stdint.h>
#include <limits.h>


void TwoPhase::nwComputeQuantities(nwMyMesh & m_nw) {
    
    double iVol_global = 1.0/double(Volume);
    double iVol_nodal = 1.0/double(nodeVolume);
    
    /* First steps are to delete any duplicate faces and vertices */
    int beforeVN, beforeFN, nwAfterVN, nwAfterFN;
    beforeVN = nwAfterVN = beforeFN = nwAfterFN = 0;
    beforeVN = m_nw.VN(); beforeFN = m_nw.FN();


    vcg::tri::Clean<nwMyMesh>::RemoveDuplicateFace(m_nw);
    vcg::tri::Clean<nwMyMesh>::RemoveDuplicateVertex(m_nw,false);

    size_t count = 0;

    /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
    int degenface = vcg::tri::Clean<nwMyMesh>::RemoveDegenerateFace(m_nw);
    vcg::tri::UpdateTopology<nwMyMesh>::FaceFace(m_nw);
    int degenvertex = vcg::tri::Clean<nwMyMesh>::RemoveDegenerateVertex(m_nw);
    vcg::tri::UpdateTopology<nwMyMesh>::VertexFace(m_nw);

    /* FINISH REPAIRING STEPS */
    // int splitNonManifoldVertex = vcg::tri::Clean<nwMyMesh>::SplitNonManifoldVertex(m_nw, 5);

    /* Need to remove any non 2-manifold faces */
    int removeNonManifoldFace = vcg::tri::Clean<nwMyMesh>::RemoveNonManifoldFace(m_nw);

    std::pair<int,int> nw_CCs;
    nw_CCs =  vcg::tri::Clean<nwMyMesh>::RemoveSmallConnectedComponentsDiameter(m_nw, 0.0 /* seems to be a good threshold based upon visual experiments in MeshLab */);
    nw_CCs.first = nw_CCs.first - nw_CCs.second;

    //  std::cout << Dm->rank() << ") CC=" << nw_CCs.first << std::endl;

    int removeUnreferencedVertex = vcg::tri::Clean<nwMyMesh>::RemoveUnreferencedVertex(m_nw, true);

    removeNonManifoldFace = vcg::tri::Clean<nwMyMesh>::RemoveNonManifoldFace(m_nw);

    degenvertex = vcg::tri::Clean<nwMyMesh>::RemoveDegenerateVertex(m_nw);

    // Final check for Euler
    vcg::tri::Clean<nwMyMesh>::RemoveDuplicateFace(m_nw);
    vcg::tri::Clean<nwMyMesh>::RemoveDuplicateVertex(m_nw,false);

    vcg::tri::Allocator<nwMyMesh>::CompactEveryVector(m_nw);

    vcg::tri::UpdateCurvatureFitting<nwMyMesh>::computeCurvature(m_nw);

    double area = vcg::tri::Stat<nwMyMesh>::ComputeMeshArea(m_nw);
    awn = area;
    double area_global = 0.0;
    MPI_Allreduce(&area,&area_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    
    
//    double volen = -vcg::tri::Stat<nwMyMesh>::ComputeMeshVolume(m_nw);
//    double volen_global = 0;
//   //
//
//    MPI_Allreduce(&volen,&volen_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
//
//    volen_global *= iVol_global;
//    volen *= iVol_nodal;
//
//     std::cout << "volen=" << volen << " volen_global=" << volen_global << std::endl;
//
//    an = volen;
//    an_global = volen_global;
    
    
//
//
//    volen_global *= iVol_global;
//    volen *= iVol_nodal;
//
//    an_global = volen_global;
//    an = volen;
   

    /* mean and Gaussian curvature calculations */
    double k1 = 0.0;
    double k2 = 0.0;
    double k1_global = 0;
    double k2_global = 0;

    int k1nan = 0;
    int k2nan = 0;
    double tmp = 0;
    for(auto vi = m_nw.vert.begin(); vi != m_nw.vert.end(); ++vi) {
        tmp = (*vi).K1();
        if (isnan(tmp)) {k1nan+=1;}
        else k1 += tmp;
        tmp = (*vi).K2();
        if (isnan(tmp)) {k2nan+=1;}
        else k2 += tmp;
    }
    // std::cout << "(" << Dm->rank() << ") nw: k1nan=" << k1nan << " k2nan=" << k2nan << " k1,k2=" << k1 << " " << k2 <<  std::endl;

    double hval = 0;
    double gval = 0;

    hval = -0.5*(k1+k2);
    gval = k1*k2;

    // avg k1 and k2 on a process
    //(1a)
    if (m_nw.VN() > 0) k1 = 2*hval/double(m_nw.VN()); else k1 = 0;
    if (m_nw.VN() > 0) k2 = gval/double(m_nw.VN())/double(m_nw.VN()); else k2 = 0;

    // (1b)
    //    k1 = hval;
    //    k2 = gval;

    Jwn = k1;
    Kwn = k2;

    // std::cout << Dm->rank() << ") Jwn=" << Jwn << std::endl;

    JwnW = -k1;
    KwnW = -k2;

    awn_global = area_global;// /=double(npx*npy*npz);

    MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);


    k1_global/=double(npx*npy*npz);
    k2_global/=double(npx*npy*npz);

    Jwn_global = k1_global;
    Kwn_global = k2_global;

    JwnW_global = -k1_global;
    KwnW_global = -k2_global;
    

    

    
}







void TwoPhase::nsComputeQuantities(nsMyMesh & m_ns) {
    
    /* First steps are to delete any duplicate faces and vertices */
    int beforeVN, beforeFN, nsAfterVN, nsAfterFN; beforeVN = nsAfterVN = beforeFN = nsAfterFN = 0;  beforeVN = m_ns.VN(); beforeFN = m_ns.FN();
    vcg::tri::Clean<nsMyMesh>::RemoveDuplicateFace(m_ns);
    vcg::tri::Clean<nsMyMesh>::RemoveDuplicateVertex(m_ns,false);

    // cout << "beforeVN=" << beforeVN << endl;

   vcg::tri::Clean<nsMyMesh>::RemoveDuplicateFace(m_ns);
   vcg::tri::Clean<nsMyMesh>::RemoveDuplicateVertex(m_ns,false);

   size_t count = 0;

   /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
   int degenface = vcg::tri::Clean<nsMyMesh>::RemoveDegenerateFace(m_ns);
   vcg::tri::UpdateTopology<nsMyMesh>::FaceFace(m_ns);
   int degenvertex = vcg::tri::Clean<nsMyMesh>::RemoveDegenerateVertex(m_ns);
   vcg::tri::UpdateTopology<nsMyMesh>::VertexFace(m_ns);

   /* FINISH REPAIRING STEPS */
   // int splitNonManifoldVertex = vcg::tri::Clean<nsMyMesh>::SplitNonManifoldVertex(m_ns, 5);

   /* Need to remove any non 2-manifold faces */
   int removeNonManifoldFace = vcg::tri::Clean<nsMyMesh>::RemoveNonManifoldFace(m_ns);

   std::pair<int,int> ns_CCs;
   ns_CCs =  vcg::tri::Clean<nsMyMesh>::RemoveSmallConnectedComponentsDiameter(m_ns, 0.0 /* seems to be a good threshold based upon visual experiments in MeshLab */);
   ns_CCs.first = ns_CCs.first - ns_CCs.second;

   //  std::cout << Dm->rank() << ") CC=" << ns_CCs.first << std::endl;

   int removeUnreferencedVertex = vcg::tri::Clean<nsMyMesh>::RemoveUnreferencedVertex(m_ns, true);

   removeNonManifoldFace = vcg::tri::Clean<nsMyMesh>::RemoveNonManifoldFace(m_ns);

   degenvertex = vcg::tri::Clean<nsMyMesh>::RemoveDegenerateVertex(m_ns);

   // Final check for Euler
   vcg::tri::Clean<nsMyMesh>::RemoveDuplicateFace(m_ns);
   vcg::tri::Clean<nsMyMesh>::RemoveDuplicateVertex(m_ns,false);

   vcg::tri::Allocator<nsMyMesh>::CompactEveryVector(m_ns);

    /* Hope that the damn mesh is finally OK... */
    vcg::tri::UpdateCurvatureFitting<nsMyMesh>::computeCurvature(m_ns);


    double area = vcg::tri::Stat<nsMyMesh>::ComputeMeshArea(m_ns);
    ans = area;
    double area_global = 0.0;
    MPI_Allreduce(&area,&area_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    /* mean and Gaussian curvature calculations */
    double k1 = 0.0;
    double k2 = 0.0;
    double k1_global = 0;
    double k2_global = 0;

    int k1nan = 0;
    int k2nan = 0;
    double tmp = 0;
    for(auto vi = m_ns.vert.begin(); vi != m_ns.vert.end(); ++vi) {
        tmp = (*vi).K1();
        if (isnan(tmp)) {k1nan+=1;}
        else k1 += tmp;
        tmp = (*vi).K2();
        if (isnan(tmp)) {k2nan+=1;}
        else k2 += tmp;
    }

    // std::cout << "(" << Dm->rank() << ") NS: k1nan=" << k1nan << " k2nan=" << k2nan << " k1,k2=" << k1 << " " << k2 <<  std::endl;

    double hval = 0;
    double gval = 0;

    hval = -0.5*(k1+k2);
    gval = k1*k2;

    // avg k1 and k2 on a process
    if (m_ns.VN() > 0) k1 = 2*hval/double(m_ns.VN()); else k1 = 0;
    if (m_ns.VN() > 0) k2 = gval/double(m_ns.VN())/double(m_ns.VN()); else k2 = 0;

    Jns = k1;
    Kns = k2;

    JnsS = -k1;
    KnsS = -k2;

    ans_global = area_global;// /=double(npx*npy*npz);

    MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    k1_global/=double(npx*npy*npz);
    k2_global/=double(npx*npy*npz);

    Jns_global = k1_global;
    Kns_global = k2_global;

    JnsS_global = -k1_global;
    KnsS_global = -k2_global;


    
    
    
}





void TwoPhase::wsComputeQuantities(wsMyMesh & m_ws) {
    double iVol_global = 1.0/double(Volume);
    
    /* First steps are to delete any duplicate faces and vertices */
    int beforeVN, beforeFN, wsAfterVN, wsAfterFN;
    beforeVN = wsAfterVN = beforeFN = wsAfterFN = 0;
    beforeVN = m_ws.VN(); beforeFN = m_ws.FN();


  vcg::tri::Clean<wsMyMesh>::RemoveDuplicateFace(m_ws);
    vcg::tri::Clean<wsMyMesh>::RemoveDuplicateVertex(m_ws,false);

    size_t count = 0;

    /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
    int degenface = vcg::tri::Clean<wsMyMesh>::RemoveDegenerateFace(m_ws);
    vcg::tri::UpdateTopology<wsMyMesh>::FaceFace(m_ws);
    int degenvertex = vcg::tri::Clean<wsMyMesh>::RemoveDegenerateVertex(m_ws);
    vcg::tri::UpdateTopology<wsMyMesh>::VertexFace(m_ws);

    /* FINISH REPAIRING STEPS */
    // int splitNonManifoldVertex = vcg::tri::Clean<wsMyMesh>::SplitNonManifoldVertex(m_ws, 5);

    /* Need to remove any non 2-manifold faces */
    int removeNonManifoldFace = vcg::tri::Clean<wsMyMesh>::RemoveNonManifoldFace(m_ws);

    std::pair<int,int> ws_CCs;
    ws_CCs =  vcg::tri::Clean<wsMyMesh>::RemoveSmallConnectedComponentsDiameter(m_ws, 0.0 /* seems to be a good threshold based upon visual experiments in MeshLab */);
    ws_CCs.first = ws_CCs.first - ws_CCs.second;

    //  std::cout << Dm->rank() << ") CC=" << ws_CCs.first << std::endl;

    int removeUnreferencedVertex = vcg::tri::Clean<wsMyMesh>::RemoveUnreferencedVertex(m_ws, true);

    removeNonManifoldFace = vcg::tri::Clean<wsMyMesh>::RemoveNonManifoldFace(m_ws);

    degenvertex = vcg::tri::Clean<wsMyMesh>::RemoveDegenerateVertex(m_ws);

    // Final check for Euler
    vcg::tri::Clean<wsMyMesh>::RemoveDuplicateFace(m_ws);
    vcg::tri::Clean<wsMyMesh>::RemoveDuplicateVertex(m_ws,false);

    vcg::tri::Allocator<wsMyMesh>::CompactEveryVector(m_ws);

    vcg::tri::UpdateCurvatureFitting<wsMyMesh>::computeCurvature(m_ws);



    double area = vcg::tri::Stat<wsMyMesh>::ComputeMeshArea(m_ws);
    awn = area;
    double area_global = 0.0;
    MPI_Allreduce(&area,&area_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    /* mean and Gaussian curvature calculations */
    double k1 = 0.0;
    double k2 = 0.0;
    double k1_global = 0;
    double k2_global = 0;

    int k1nan = 0;
    int k2nan = 0;
    double tmp = 0;
    for(auto vi = m_ws.vert.begin(); vi != m_ws.vert.end(); ++vi) {
        tmp = (*vi).K1();
        if (isnan(tmp)) {k1nan+=1;}
        else k1 += tmp;
        tmp = (*vi).K2();
        if (isnan(tmp)) {k2nan+=1;}
        else k2 += tmp;
    }
    // std::cout << "(" << Dm->rank() << ") ws: k1nan=" << k1nan << " k2nan=" << k2nan << " k1,k2=" << k1 << " " << k2 <<  std::endl;

    double hval = 0;
    double gval = 0;

    hval = -0.5*(k1+k2);
    gval = k1*k2;

    // avg k1 and k2 on a process
    //(1a)
    if (m_ws.VN() > 0) k1 = 2*hval/double(m_ws.VN()); else k1 = 0;
    if (m_ws.VN() > 0) k2 = gval/double(m_ws.VN())/double(m_ws.VN()); else k2 = 0;

    // (1b)
    //    k1 = hval;
    //    k2 = gval;

    Jws = k1;
    Kws = k2;

    // std::cout << Dm->rank() << ") Jwn=" << Jwn << std::endl;

    Jws = -k1;
    Kws = -k2;

    aws_global = area_global;// /=double(npx*npy*npz);

    MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);


    k1_global/=double(npx*npy*npz);
    k2_global/=double(npx*npy*npz);

    Jws_global = k1_global;
    Kws_global = k2_global;

    Jws_global = -k1_global;
    Kws_global = -k2_global;
    
}




void TwoPhase::sComputeQuantities(sMyMesh & m_s) {}



void TwoPhase::ComputeAccurateNWInterfaceCurvatures(nwMyMesh & m_nw) {}




void TwoPhase::wComputeQuantities(wMyMesh & m_w) {}


void TwoPhase::nComputeQuantities(nMyMesh & m_n) {}
