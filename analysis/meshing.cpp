#include "analysis/TwoPhase.h"



#include <stdint.h>
#include <limits.h>

#include <string>


void TwoPhase::nwComputeQuantities(nwMyMesh & m_nw) {
    
    double iVol_global = 1.0/double(Volume);
    double iVol_nodal = 1.0/double(nodeVolume);
    
    /* First steps are to delete any duplicate faces and vertices */
    int beforeVN, beforeFN, nwAfterVN, nwAfterFN;
    beforeVN = nwAfterVN = beforeFN = nwAfterFN = 0;
    beforeVN = m_nw.VN(); beforeFN = m_nw.FN();


//    vcg::tri::Clean<nwMyMesh>::RemoveDuplicateFace(m_nw);
//    vcg::tri::Clean<nwMyMesh>::RemoveDuplicateVertex(m_nw,false);

    size_t count = 0;

    /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
    int degenface = vcg::tri::Clean<nwMyMesh>::RemoveDegenerateFace(m_nw);
    vcg::tri::UpdateTopology<nwMyMesh>::FaceFace(m_nw);
    int degenvertex = vcg::tri::Clean<nwMyMesh>::RemoveDegenerateVertex(m_nw);
    vcg::tri::UpdateTopology<nwMyMesh>::VertexFace(m_nw);

    /* FINISH REPAIRING STEPS */
     int splitNonManifoldVertex = vcg::tri::Clean<nwMyMesh>::SplitNonManifoldVertex(m_nw, 5);

    // Need to remove any non 2-manifold faces
    int removeNonManifoldFace = vcg::tri::Clean<nwMyMesh>::RemoveNonManifoldFace(m_nw);

    std::pair<int,int> nw_CCs;
    nw_CCs =  vcg::tri::Clean<nwMyMesh>::RemoveSmallConnectedComponentsDiameter(m_nw, 0.0 /* seems to be a good threshold based upon visual experiments in MeshLab */);
    nw_CCs.first = nw_CCs.first - nw_CCs.second;

    //  std::cout << Dm->rank() << ") CC=" << nw_CCs.first << std::endl;

    int removeUnreferencedVertex = vcg::tri::Clean<nwMyMesh>::RemoveUnreferencedVertex(m_nw, true);

//    removeNonManifoldFace = vcg::tri::Clean<nwMyMesh>::RemoveNonManifoldFace(m_nw);

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

    hval = 0.5*(k1+k2);
    gval = k1*k2;
    
    if (m_nw.VN() > 0) k1 = 2*hval/double(m_nw.VN()); else k1 = 0;
    if (m_nw.VN() > 0) k2 = gval/double(m_nw.VN())/double(m_nw.VN()); else k2 = 0;


    Jwn = k1;
    Kwn = k2;

    awn_global = area_global;// /=double(npx*npy*npz);

    MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);


    k1_global/=double(npx*npy*npz);
    k2_global/=double(npx*npy*npz);

    Jwn_global = k1_global;
    Kwn_global = k2_global;
    
//    std::cout << "Jwn=" << Jwn << std::endl;
    vcg::tri::UpdateNormal<nwMyMesh>::PerFaceNormalized(m_nw);
    

    int norm_count = 0;
    double k1sum, k2sum;
    double geo_area = 0;
    k1sum = k2sum = 0.0;
    int vertex_sum = 0;
    nwMyMesh::VertexPointer ivp[3];

    nwOwnO2(0) = nwOwnO2(1) = nwOwnO2(2) = 0;
    for(nwMyMesh::FaceIterator fi=m_nw.face.begin();fi!=m_nw.face.end();++fi) if(!(*fi).IsD())
    {
        double nx = (*fi).N()[0];
        double ny = (*fi).N()[1];
        double nz = (*fi).N()[2];

        double normalMag = nx*nx + ny*ny + nz*nz;
        if ((normalMag - 1.0) > 1E-14) { count++; }

        double geo_area = vcg::DoubleArea(*fi)/2.0;

        nwOwnO2(0) += nx*geo_area;
        nwOwnO2(1) += ny*geo_area;
        nwOwnO2(2) += nz*geo_area;
    }

    MPI_Allreduce(&nwOwnO2(0),&nwOwnO_global2(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    // Dont scale the normals since the averaging domain will be larger
    nwOwnO_global2(0) *= iVol_global;
    nwOwnO_global2(1) *= iVol_global;
    nwOwnO_global2(2) *= iVol_global;
    
    nwOwnO2(0) *= iVol_nodal;
    nwOwnO2(1) *= iVol_nodal;
    nwOwnO2(2) *= iVol_nodal;
    
   
    
  
    
    
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

    ans_global = area_global;// /=double(npx*npy*npz);

    MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    k1_global/=double(npx*npy*npz);
    k2_global/=double(npx*npy*npz);

    Jns_global = k1_global;
    Kns_global = k2_global;


    nwOnsO2(0) = nwOnsO2(1) = nwOnsO2(2) = 0;
    for(nsMyMesh::FaceIterator fi=m_ns.face.begin();fi!=m_ns.face.end();++fi) if(!(*fi).IsD())
    {
        double nx = (*fi).N()[0];
        double ny = (*fi).N()[1];
        double nz = (*fi).N()[2];

        double normalMag = nx*nx + ny*ny + nz*nz;
        if ((normalMag - 1.0) > 1E-14) { count++; }

        double geo_area = vcg::DoubleArea(*fi)/2.0;

        nwOnsO2(0) += nx*geo_area;
        nwOnsO2(1) += ny*geo_area;
        nwOnsO2(2) += nz*geo_area;
    }
    
  
    
    
}



void TwoPhase::nwsComputeQuantities(nwsMyMesh & m_nws) {
    
    /* First steps are to delete any duplicate faces and vertices */
      int beforeVN, beforeFN, nsAfterVN, nsAfterFN; beforeVN = nsAfterVN = beforeFN = nsAfterFN = 0;  beforeVN = m_nws.VN(); beforeFN = m_nws.FN();
      vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateFace(m_nws);
      vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateVertex(m_nws,false);

      // cout << "beforeVN=" << beforeVN << endl;

     vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateFace(m_nws);
     vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateVertex(m_nws,false);

     size_t count = 0;

     /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
     int degenface = vcg::tri::Clean<nwsMyMesh>::RemoveDegenerateFace(m_nws);
     vcg::tri::UpdateTopology<nwsMyMesh>::FaceFace(m_nws);
     int degenvertex = vcg::tri::Clean<nwsMyMesh>::RemoveDegenerateVertex(m_nws);
     vcg::tri::UpdateTopology<nwsMyMesh>::VertexFace(m_nws);

     /* FINISH REPAIRING STEPS */
     // int splitNonManifoldVertex = vcg::tri::Clean<nwsMyMesh>::SplitNonManifoldVertex(m_nws, 5);

     /* Need to remove any non 2-manifold faces */
     int removeNonManifoldFace = vcg::tri::Clean<nwsMyMesh>::RemoveNonManifoldFace(m_nws);

     std::pair<int,int> ns_CCs;
     ns_CCs =  vcg::tri::Clean<nwsMyMesh>::RemoveSmallConnectedComponentsDiameter(m_nws, 0.0 /* seems to be a good threshold based upon visual experiments in MeshLab */);
     ns_CCs.first = ns_CCs.first - ns_CCs.second;

     //  std::cout << Dm->rank() << ") CC=" << ns_CCs.first << std::endl;

     int removeUnreferencedVertex = vcg::tri::Clean<nwsMyMesh>::RemoveUnreferencedVertex(m_nws, true);

     removeNonManifoldFace = vcg::tri::Clean<nwsMyMesh>::RemoveNonManifoldFace(m_nws);

     degenvertex = vcg::tri::Clean<nwsMyMesh>::RemoveDegenerateVertex(m_nws);

     // Final check for Euler
     vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateFace(m_nws);
     vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateVertex(m_nws,false);

     vcg::tri::Allocator<nwsMyMesh>::CompactEveryVector(m_nws);

    vcg::tri::UpdateTopology<nwsMyMesh>::AllocateEdge(m_nws);
   
    
    
  
    double RadAngle = 0;
    double DegAngle = 0;
    size_t DegCount = 0;
    for(nwsMyMesh::FaceIterator fi=m_nws.face.begin();fi!=m_nws.face.end();++fi) if(!(*fi).IsD())
    {
        for(int z=0;z<(*fi).VN();++z)
        {
            RadAngle = std::abs(vcg::face::DihedralAngleRad(*fi,z)); // adjust incorrect sign depending on phase affinity later
            RadAngle *= (180 /3.14159265359);
            if (RadAngle > 1E-4) {
              //  std::cout << "z=" << z << " DegreeAngle=" << RadAngle << std::endl;

                DegAngle += RadAngle;
                DegCount++;
            }
           // }
        }
        
    }
    
    if (DegCount > 0) {
        DegAngle /= double(DegCount);
        efawns += DegAngle;
        efawns_count++;
    } else {
        efawns += 0.0;
        efawns_count += 0;
    }
    
    
        
   
//    for(nwsMyMesh::EdgeIterator ei=m_nws.edge.begin(); ei!=m_nws.edge.end();++ei) if(!(*ei).IsD())
//    {
//        double nx1 = (*fi).N()[0];
//        double ny1 = (*fi).N()[1];
//        double nz1 = (*fi).N()[2];
//
//        double normalMag = nx1*nx1 + ny1*ny1 + nz1*nz1;
//        if ((normalMag - 1.0) > 1E-14) { count++; }
//        nx1 /= normalMag;
//        ny1 /= normalMag;
//        nz1 /= normalMag;
//
//    }
//
//    nwOwsO2(0) = nwOwsO2(1) = nwOwsO2(2) = 0;
//    for(nwsMyMesh::FaceIterator fi=m_nws.face.begin();fi!=m_nws.face.end();++fi) if(!(*fi).IsD())
//    {
//        double nx1 = (*fi).N()[0];
//        double ny1 = (*fi).N()[1];
//        double nz1 = (*fi).N()[2];
//
//        double normalMag = nx1*nx1 + ny1*ny1 + nz1*nz1;
//        if ((normalMag - 1.0) > 1E-14) { count++; }
//        nx1 /= normalMag;
//        ny1 /= normalMag;
//        nz1 /= normalMag;
//
//    }
  
    
    
    
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
    aws = area;
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

    double hval = 0;
    double gval = 0;

    hval = -0.5*(k1+k2);
    gval = k1*k2;
    
    if (m_ws.VN() > 0) k1 = 2*hval/double(m_ws.VN()); else k1 = 0;
    if (m_ws.VN() > 0) k2 = gval/double(m_ws.VN())/double(m_ws.VN()); else k2 = 0;

    Jws = k1;
    Kws = k2;
    
    aws_global = area_global;// /=double(npx*npy*npz);

    MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    k1_global/=double(npx*npy*npz);
    k2_global/=double(npx*npy*npz);

    Jws_global = k1_global;
    Kws_global = k2_global;

    
    nwOwsO2(0) = nwOwsO2(1) = nwOwsO2(2) = 0;
      for(wsMyMesh::FaceIterator fi=m_ws.face.begin();fi!=m_ws.face.end();++fi) if(!(*fi).IsD())
      {
          double nx = (*fi).N()[0];
          double ny = (*fi).N()[1];
          double nz = (*fi).N()[2];

          double normalMag = nx*nx + ny*ny + nz*nz;
          if ((normalMag - 1.0) > 1E-14) { count++; }

          double geo_area = vcg::DoubleArea(*fi)/2.0;

          nwOwsO2(0) += nx*geo_area;
          nwOwsO2(1) += ny*geo_area;
          nwOwsO2(2) += nz*geo_area;
      }
    
    
   
    
}




void TwoPhase::sComputeQuantities(sMyMesh & m_s) {
    
    double iVol_global = 1.0/double(Volume);

      /* First steps are to delete any duplicate faces and vertices */
      int beforeVN, beforeFN, sAfterVN, sAfterFN;
      beforeVN = sAfterVN = beforeFN = sAfterFN = 0;
      beforeVN = m_s.VN(); beforeFN = m_s.FN();


      vcg::tri::Clean<sMyMesh>::RemoveDuplicateFace(m_s);
      vcg::tri::Clean<sMyMesh>::RemoveDuplicateVertex(m_s,true);

      size_t count = 0;

      /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
      vcg::tri::Clean<sMyMesh>::RemoveDegenerateFace(m_s);
      vcg::tri::UpdateTopology<sMyMesh>::FaceFace(m_s);
      vcg::tri::Clean<sMyMesh>::RemoveDegenerateVertex(m_s);
      vcg::tri::UpdateTopology<sMyMesh>::VertexFace(m_s);

      /* FINISH REPAIRING STEPS */
      int splitNonManifoldVertex = vcg::tri::Clean<sMyMesh>::SplitNonManifoldVertex(m_s, 5);

      /* Need to remove any non 2-manifold faces */
    //  vcg::tri::Clean<sMyMesh>::RemoveNonManifoldFace(m_s);

      std::pair<int,int> s_CCs;
      s_CCs =  vcg::tri::Clean<sMyMesh>::RemoveSmallConnectedComponentsDiameter(m_s, 3 /* seems to be a good threshold based upon visual experiments in MeshLab */);
      s_CCs.first = s_CCs.first - s_CCs.second;

      //  std::cout << Dm->rank() << ") CC=" << s_CCs.first << std::endl;

      vcg::tri::Clean<sMyMesh>::RemoveUnreferencedVertex(m_s);

      vcg::tri::Clean<sMyMesh>::RemoveNonManifoldFace(m_s);

      vcg::tri::Clean<sMyMesh>::RemoveDegenerateVertex(m_s);

      // Final check for Euler
      vcg::tri::Clean<sMyMesh>::RemoveDuplicateFace(m_s);
      vcg::tri::Clean<sMyMesh>::RemoveDuplicateVertex(m_s);

      vcg::tri::Allocator<sMyMesh>::CompactEveryVector(m_s);

//      vcg::tri::UpdateCurvatureFitting<sMyMesh>::computeCurvature(m_s);

//      double area = vcg::tri::Stat<sMyMesh>::ComputeMeshArea(m_s);
//      awn = area;
//      double area_global = 0.0;
//      MPI_Allreduce(&area,&area_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
//
//      /* mean and Gaussian curvature calculations */
//      double k1 = 0.0;
//      double k2 = 0.0;
//      double k1_global = 0;
//      double k2_global = 0;
//
//      int k1nan = 0;
//      int k2nan = 0;
//      double tmp = 0;
//      for(auto vi = m_s.vert.begin(); vi != m_s.vert.end(); ++vi) {
//          tmp = (*vi).K1();
//          if (isnan(tmp)) {k1nan+=1;}
//          else k1 += tmp;
//          tmp = (*vi).K2();
//          if (isnan(tmp)) {k2nan+=1;}
//          else k2 += tmp;
//      }
//      // std::cout << "(" << Dm->rank() << ") s: k1nan=" << k1nan << " k2nan=" << k2nan << " k1,k2=" << k1 << " " << k2 <<  std::endl;
//
//      double hval = 0;
//      double gval = 0;
//
//      hval = -0.5*(k1+k2);
//      gval = k1*k2;

      // avg k1 and k2 on a process
      //(1a)
//      if (m_s.VN() > 0) k1 = 2*hval/double(m_s.VN()); else k1 = 0;
//      if (m_s.VN() > 0) k2 = gval/double(m_s.VN())/double(m_s.VN()); else k2 = 0;
//
//       (1b)
//          k1 = hval;
//          k2 = gval;
//
//      Js = k1;
//      Ks = k2;
//
//      // std::cout << Dm->rank() << ") Jwn=" << Jwn << std::endl;
//
//      Js = -k1;
//      Ks = -k2;
//
//      awn2_global = area_global;// /=double(npx*npy*npz);
//
//      MPI_Allreduce(&k1,&k1_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
//      MPI_Allreduce(&k2,&k2_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
//
//
//      k1_global/=double(npx*npy*npz);
//      k2_global/=double(npx*npy*npz);
//
//      Js_global = k1_global;
//      Ks_global = k2_global;
//
//      Js_global = -k1_global;
//      Ks_global = -k2_global;

   //   if (0 == Dm->rank()) std::cout << "Jwn_global=" << Jwn_global << " Kwn_global=" << Kwn_global << std::endl;

      vcg::tri::UpdateNormal<sMyMesh>::PerFaceNormalized(m_s);
//
//    double returnEuler = 0;
//
      double EulerChar,nvert,nface;
      double ii, jj, kk;
      int v_count = 0;
      int f_count = 0;
      sEuler = 0.0;
//
//      vcg::tri::UpdateSelection<sMyMesh>::VertexClear(m_s);
//      vcg::tri::UpdateSelection<sMyMesh>::FaceClear(m_s);
//
    size_t wsFacesCount = 0;
    size_t nfaces, nvertices;
    nfaces = nvertices = 0;
//
////     std::cout << " nvert sCompute=" << m_s.VN() << std::endl;
    
    int VIB;
      // make a box
      // select vertices and faces within the box
      // perform Euler characteristic measurement
//      for (int z=1; z<Nz-1; z++) { // kmin kmax
//          for (int y=1; y<Ny-1; y++) {
//              for (int x=1; x<Nx-1; x++) {
//
//                  double i = double(x);//+ 1.0*Dm->iproc()*(Nx-2);
//                  double j = double(y);//+ 1.0*Dm->jproc()*(Ny-2);
//                  double k = double(z);//+ 1.0*Dm->kproc()*(Nz-2);
//
//                  vcg::tri::UpdateSelection<sMyMesh>::FaceClear(m_s);
//                  vcg::tri::UpdateSelection<sMyMesh>::VertexClear(m_s);
//
//                  int n = i + j*Nx + k*Nx*Ny;
//                  v_count = 0;
//                  f_count = 0;
//                  ii = i+1.0;
//                  jj = j+1.0;
//                  kk = k+1.0;
//
//                  vcg::Point3<double> BottomCorner(i,j,k);
//                  vcg::Point3<double> TopCorner(ii,jj,kk);
//                  vcg::Box3<double> RankBox(BottomCorner,TopCorner);
//                  v_count = vcg::tri::UpdateSelection<sMyMesh>::VertexInBox(m_s, RankBox, false);
//                  f_count = vcg::tri::UpdateSelection<sMyMesh>::FaceFromVertexStrict(m_s, false);
//                 // v_count = vcg::tri::UpdateSelection<sMyMesh>::VertexFromFaceStrict(m_s, false);
//
//                  // v_count = vcg::tri::UpdateSelection<sMyMesh>::VertexFromFaceLoose(m_s, false);
//
//
//                  if (f_count > 0) {
//
//                      std::cout << "* n=" << n << "(" << i << "," << j << "," << k << ") nface=" << f_count << " nvert=" << v_count << std::endl;
//                          // Number of faces is number of triangles
//                      nface = double(f_count);
//                      nfaces += f_count;
//                      // Each vertex is shared by four cubes
//                      nvert = double(v_count);
//                      nvertices += v_count;
//                      // Subtract shared sides to avoid double counting
//                      double nside_extern = double(v_count);
//                      double nside_intern = double(v_count)-3.0;
//                      EulerChar=0.0;
//                      if (nvert > 0)    EulerChar = (0.25*nvert - nside_intern - 0.5*nside_extern + nface);
//                      else EulerChar = 0.0;
//                      sEuler += EulerChar;
//                  }
//
//
//
//              }
//          }
//      }

}



void TwoPhase::ComputeAccurateNWInterfaceCurvatures(nwMyMesh & m_nw) {
    
    
    int total_ranks = Dm->nprocz() * Dm->nprocy() * Dm->nprocx();
    
    // Determine where the neighbor ranks are for periodic conditions
    bool Periodic_x, Periodic_X, Periodic_y, Periodic_Y, Periodic_z, Periodic_Z;
    Periodic_x = Periodic_X = Periodic_y = Periodic_Y = Periodic_z = Periodic_Z = false;
    
    if (Dm->rank_x() != Dm->rank()-1) Periodic_x = true;
    if (Dm->rank_X() != Dm->rank()+1) Periodic_X = true;
    
    if (Dm->rank_y() != Dm->rank()-Dm->nprocx()) Periodic_y = true;
    if (Dm->rank_Y() != Dm->rank()+Dm->nprocx()) Periodic_Y = true;
    
    if (Dm->rank_z() != Dm->rank()-Dm->nprocx()*Dm->nprocy()) Periodic_z = true;
    if (Dm->rank_Z() != Dm->rank()+Dm->nprocx()*Dm->nprocy()) Periodic_Z = true;
    
    //    if (Dm->rank() == 7) std::cout << "y=" << Periodic_y << " Y=" << Periodic_Y << std::endl;
    //    if (Dm->rank() == 7) std::cout << "z=" << Periodic_z << " Z=" << Periodic_Z << std::endl;
    
    // Create Eigen matrices:
    Eigen::MatrixXi nwFaces;
    Eigen::MatrixXd nwVertices;
    vcg::tri::MeshToMatrix<nwMyMesh>::GetTriMeshData(m_nw,nwFaces,nwVertices);
    nwFaces.transposeInPlace();
    nwVertices.transposeInPlace();
    
    //
    //    Gradient computation
    //    MatrixXd V;
    //    MatrixXi F;
    
    //    // Load a mesh in OFF format
    //    igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", V, F);
    
    // Read scalar function values from a file, U: #V by 1
    //  VectorXd U; // Populate with values on vertex
    
    //    igl::readDMAT(TUTORIAL_SHARED_PATH "/cheburashka-scalar.dmat",U);
    
    // Compute gradient operator: #F*3 by #V
    //    SparseMatrix<double> G;
    //    igl::grad(V,F,G);
    //
    //    // Compute gradient of U
    //    MatrixXd GU = Map<const MatrixXd>((G*U).eval().data(),F.rows(),3);
    //    // Compute gradient magnitude
    //    const VectorXd GU_mag = GU.rowwise().norm();
    
    int *SendFaceNumber; SendFaceNumber = new int[6]; for (int i=0; i < 6; i++ ) { SendFaceNumber[i] = nwFaces.rows() * nwFaces.cols(); }
    int *SendVertexNumber; SendVertexNumber = new int[6]; for (int i=0; i < 6; i++ ) { SendVertexNumber[i] = nwVertices.rows() * nwVertices.cols(); }
    
    int *RecvFaceNumber; RecvFaceNumber = new int[6]; for (int i=0; i < 6; i++ ) { RecvFaceNumber[i] = -i; }
    int *RecvVertexNumber; RecvVertexNumber = new int[6]; for (int i=0; i < 6; i++ ) { RecvVertexNumber[i] = -i; }
    
    MPI_Request req1[12];
    MPI_Request req2[12];
    
    MPI_Status stat1[12];
    MPI_Status stat2[12];
    
    MPI_Isend(&SendFaceNumber[0],1,MPI_INT, Dm->rank_x(), 12, Dm->Comm,&req1[0]);
    MPI_Isend(&SendFaceNumber[1],1,MPI_INT, Dm->rank_X(), 12, Dm->Comm,&req1[1]);
    MPI_Isend(&SendFaceNumber[2],1,MPI_INT, Dm->rank_y(), 12, Dm->Comm,&req1[2]);
    MPI_Isend(&SendFaceNumber[3],1,MPI_INT, Dm->rank_Y(), 12, Dm->Comm,&req1[3]);
    MPI_Isend(&SendFaceNumber[4],1,MPI_INT, Dm->rank_z(), 12, Dm->Comm,&req1[4]);
    MPI_Isend(&SendFaceNumber[5],1,MPI_INT, Dm->rank_Z(), 12, Dm->Comm,&req1[5]);
    
    MPI_Irecv(&RecvFaceNumber[0],1,MPI_INT, Dm->rank_x(),12,Dm->Comm,&req2[0]);
    MPI_Irecv(&RecvFaceNumber[1],1,MPI_INT, Dm->rank_X(),12,Dm->Comm,&req2[1]);
    MPI_Irecv(&RecvFaceNumber[2],1,MPI_INT, Dm->rank_y(),12,Dm->Comm,&req2[2]);
    MPI_Irecv(&RecvFaceNumber[3],1,MPI_INT, Dm->rank_Y(),12,Dm->Comm,&req2[3]);
    MPI_Irecv(&RecvFaceNumber[4],1,MPI_INT, Dm->rank_z(),12,Dm->Comm,&req2[4]);
    MPI_Irecv(&RecvFaceNumber[5],1,MPI_INT, Dm->rank_Z(),12,Dm->Comm,&req2[5]);
    
    MPI_Isend(&SendVertexNumber[0],1,MPI_INT, Dm->rank_x(), 13, Dm->Comm,&req1[6]);
    MPI_Isend(&SendVertexNumber[1],1,MPI_INT, Dm->rank_X(), 13, Dm->Comm,&req1[7]);
    MPI_Isend(&SendVertexNumber[2],1,MPI_INT, Dm->rank_y(), 13, Dm->Comm,&req1[8]);
    MPI_Isend(&SendVertexNumber[3],1,MPI_INT, Dm->rank_Y(), 13, Dm->Comm,&req1[9]);
    MPI_Isend(&SendVertexNumber[4],1,MPI_INT, Dm->rank_z(), 13, Dm->Comm,&req1[10]);
    MPI_Isend(&SendVertexNumber[5],1,MPI_INT, Dm->rank_Z(), 13, Dm->Comm,&req1[11]);
    
    MPI_Irecv(&RecvVertexNumber[0],1,MPI_INT, Dm->rank_x(),13,Dm->Comm,&req2[6]);
    MPI_Irecv(&RecvVertexNumber[1],1,MPI_INT, Dm->rank_X(),13,Dm->Comm,&req2[7]);
    MPI_Irecv(&RecvVertexNumber[2],1,MPI_INT, Dm->rank_y(),13,Dm->Comm,&req2[8]);
    MPI_Irecv(&RecvVertexNumber[3],1,MPI_INT, Dm->rank_Y(),13,Dm->Comm,&req2[9]);
    MPI_Irecv(&RecvVertexNumber[4],1,MPI_INT, Dm->rank_z(),13,Dm->Comm,&req2[10]);
    MPI_Irecv(&RecvVertexNumber[5],1,MPI_INT, Dm->rank_Z(),13,Dm->Comm,&req2[11]);
    
    MPI_Waitall(12,req1,stat1);
    MPI_Waitall(12,req2,stat2);
    MPI_Barrier(Dm->Comm);
    
    //    if (Dm->rank() == 0) {
    //        for (int i=0; i < 6; i++ ) { std::cout << Dm->rank() << ") recvd nums=" << RecvFaceNumber[i] << std::endl; }
    //     }
    
    int * NFaces_x; if (RecvFaceNumber[0] > 0) NFaces_x = new int[RecvFaceNumber[0]];
    int * NFaces_X; if (RecvFaceNumber[1] > 0) NFaces_X = new int[RecvFaceNumber[1]];
    int * NFaces_y; if (RecvFaceNumber[2] > 0) NFaces_y = new int[RecvFaceNumber[2]];
    int * NFaces_Y; if (RecvFaceNumber[3] > 0) NFaces_Y = new int[RecvFaceNumber[3]];
    int * NFaces_z; if (RecvFaceNumber[4] > 0) NFaces_z = new int[RecvFaceNumber[4]];
    int * NFaces_Z; if (RecvFaceNumber[5] > 0) NFaces_Z = new int[RecvFaceNumber[5]];
    
    double * NVertices_x; if (RecvVertexNumber[0] > 0) NVertices_x = new double[RecvVertexNumber[0]];
    double * NVertices_X; if (RecvVertexNumber[1] > 0) NVertices_X = new double[RecvVertexNumber[1]];
    double * NVertices_y; if (RecvVertexNumber[2] > 0) NVertices_y = new double[RecvVertexNumber[2]];
    double * NVertices_Y; if (RecvVertexNumber[3] > 0) NVertices_Y = new double[RecvVertexNumber[3]];
    double * NVertices_z; if (RecvVertexNumber[4] > 0) NVertices_z = new double[RecvVertexNumber[4]];
    double * NVertices_Z; if (RecvVertexNumber[5] > 0) NVertices_Z = new double[RecvVertexNumber[5]];
    
    MPI_Request req3[12];
    MPI_Request req4[12];
    
    MPI_Status stat3[12];
    MPI_Status stat4[12];
    
    MPI_Isend(nwFaces.data(), SendFaceNumber[0], MPI_INT, Dm->rank_x(), 14, Dm->Comm,&req3[0]);
    MPI_Isend(nwFaces.data(), SendFaceNumber[1], MPI_INT, Dm->rank_X(), 14, Dm->Comm,&req3[1]);
    MPI_Isend(nwFaces.data(), SendFaceNumber[2], MPI_INT, Dm->rank_y(), 14, Dm->Comm,&req3[2]);
    MPI_Isend(nwFaces.data(), SendFaceNumber[3], MPI_INT, Dm->rank_Y(), 14, Dm->Comm,&req3[3]);
    MPI_Isend(nwFaces.data(), SendFaceNumber[4], MPI_INT, Dm->rank_z(), 14, Dm->Comm,&req3[4]);
    MPI_Isend(nwFaces.data(), SendFaceNumber[5], MPI_INT, Dm->rank_Z(), 14, Dm->Comm,&req3[5]);
    
    MPI_Irecv(NFaces_x, RecvFaceNumber[0], MPI_INT, Dm->rank_x(), 14, Dm->Comm,&req4[0]);
    MPI_Irecv(NFaces_X, RecvFaceNumber[1], MPI_INT, Dm->rank_X(), 14, Dm->Comm,&req4[1]);
    MPI_Irecv(NFaces_y, RecvFaceNumber[2], MPI_INT, Dm->rank_y(), 14, Dm->Comm,&req4[2]);
    MPI_Irecv(NFaces_Y, RecvFaceNumber[3], MPI_INT, Dm->rank_Y(), 14, Dm->Comm,&req4[3]);
    MPI_Irecv(NFaces_z, RecvFaceNumber[4], MPI_INT, Dm->rank_z(), 14, Dm->Comm,&req4[4]);
    MPI_Irecv(NFaces_Z, RecvFaceNumber[5], MPI_INT, Dm->rank_Z(), 14, Dm->Comm,&req4[5]);
    
    MPI_Isend(nwVertices.data(), SendVertexNumber[0], MPI_DOUBLE, Dm->rank_x(), 15, Dm->Comm,&req3[6]);
    MPI_Isend(nwVertices.data(), SendVertexNumber[1], MPI_DOUBLE, Dm->rank_X(), 15, Dm->Comm,&req3[7]);
    MPI_Isend(nwVertices.data(), SendVertexNumber[2], MPI_DOUBLE, Dm->rank_y(), 15, Dm->Comm,&req3[8]);
    MPI_Isend(nwVertices.data(), SendVertexNumber[3], MPI_DOUBLE, Dm->rank_Y(), 15, Dm->Comm,&req3[9]);
    MPI_Isend(nwVertices.data(), SendVertexNumber[4], MPI_DOUBLE, Dm->rank_z(), 15, Dm->Comm,&req3[10]);
    MPI_Isend(nwVertices.data(), SendVertexNumber[5], MPI_DOUBLE, Dm->rank_Z(), 15, Dm->Comm,&req3[11]);
    
    MPI_Irecv(NVertices_x, RecvVertexNumber[0], MPI_DOUBLE, Dm->rank_x(), 15, Dm->Comm,&req4[6]);
    MPI_Irecv(NVertices_X, RecvVertexNumber[1], MPI_DOUBLE, Dm->rank_X(), 15, Dm->Comm,&req4[7]);
    MPI_Irecv(NVertices_y, RecvVertexNumber[2], MPI_DOUBLE, Dm->rank_y(), 15, Dm->Comm,&req4[8]);
    MPI_Irecv(NVertices_Y, RecvVertexNumber[3], MPI_DOUBLE, Dm->rank_Y(), 15, Dm->Comm,&req4[9]);
    MPI_Irecv(NVertices_z, RecvVertexNumber[4], MPI_DOUBLE, Dm->rank_z(), 15, Dm->Comm,&req4[10]);
    MPI_Irecv(NVertices_Z, RecvVertexNumber[5], MPI_DOUBLE, Dm->rank_Z(), 15, Dm->Comm,&req4[11]);
    
    MPI_Waitall(12,req3,stat3);
    MPI_Waitall(12,req4,stat4);
    MPI_Barrier(Dm->Comm);
    
    
    //        if (Dm->rank() == 0) {
    //            for (int i=0; i < 6; i++ ) { std::cout << Dm->rank() << ") recvd nums=" << RecvVertexNumber[i] << std::endl; }
    //         }
    
    //    // Create a new mesh with the Eigen data:
    Point RootPtA,RootPtB,RootPtC;
    nwRootMesh m_root_nw;
    nwRootMesh::VertexIterator root_vi;
    nwRootMesh::FaceIterator root_fi;
    
    
    for (int aa = 0; aa < RecvVertexNumber[0]; aa+=9 ) {
        RootPtA.x = NVertices_x[aa];
        RootPtA.y = NVertices_x[aa+1];
        RootPtA.z = NVertices_x[aa+2];
        
        RootPtB.x = NVertices_x[aa+3];
        RootPtB.y = NVertices_x[aa+4];
        RootPtB.z = NVertices_x[aa+5];
        
        RootPtC.x = NVertices_x[aa+6];
        RootPtC.y = NVertices_x[aa+7];
        RootPtC.z = NVertices_x[aa+8];
        
        if (Periodic_x) {
            RootPtA.x -= 1.0*Dm->nprocx()*(Nx-2);
            RootPtB.x -= 1.0*Dm->nprocx()*(Nx-2);
            RootPtC.x -= 1.0*Dm->nprocx()*(Nx-2);
        }
        
        nwRootMesh::VertexPointer root_ivp[3];
        
        root_vi = vcg::tri::Allocator<nwRootMesh>::AddVertices(m_root_nw,3);
        root_fi = vcg::tri::Allocator<nwRootMesh>::AddFaces(m_root_nw,1);
        root_ivp[0]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtA.x, RootPtA.y, RootPtA.z); ++root_vi;
        root_ivp[1]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtB.x, RootPtB.y, RootPtB.z); ++root_vi;
        root_ivp[2]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtC.x, RootPtC.y, RootPtC.z); ++root_vi;
        
        root_fi->V(0)=root_ivp[0];
        root_fi->V(1)=root_ivp[1];
        root_fi->V(2)=root_ivp[2];
    }
    
    for (int aa = 0; aa < RecvVertexNumber[1]; aa+=9 ) {
        RootPtA.x = NVertices_X[aa];
        RootPtA.y = NVertices_X[aa+1];
        RootPtA.z = NVertices_X[aa+2];
        
        RootPtB.x = NVertices_X[aa+3];
        RootPtB.y = NVertices_X[aa+4];
        RootPtB.z = NVertices_X[aa+5];
        
        RootPtC.x = NVertices_X[aa+6];
        RootPtC.y = NVertices_X[aa+7];
        RootPtC.z = NVertices_X[aa+8];
        
        if (Periodic_X) {
            RootPtA.x += 1.0*Dm->nprocx()*(Nx-2);
            RootPtB.x += 1.0*Dm->nprocx()*(Nx-2);
            RootPtC.x += 1.0*Dm->nprocx()*(Nx-2);
        }
        
        
        nwRootMesh::VertexPointer root_ivp[3];
        
        root_vi = vcg::tri::Allocator<nwRootMesh>::AddVertices(m_root_nw,3);
        root_fi = vcg::tri::Allocator<nwRootMesh>::AddFaces(m_root_nw,1);
        root_ivp[0]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtA.x, RootPtA.y, RootPtA.z); ++root_vi;
        root_ivp[1]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtB.x, RootPtB.y, RootPtB.z); ++root_vi;
        root_ivp[2]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtC.x, RootPtC.y, RootPtC.z); ++root_vi;
        
        root_fi->V(0)=root_ivp[0];
        root_fi->V(1)=root_ivp[1];
        root_fi->V(2)=root_ivp[2];
    }
    
    for (int aa = 0; aa < RecvVertexNumber[2]; aa+=9 ) {
        RootPtA.x = NVertices_y[aa];
        RootPtA.y = NVertices_y[aa+1];
        RootPtA.z = NVertices_y[aa+2];
        
        RootPtB.x = NVertices_y[aa+3];
        RootPtB.y = NVertices_y[aa+4];
        RootPtB.z = NVertices_y[aa+5];
        
        RootPtC.x = NVertices_y[aa+6];
        RootPtC.y = NVertices_y[aa+7];
        RootPtC.z = NVertices_y[aa+8];
        
        if (Periodic_y) {
            RootPtA.y -= 1.0*Dm->nprocy()*(Ny-2);
            RootPtB.y -= 1.0*Dm->nprocy()*(Ny-2);
            RootPtC.y -= 1.0*Dm->nprocy()*(Ny-2);
        }
        
        nwRootMesh::VertexPointer root_ivp[3];
        
        root_vi = vcg::tri::Allocator<nwRootMesh>::AddVertices(m_root_nw,3);
        root_fi = vcg::tri::Allocator<nwRootMesh>::AddFaces(m_root_nw,1);
        root_ivp[0]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtA.x, RootPtA.y, RootPtA.z); ++root_vi;
        root_ivp[1]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtB.x, RootPtB.y, RootPtB.z); ++root_vi;
        root_ivp[2]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtC.x, RootPtC.y, RootPtC.z); ++root_vi;
        
        root_fi->V(0)=root_ivp[0];
        root_fi->V(1)=root_ivp[1];
        root_fi->V(2)=root_ivp[2];
    }
    
    for (int aa = 0; aa < RecvVertexNumber[3]; aa+=9 ) {
        RootPtA.x = NVertices_Y[aa];
        RootPtA.y = NVertices_Y[aa+1];
        RootPtA.z = NVertices_Y[aa+2];
        
        RootPtB.x = NVertices_Y[aa+3];
        RootPtB.y = NVertices_Y[aa+4];
        RootPtB.z = NVertices_Y[aa+5];
        
        RootPtC.x = NVertices_Y[aa+6];
        RootPtC.y = NVertices_Y[aa+7];
        RootPtC.z = NVertices_Y[aa+8];
        
        if (Periodic_Y) {
            RootPtA.y += 1.0*Dm->nprocy()*(Ny-2);
            RootPtB.y += 1.0*Dm->nprocy()*(Ny-2);
            RootPtC.y += 1.0*Dm->nprocy()*(Ny-2);
        }
        
        nwRootMesh::VertexPointer root_ivp[3];
        
        root_vi = vcg::tri::Allocator<nwRootMesh>::AddVertices(m_root_nw,3);
        root_fi = vcg::tri::Allocator<nwRootMesh>::AddFaces(m_root_nw,1);
        root_ivp[0]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtA.x, RootPtA.y, RootPtA.z); ++root_vi;
        root_ivp[1]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtB.x, RootPtB.y, RootPtB.z); ++root_vi;
        root_ivp[2]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtC.x, RootPtC.y, RootPtC.z); ++root_vi;
        
        root_fi->V(0)=root_ivp[0];
        root_fi->V(1)=root_ivp[1];
        root_fi->V(2)=root_ivp[2];
    }
    
    for (int aa = 0; aa < RecvVertexNumber[4]; aa+=9 ) {
        RootPtA.x = NVertices_z[aa];
        RootPtA.y = NVertices_z[aa+1];
        RootPtA.z = NVertices_z[aa+2];
        
        RootPtB.x = NVertices_z[aa+3];
        RootPtB.y = NVertices_z[aa+4];
        RootPtB.z = NVertices_z[aa+5];
        
        RootPtC.x = NVertices_z[aa+6];
        RootPtC.y = NVertices_z[aa+7];
        RootPtC.z = NVertices_z[aa+8];
        
        if (Periodic_z) {
            RootPtA.z -= 1.0*Dm->nprocz()*(Nz-2);
            RootPtB.z -= 1.0*Dm->nprocz()*(Nz-2);
            RootPtC.z -= 1.0*Dm->nprocz()*(Nz-2);
        }
        
        nwRootMesh::VertexPointer root_ivp[3];
        
        root_vi = vcg::tri::Allocator<nwRootMesh>::AddVertices(m_root_nw,3);
        root_fi = vcg::tri::Allocator<nwRootMesh>::AddFaces(m_root_nw,1);
        root_ivp[0]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtA.x, RootPtA.y, RootPtA.z); ++root_vi;
        root_ivp[1]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtB.x, RootPtB.y, RootPtB.z); ++root_vi;
        root_ivp[2]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtC.x, RootPtC.y, RootPtC.z); ++root_vi;
        
        root_fi->V(0)=root_ivp[0];
        root_fi->V(1)=root_ivp[1];
        root_fi->V(2)=root_ivp[2];
    }
    
    for (int aa = 0; aa < RecvVertexNumber[5]; aa+=9 ) {
        RootPtA.x = NVertices_Z[aa];
        RootPtA.y = NVertices_Z[aa+1];
        RootPtA.z = NVertices_Z[aa+2];
        
        RootPtB.x = NVertices_Z[aa+3];
        RootPtB.y = NVertices_Z[aa+4];
        RootPtB.z = NVertices_Z[aa+5];
        
        RootPtC.x = NVertices_Z[aa+6];
        RootPtC.y = NVertices_Z[aa+7];
        RootPtC.z = NVertices_Z[aa+8];
        
        if (Periodic_Z) {
            RootPtA.z += 1.0*Dm->nprocz()*(Nz-2);
            RootPtB.z += 1.0*Dm->nprocz()*(Nz-2);
            RootPtC.z += 1.0*Dm->nprocz()*(Nz-2);
        }
        
        nwRootMesh::VertexPointer root_ivp[3];
        
        root_vi = vcg::tri::Allocator<nwRootMesh>::AddVertices(m_root_nw,3);
        root_fi = vcg::tri::Allocator<nwRootMesh>::AddFaces(m_root_nw,1);
        root_ivp[0]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtA.x, RootPtA.y, RootPtA.z); ++root_vi;
        root_ivp[1]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtB.x, RootPtB.y, RootPtB.z); ++root_vi;
        root_ivp[2]=&*root_vi;  root_vi->P()=nwRootMesh::CoordType ( RootPtC.x, RootPtC.y, RootPtC.z); ++root_vi;
        
        root_fi->V(0)=root_ivp[0];
        root_fi->V(1)=root_ivp[1];
        root_fi->V(2)=root_ivp[2];
    }
    
    vcg::tri::Append<nwRootMesh, nwMyMesh>::Mesh(m_root_nw, m_nw, false  , false );
    
    
    
    // SECTION: Clean mesh and compute windows
    // Mesh clean
    vcg::tri::Clean<nwRootMesh>::RemoveDuplicateFace(m_root_nw);
    vcg::tri::Clean<nwRootMesh>::RemoveDuplicateVertex(m_root_nw,false);
    

    int degenface = vcg::tri::Clean<nwRootMesh>::RemoveDegenerateFace(m_root_nw);
    vcg::tri::UpdateTopology<nwRootMesh>::FaceFace(m_root_nw);
    int degenvertex = vcg::tri::Clean<nwRootMesh>::RemoveDegenerateVertex(m_root_nw);
    vcg::tri::UpdateTopology<nwRootMesh>::VertexFace(m_root_nw);
    
    // Need to remove any non 2-manifold faces
    int removeNonManifoldFace = vcg::tri::Clean<nwRootMesh>::RemoveNonManifoldFace(m_root_nw);
    
    std::pair<int,int> ws_CCs;
    ws_CCs =  vcg::tri::Clean<nwRootMesh>::RemoveSmallConnectedComponentsDiameter(m_root_nw, 3.5);
    ws_CCs.first = ws_CCs.first - ws_CCs.second;
    
    
    int removeUnreferencedVertex = vcg::tri::Clean<nwRootMesh>::RemoveUnreferencedVertex(m_root_nw, true);
    
    removeNonManifoldFace = vcg::tri::Clean<nwRootMesh>::RemoveNonManifoldFace(m_root_nw);
    
    degenvertex = vcg::tri::Clean<nwRootMesh>::RemoveDegenerateVertex(m_root_nw);
    
    vcg::tri::Clean<nwRootMesh>::RemoveDuplicateVertex(m_root_nw,false);
    
    // Last step before compact: Can be determined using the meshlab help: Compute principal curvature directiows
    removeUnreferencedVertex = vcg::tri::Clean<nwRootMesh>::RemoveUnreferencedVertex(m_root_nw, true);
    
    vcg::tri::UpdateFlags<nwRootMesh>::VertexSetB(m_root_nw);
    vcg::tri::UpdateSelection<nwRootMesh>::VertexFromBorderFlag(m_root_nw,false  );
    
    vcg::tri::UpdateTopology<nwRootMesh>::AllocateEdge(m_root_nw);
    
    // Final step: Faces, vertices, and edges must be compact for curvature calculation
    vcg::tri::Allocator<nwRootMesh>::CompactEveryVector(m_root_nw);
    
    vcg::tri::UpdateCurvatureFitting<nwRootMesh>::computeCurvature(m_root_nw);
    
    vcg::tri::UpdateSelection<nwRootMesh>::VertexClear(m_root_nw);
    vcg::tri::UpdateSelection<nwRootMesh>::FaceClear(m_root_nw);
    
    double rank_mean_curvature = 0;
    double rank_Gaussian_curvature = 0;
    double rank_mean_curvature_global = 0;
    double rank_Gaussian_curvature_global = 0;
    
    // Save curvatures in new registers
    
    int x_rank_min;
    int y_rank_min;
    int z_rank_min;
    int x_rank_max;
    int y_rank_max;
    int z_rank_max;
    int count = 0;
    
    x_rank_min = Dm->iproc() *(Nx-2) + 1;
    y_rank_min = Dm->jproc() *(Ny-2) + 1;
    z_rank_min = Dm->kproc() *(Nz-2) + 1;
    x_rank_max = (Dm->iproc()+1) *(Nx-2) + 1;
    y_rank_max = (Dm->jproc()+1) *(Ny-2) + 1;
    z_rank_max = (Dm->kproc()+1) *(Nz-2) + 1;
    vcg::Point3<double> BottomCorner(x_rank_min,y_rank_min,z_rank_min);
    vcg::Point3<double> TopCorner(x_rank_max,y_rank_max,z_rank_max);
    vcg::Box3<double> RankBox(BottomCorner,TopCorner);
    
    
    vcg::tri::UpdateSelection<nwRootMesh>::VertexInBox(m_root_nw, RankBox, true);
    vcg::tri::UpdateSelection<nwRootMesh>::FaceFromVertexStrict(m_root_nw, true);
    
    
    
    // mean and Gaussian curvature calculatiows
    double k1 = 0.0;
    double k2 = 0.0;
    double k1_global = 0;
    double k2_global = 0;
    
    int k1nan = 0;
    int k2nan = 0;
    double tmp = 0;
    size_t v_count = 0;
    for(auto vi = m_root_nw.vert.begin(); vi != m_root_nw.vert.end(); ++vi) if((*vi).IsS()) {
        tmp = (*vi).K1();
        if (isnan(tmp)) {k1nan+=1;}
        else k1 += tmp;
        tmp = (*vi).K2();
        if (isnan(tmp)) {k2nan+=1;}
        else k2 += tmp;
        v_count++;
    }
    //   std::cout << "(" << Dm->rank() << ") WS: k1nan=" << k1nan << " k2nan=" << k2nan << " k1,k2=" << k1 << " " << k2 <<  std::endl;
    
    double hval = 0;
    double gval = 0;
    
    hval = -0.5*(k1+k2);
    gval = k1*k2;
    
    // avg k1 and k2 on a process
    (v_count > 0) ? k1 = 2.0*hval/double(v_count) : 0;
    (v_count > 0) ? k2 = gval/double(v_count)/double(v_count) : 0;
    
    rank_mean_curvature= k1;
    rank_Gaussian_curvature = k2;
    
    Jwn2 = rank_mean_curvature;
    Kwn2 = rank_Gaussian_curvature;
    
    
    MPI_Allreduce(&rank_mean_curvature,&rank_mean_curvature_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&rank_Gaussian_curvature,&rank_Gaussian_curvature_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    rank_mean_curvature_global/=double(total_ranks);
    rank_Gaussian_curvature_global/=double(total_ranks);
    
    Jwn2_global = rank_mean_curvature_global;
    Kwn2_global = rank_Gaussian_curvature_global;
    

    delete[] SendFaceNumber;
    delete[] SendVertexNumber;
      delete[] RecvFaceNumber;
     delete[] RecvVertexNumber;
 
    delete[] NFaces_x;
    delete[] NFaces_X;
    delete[] NFaces_y;
    delete[] NFaces_Y;
    delete[] NFaces_z;
    delete[] NFaces_Z;
      
    delete[] NVertices_x;
    delete[] NVertices_X;
    delete[] NVertices_y;
    delete[] NVertices_Y;
    delete[] NVertices_z;
    delete[] NVertices_Z;
    
    
    
}




void TwoPhase::wComputeQuantities(wMyMesh & m_w) {
    
    /* First steps are to delete any duplicate faces and vertices */
    int beforeVN, beforeFN, wAfterVN, wAfterFN; beforeVN = wAfterVN = beforeFN = wAfterFN = 0;  beforeVN = m_w.VN(); beforeFN = m_w.FN();
    vcg::tri::Clean<wMyMesh>::RemoveDuplicateFace(m_w);
    vcg::tri::Clean<wMyMesh>::RemoveDuplicateVertex(m_w,false);
    
    /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
    int degenface = vcg::tri::Clean<wMyMesh>::RemoveDegenerateFace(m_w);
    vcg::tri::UpdateTopology<wMyMesh>::FaceFace(m_w);
    int degenvertex = vcg::tri::Clean<wMyMesh>::RemoveDegenerateVertex(m_w);
    vcg::tri::UpdateTopology<wMyMesh>::VertexFace(m_w);
    
    /* Need to remove any non 2-manifold faces */
    int removeNonManifoldFace = vcg::tri::Clean<wMyMesh>::RemoveNonManifoldFace(m_w);
    
    std::pair<int,int> w_CCs;
    w_CCs =  vcg::tri::Clean<wMyMesh>::RemoveSmallConnectedComponentsSize(m_w, 0);
    w_CCs.first = w_CCs.first - w_CCs.second;
    
    
    
    /* Last step before compact: Can be determined using the meshlab help: Compute principal curvature directiow*/
    int removeUnreferencedVertex = vcg::tri::Clean<wMyMesh>::RemoveUnreferencedVertex(m_w, true);
    
    int wGlobal_CCs = 0;
    vcg::tri::UpdateFlags<wMyMesh>::VertexSetB(m_w);
    int wBorderVert = vcg::tri::UpdateSelection<wMyMesh>::VertexFromBorderFlag(m_w,false /* keep the selection */);
    int w_dS = 0;
    int wGlobal_dS = 0;
    if (wBorderVert > 0) w_dS = 1;
    //  printf("%i wBorderVert=%d\n",Dm->rank(),wBorderVert);
    MPI_Allreduce(&w_dS,&wGlobal_dS,1,MPI_INT,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&w_CCs.first,&wGlobal_CCs,1,MPI_INT,MPI_SUM,Dm->Comm);
    
    /* Get all edges */
    vcg::tri::UpdateTopology<wMyMesh>::AllocateEdge(m_w);
    
    /* Final step: Faces, vertices, and edges must be compact for curvature calculation */
    vcg::tri::Allocator<wMyMesh>::CompactEveryVector(m_w);
    
    bool IsOrientable = false;
    bool IsOriented = false;
    //vcg::tri::Clean<wMyMesh>::OrientCoherentlyMesh(m_w, IsOriented, IsOrientable);
    
    wAfterVN = m_w.VN(); wAfterFN = m_w.FN();
    /* Get Edge count for Euler characteristic */
    size_t wEN = m_w.EN();
    
    
    
    int z,y,x;
    
    if (wAfterFN > 0) {
        /* Allocate new mesh and vertex iterator */
        wMyMesh::VertexIterator vi;
        /* Create convex hull mesh object */
        chMyMesh CH_w;
        
        Point Q;
        
        const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
        
        // If external boundary conditions are set, do not average over the inlet
        int kmin=1; int kmax=Nz-1;
        
        /* Add corner points to surface for volume calculations */
        for (int k=kmin; k<kmax; k++){ // kmin kmax
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++) {
                    for (int p = 0; p < 8; p++) {
                        int m = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                        if ( SDn(m) < 0 ) {
                            z = m/(Nx*Ny); y = (m-Nx*Ny*z)/Nx; x = m-Nx*Ny*z-Nx*y;
                            Q.z = z; Q.z += 1.0*Dm->kproc()*(Nx-2);
                            Q.y = y; Q.y += 1.0*Dm->jproc()*(Nx-2);
                            Q.x = x; Q.x += 1.0*Dm->iproc()*(Nx-2);
                            
                            vi = vcg::tri::Allocator<wMyMesh>::AddVertices(m_w,1);
                            vi->P()=wMyMesh::CoordType ( Q.x, Q.y, Q.z); vi++;
                        }
                    }
                }
            }
        }
        
        beforeVN = wAfterVN = beforeFN = wAfterFN = 0;
        beforeVN = m_w.VN(); beforeFN = m_w.FN();
        vcg::tri::Clean<wMyMesh>::RemoveDuplicateFace(m_w);
        vcg::tri::Clean<wMyMesh>::RemoveDuplicateVertex(m_w,false);
        
        /* Generate convex hull mesh */
        vcg::tri::ConvexHull<wMyMesh,chMyMesh>::ComputeConvexHull(m_w,CH_w);
        
        /* For visualization purposes */
        vcg::tri::Append<wMyMesh, chMyMesh>::Mesh(m_w, CH_w, false /* selected */, false /* adjFlag */);
        
        /* Compute volume fraction of phase */
        vcg::tri::Inertia<chMyMesh> CH_inertia(CH_w);
        double CH_vol = CH_inertia.Mass();
        
        vcg::tri::UpdateNormal<chMyMesh>::PerFaceNormalized(CH_w);
        
        //   std::cout << "w volume=" << CH_vol << std::endl;
    }
}


void TwoPhase::nComputeQuantities(nMyMesh & m_n) {
    
   int beforeVN, beforeFN, nAfterVN, nAfterFN;
     beforeVN = nAfterVN = beforeFN = nAfterFN = 0;
     beforeVN = m_n.VN(); beforeFN = m_n.FN();


    vcg::tri::Clean<nMyMesh>::RemoveDuplicateFace(m_n);
     vcg::tri::Clean<nMyMesh>::RemoveDuplicateVertex(m_n,false);

     size_t count = 0;

     /* Before Updating FaceFace and VertexFace components, remove any degenerate faces or vertices */
     int degenface = vcg::tri::Clean<nMyMesh>::RemoveDegenerateFace(m_n);
     vcg::tri::UpdateTopology<nMyMesh>::FaceFace(m_n);
     int degenvertex = vcg::tri::Clean<nMyMesh>::RemoveDegenerateVertex(m_n);
     vcg::tri::UpdateTopology<nMyMesh>::VertexFace(m_n);

     /* FINISH REPAIRING STEPS */
     // int splitNonManifoldVertex = vcg::tri::Clean<nMyMesh>::SplitNonManifoldVertex(m_n, 5);

     /* Need to remove any non 2-manifold faces */
     int removeNonManifoldFace = vcg::tri::Clean<nMyMesh>::RemoveNonManifoldFace(m_n);

     std::pair<int,int> n_CCs;
     n_CCs =  vcg::tri::Clean<nMyMesh>::RemoveSmallConnectedComponentsDiameter(m_n, 3.0 /* seems to be a good threshold based upon visual experiments in MeshLab */);
     n_CCs.first = n_CCs.first - n_CCs.second;

     //  std::cout << Dm->rank() << ") CC=" << n_CCs.first << std::endl;

     int removeUnreferencedVertex = vcg::tri::Clean<nMyMesh>::RemoveUnreferencedVertex(m_n, true);

     removeNonManifoldFace = vcg::tri::Clean<nMyMesh>::RemoveNonManifoldFace(m_n);

     degenvertex = vcg::tri::Clean<nMyMesh>::RemoveDegenerateVertex(m_n);

     // Final check for Euler
     vcg::tri::Clean<nMyMesh>::RemoveDuplicateFace(m_n);
     vcg::tri::Clean<nMyMesh>::RemoveDuplicateVertex(m_n,false);

     vcg::tri::Allocator<nMyMesh>::CompactEveryVector(m_n);

    

    // clear an:
    an = 0;
    if (nAfterFN > 0) {
        double CH_vol = 0;
        
        double x,y,z;
//        const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
        
        // If external boundary conditions are set, do not average over the inlet
       
        /* Allocate new mesh and vertex iterator */
        nMyMesh::VertexIterator vi;
        /* Create convex hull mesh object */
        chMyMesh CH_n;
        
        Point Q;
        
        const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
        
        // If external boundary conditions are set, do not average over the inlet
    //    int kmin=1; int kmax=Nz-1;
        int kmin=1; int kmax=Nz-1;
        if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=3;
        if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-3;
        
        /* Add corner points to surface for volume calculations */
        for (int k=kmin; k<kmax; k++){ // kmin kmax
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++) {
                    for (int p = 0; p < 8; p++) {
                        int m = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                        if ( SDw(m) > 0 && SDs(m) > 0 ) {
                            z = m/(Nx*Ny); y = (m-Nx*Ny*z)/Nx; x = m-Nx*Ny*z-Nx*y;
                            Q.z = z; Q.z += 1.0*Dm->kproc()*(Nx-2);
                            Q.y = y; Q.y += 1.0*Dm->jproc()*(Nx-2);
                            Q.x = x; Q.x += 1.0*Dm->iproc()*(Nx-2);
                            
                            vi = vcg::tri::Allocator<nMyMesh>::AddVertices(m_n,1);
                            vi->P()=nMyMesh::CoordType ( Q.x, Q.y, Q.z); vi++;
                        }
                    }
                }
            }
        }
        
        beforeVN = nAfterVN = beforeFN = nAfterFN = 0;
        beforeVN = m_n.VN(); beforeFN = m_n.FN();
        vcg::tri::Clean<nMyMesh>::RemoveDuplicateFace(m_n);
        vcg::tri::Clean<nMyMesh>::RemoveDuplicateVertex(m_n,false);
        
        /* Generate convex hull mesh */
        vcg::tri::ConvexHull<nMyMesh,chMyMesh>::ComputeConvexHull(m_n,CH_n);
        
        /* For visualization purposes */
        vcg::tri::Append<nMyMesh, chMyMesh>::Mesh(m_n, CH_n, false /* selected */, false /* adjFlag */);
        
        /* Compute volume fraction of phase */
        vcg::tri::Inertia<chMyMesh> CH_inertia(CH_n);
        CH_vol = CH_inertia.Mass();
        an += CH_vol;
    }
    
    
   
       
}



void TwoPhase::CreateIntegrationBoundaryRegionOne() {
    // Manually enter the integration field data
    
    IntegrationRegionOne.resize(Nx,Ny,Nz); IntegrationRegionOne.fill(0);
}

void TwoPhase::CreateIntegrationBoundaryRegionTwo(DoubleArray & Field) {
   
   
}


