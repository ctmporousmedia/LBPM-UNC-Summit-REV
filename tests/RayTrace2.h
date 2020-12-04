#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <array>
#include <math.h>
#include "common/Array.h"

#include "analysis/TwoPhase.h"

using namespace std;


int poreCount(char * id,
              int Nx, int Ny, int Nz) {
    size_t count = 0;
    int dat;
    size_t m;
    for (int k=1;k<Nz-1;k++){ // printf("k=%d\n",k);
        for (int j=1;j<Ny-1;j++){
            for (int i=1;i<Nx-1;i++){
                m = i + j*Nx + k*Nx*Ny;
                dat = id[m];
                //if (dat > 0) count++;
                count++;
            }
        }
    }
    return count;
}



static inline void fgetl( char * str,
                         int num, FILE * stream ) {
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}



void SetIdAndCount(char * id,
                   size_t & count, double & hx, int & Nx, int & Ny, int & Nz, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz,int rank) {} // Not used

class RayTrace {
public:
    RayTrace() {
        ray_origin = new double[3];
        ray_direction = new double[3];
        geometry = 0;
        
        /* s: geometry = 2 */
        plane_s = new double[3];
        /* */
        
        
        
        
        /* sp: geometry = 1 */
        sphere_center = new double[3];
        sphere_radius = 0;
        o_minus_c = new double[3];
    }
    
    /* Setters and Getters */
    void setSphereRadius(double & rad) {
        sphere_radius = rad;
    }
    void setRayOrigin(double *  tmp) {
        ray_origin[0] = tmp[0];
        ray_origin[1] = tmp[1];
        ray_origin[2] = tmp[2];
    }
    void setSphereCenter(double *  tmp) {
        sphere_center[0] = tmp[0];
        sphere_center[1] = tmp[1];
        sphere_center[2] = tmp[2];
    }
    double * getRayOrigin() {
        return ray_origin;
    }
    void setRayDirection(double *  tmp) {
        ray_direction[0] = tmp[0];
        ray_direction[1] = tmp[1];
        ray_direction[2] = tmp[2];
    }
    
    
    void setPlane_s(double A, double B, double C, double D) {
        plane_s[0] = A;
        plane_s[1] = B;
        plane_s[2] = C;
        D_s = D;
    }
    
    
    double * getRayDirection_s() { return ray_direction; }
    double getRayTraceDistance_s() {
        return ray_trace_s(ray_origin, plane_s, D_s);
    }
    
    double * getRayDirection_id() { return ray_direction; }
    double getRayTraceDistance_id() {
        return ray_trace_id(ray_origin, plane_s, D_s);
    }
    
    double getRayTraceDistance_s_q() {
        return ray_trace_s_q(ray_origin, ray_direction, plane_s, D_s);
    }
    
    
    double * getRayDirection_sp() { return ray_direction; }
    double getRayTraceDistance_sp() {
        return ray_trace_sp(ray_origin, sphere_center, ray_direction, sphere_radius);
    }
    
private:
    // 1 periodic sphere pack, 2 slab, 3 parallel plates, 4 square tube, 5 capillary tube
    int geometry;
    
    /* Operators needed to compute the ray_trace function */
    double dot(double *  vec1, double *  vec2);
    double mag(double *  vec1 /* Get || . ||  of this vector */);
    double magSq(double *  vec1 /* Get || . ||^2 of this vector */);
    void subtract_vectors(double *  vec1 /* ray_origin */, double *  vec2 /* cell_center */,double *  vec3 /* Return vector */);
    
    double ray_trace_id(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* plane translation */);
    double ray_trace_s(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* plane translation */);
    double ray_trace_s_q(double *  ray_origin /* ray_origin */, double * ray_direction,  double *  plane_s /* plane_s */, double D_s /* plane translation */);
    
    double dist_sp(double & term1 /* u_dot_o_minus_c */, double & term2 /* under_sqrt */);
    
    double ray_trace_sp(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */);
    
    
    
private:
    
    /* general ray information */
    double * ray_origin;
    double * ray_direction;
    
    /* sphere pack */
    double * sphere_center;
    double sphere_radius;
    
    /* slab */
    double * plane_s;
    double D_s;
    
    
    double * o_minus_c;
    double u_dot_o_minus_c;
    double dp,dp1,dp2,dp3,dp4;
    double under_sqrt;
    
    double A,B,C,D;
};

double RayTrace::dot(double *  vec1, double *  vec2) {
    return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

double RayTrace::mag(double *  vec1 /* Get || . ||  of this vector */) {
    return sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
}

double RayTrace::magSq(double *  vec1 /* Get || . ||^2 of this vector */) {
    return (vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
}

void RayTrace::subtract_vectors(double *  vec1 /* ray_origin */, double *  vec2 /* cell_center */,
                                double *  vec3 /* Return vector */) {
    vec3[0] = vec1[0] - vec2[0];
    vec3[1] = vec1[1] - vec2[1];
    vec3[2] = vec1[2] - vec2[2];
}



double RayTrace::ray_trace_s(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* displacement */) {
    
    double tmp = mag(vec2); // Denominator of distance
    if (tmp == 0) tmp = 1.0;
    
    
    dp = dot(vec2, vec1);
    
    if (dp >0) return (-(dp - D_s)/tmp);
    else return (-(dp + D_s)/tmp);
    
}

double RayTrace::ray_trace_id(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* displacement */) {
    
    double tmp = mag(vec2); // Denominator of distance
    if (tmp == 0) tmp = 1.0;
    
    
    dp = dot(vec2, vec1);
    
    if (dp >0) return (-(dp - D_s)/tmp);
    else return (-(dp + D_s)/tmp);
    
}

double RayTrace::ray_trace_s_q(double *  lnaught /* ray_origin */, double * l /* ray direction */, double *  plane_s /* plane_s */, double D_s /* displacement */)
{
    //    double tmp = mag(vec2); if (tmp == 0) tmp = 1.0;
    //    subtract_vectors(vec2 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    //    o_minus_c[0] += 0; o_minus_c[1] += 0; o_minus_c[2] += 0;
    //
    //    double * normal; normal = new double[3];
    //    normal[0] = vec2[0]/=tmp;
    //    normal[1] = vec2[1]/=tmp;
    //    normal[2] = vec2[2]/=tmp;
    //
    //    tmp = dot(ray_direction,normal); if (tmp == 0) tmp = 1.0;
    //    dp = dot(o_minus_c, normal) / tmp;
    //
    //    delete[] normal;
    //
    //    // if (dp < 1.5) cout << "dp=" << dp << endl;
    //
    //    if (dp >0) return (-(dp - D_s)/tmp);
    //    else return (-(dp + D_s)/tmp);
    
    /* ptrs */
    double * pnaught; pnaught = new double[3];
    double * normal; normal = new double[3];
    
    /* registers */
    double p0x,p0y,p0z;
    double tmp;
    double denominator = 1;
    double numerator = 1;
    
    /* get normalized ray direction */
    tmp = mag(l); if (tmp == 0) tmp = 1.0;
    l[0] /= tmp; l[1] /= tmp;  l[2] /= tmp;
    
    tmp = mag(plane_s); if (tmp == 0) tmp = 1.0; // used if plane is rotated
    /* normal vector to the plane pp1 */
    normal[0] = plane_s[0]/tmp;
    normal[1] = plane_s[1]/tmp;
    normal[2] = plane_s[2]/tmp;
    
    
    /* distance point from world origin (needs to be parallelized!) */
    p0x = 0;  p0y = 0;  p0z = D_s;
    pnaught[0] = p0x;  pnaught[1] = p0y;  pnaught[2] = p0z;
    
    /* l0 - p0 */
    subtract_vectors(pnaught /* ray origin "lo" */,
                     lnaught /* distance from world origin point "" */,
                     o_minus_c /* Return: "l0-p0"*/);
    
    /* -(p0-l0).n */
    numerator = dot(o_minus_c,normal); if (tmp == 0) tmp = 1.0;
    /* -(p0-l0).n / (l.n) */
    denominator = dot(l,normal);
    if (denominator == 0.0) denominator = 1.0;
    
    dp1 = numerator / denominator + 0.5;  // Result for plane 1
    
    
    
    return dp1;
}


void ComputeSlabQDistances(double * qDistances, int & Nx, int & Ny, int & Nz, double & Lx,int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, double A, double B, double C, double D, double offset) {
    double D3Q19[18][3] = {{1,0,0},{-1,0,0},
        {0,1,0},{0,-1,0},
        {0,0,1},{0,0,-1},
        {1,1,0},{-1,-1,0},
        {1,-1,0},{-1,1,0},
        {1,0,1},{-1,0,-1},
        {1,0,-1},{-1,0,1},
        {0,1,1},{0,-1,-1},
        {0,1,-1},{0,-1,1}};
    
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
     if (offset > 0 && C < 0) {D-=1; }
//    if (offset > 0) offset -= 1.0;
//    D += offset;
    
    RayTrace rt;
    
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    
    /* Should read this in I guess */
    
    rt.setPlane_s(A,B,C,D);
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
    double x,y,z;
    double dat = 0;
    double a = 0;
    
    double hx = Lx/double((Nx-2)*nprocx);
    
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    imin = 0;
    imax = (Nx);
    jmin = 0;
    jmax = (Ny);
    kmin = 0;
    kmax = (Nz);
    
    double thresh = 0;
    for (int i=imin;i<imax;i++){
        for (int j=jmin;j<jmax;j++){
            for (int k=kmin;k<kmax;k++){
                m = k*Nx*Ny+j*Nx+i;
                
                /* ray origin: LBM units */
                ray_origin[0] = double(i);
                ray_origin[1] = double(j);
                ray_origin[2] = double(k);
                rt.setRayOrigin(ray_origin);
                
                /* Generate ray trace distance map */
                for (int q = 0; q < 18; q++) {
                    if (q < 6) thresh = 1; else thresh = sqrt(2);
                    
                    ray_direction[0] = D3Q19[q][0];
                    ray_direction[1] = D3Q19[q][1];
                    ray_direction[2] = D3Q19[q][2];
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_s_q();
                    if (dat <= thresh && dat > 0) qDistances[m+q*Nx*Ny*Nz] = dat; 
                }
            }
        }
    }
}


void MergeSlabQDistances(double * mergedDistances,double * qDistances1,double * qDistances2, int & Nx, int & Ny, int & Nz, double & Lx,int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank) {
    
    size_t N = Nx*Ny*Nz;
    size_t NqN = 18*Nx*Ny*Nz;
    
    
    // 1 is for the bottom plate
    // 2 is for the top plate
    for (size_t n = 0; n < NqN; n++) {
        mergedDistances[n] = std::max(qDistances1[n],qDistances2[n]);
    }
}



void ComputeIDarray(double * LIBB_id, int & Nx, int & Ny, int & Nz, double & Lx, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank,  double A, double B, double C, double D, double offset) {
    
    
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
    RayTrace rt;
    
    // Prevents ID file from being shifted up or down depending on the offset
    if (offset > 0) offset -= 1.0;
    D += offset;

    
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    
    rt.setPlane_s(A,B,C,D);
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
    double x,y,z;
    double dat = 0;
    double a = 0;
    
    double hx = Lx/double((Nx-2)*nprocx);
    
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    imin = 0;
    imax = (Nx);
    jmin = 0;
    jmax = (Ny);
    kmin = 0;
    kmax = (Nz);
    
    for (int i=imin;i<imax;i++){
        for (int j=jmin;j<jmax;j++){
            for (int k=kmin;k<kmax;k++){
                m = k*Nx*Ny+j*Nx+i;
                /* ray origin: LBM units */
                ray_origin[0] = double(i);
                ray_origin[1] = double(j);
                ray_origin[2] = double(k);
                
                rt.setRayOrigin(ray_origin);
                dat = rt.getRayTraceDistance_s();
                if (dat < LIBB_id[m]) LIBB_id[m] = dat;
            }
        }
    }
}

void MergeIDarray(char * mergedID ,double * LIBB_id1,double * LIBB_id2, int & Nx, int & Ny, int & Nz, double & Lx, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank) {
    
    double merged = 0;
    int val = 0;
    size_t N = Nx*Ny*Nz;
    for (size_t n=0; n<N; n++) {
        double id1 = LIBB_id1[n];
        double id2 = LIBB_id2[n];
        /* works for negative offsets */
        merged = std::max(id1,id2);
        if (merged <= 0.0) val = 2;
        else val = 0;
        mergedID[n] = val;
    }
}




void ComputeRayTracedDistance( double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, double A, double B, double C, double D, double offset) {
    
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
    double volume = double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
    
    RayTrace rt;
    
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    
    rt.setPlane_s(A,B,C,D);
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
    double x,y,z;
    double dat = 0;
    double a = 0;
    
    double hx = Lx/double((Nx-2)*nprocx);
    
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    imin = 0;
    imax = (Nx);
    jmin = 0;
    jmax = (Ny);
    kmin = 0;
    kmax = (Nz);
    
    for (int i=imin;i<imax;i++){
        for (int j=jmin;j<jmax;j++){
            for (int k=kmin;k<kmax;k++){
                m = k*Nx*Ny+j*Nx+i;
                /* ray origin: LBM units */
                ray_origin[0] = double(i);
                ray_origin[1] = double(j);
                ray_origin[2] = double(k);
                
                rt.setRayOrigin(ray_origin);
                dat = rt.getRayTraceDistance_s();
                if (dat < signedDistance[m]) signedDistance[m] = dat;
            }
        }
    }
}

void MergeRayTracedDistance( double * mergedRTdistance, double * rtDistance1, double * rtDistance2, int & Nx, int & Ny, int & Nz, double & Lx, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank) {
    
    size_t N = Nx*Ny*Nz;
    for (size_t n = 0; n < N; n++) {
        mergedRTdistance[n] = -std::max(rtDistance1[n],rtDistance2[n]);
    }
    
}




void ComputeMeshDistance( double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, double A, double B, double C, double D, double offset) {
    
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
    double volume = double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
    
    RayTrace rt;
    
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    
    rt.setPlane_s(A,B,C,D);
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
    double x,y,z;
    double dat = 0;
    double a = 0;
    
    double hx = Lx/double((Nx-2)*nprocx);
    
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    imin = 0;
    imax = (Nx);
    jmin = 0;
    jmax = (Ny);
    kmin = 0;
    kmax = (Nz);
    
    for (int i=imin;i<imax;i++){
        for (int j=jmin;j<jmax;j++){
            for (int k=kmin;k<kmax;k++){
                m = k*Nx*Ny+j*Nx+i;
                /* ray origin: LBM units */
                ray_origin[0] = double(i);
                ray_origin[1] = double(j);
                ray_origin[2] = double(k+0.5);
                
                rt.setRayOrigin(ray_origin);
                dat = rt.getRayTraceDistance_s();
                if (dat < signedDistance[m]) signedDistance[m] = dat;
            }
        }
    }
}


void MergeMeshDistance( double * mergedRTdistance, double * rtDistance1, double * rtDistance2, int & Nx, int & Ny, int & Nz, double & Lx, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank) {
    
    size_t N = Nx*Ny*Nz;
    for (size_t n = 0; n < N; n++) {
        mergedRTdistance[n] = -std::max(rtDistance1[n],rtDistance2[n]);
    }
    
}







double RayTrace::dist_sp(double & term1 /* "u.(o-c)" */,
                         double & term2 /* discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )" */) {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    
    double d1,d2,pos;
    d1 = 0;
    d2 = 0;
    pos = 0;
    /* No solution exists */
    if (term2 < 0) pos = 0;
    
    /* One solution exists */
    if (term2 == 0) {
        d1 = -term1;
        pos = d1;
    }
    
    /* Two solutions exist */
    if (term2 > 0) {
        d1 = -term1 + sqrt(term2);
        d2 = -term1 - sqrt(term2);
        if (d2 >= d1) std::swap(d2,d1);
        if (d2 < 0) d2 = d1;
        pos = d2;
    }
    
    return pos;
}




double RayTrace::ray_trace_sp(double *  vec1 /* ray_origin */,
                              double *  vec2 /*sphere_center */,
                              double *  vec3 /* ray direction -> will become unit vector */,
                              double rad /* analytical sphere radius */) {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // Note: using "u" for ComputeLIBBqDistances_spunit vector instead of "l"
    
    /* Ray tracing work */
    subtract_vectors(vec1 /* ray origin "o" */,
                     vec2 /* sphere center "c" */,
                     o_minus_c /* Return: "o-c"*/);
    
    /* Normalize "unit vector" */
    double tmp = mag(vec3);
    if (tmp == 0.0) tmp = 1.0; // Needed for domain (LIBB_id) generation
    vec3[0] /= tmp; vec3[1] /= tmp; vec3[2] /= tmp;
    dp = dot(vec3 /* unit vector */ ,o_minus_c); // Return: "u.(o-c)"
    under_sqrt = dp*dp - (magSq(o_minus_c) - rad*rad); // Return the discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )"
    /* Get the solution */
    return dist_sp(dp, under_sqrt);
}


void ReadSpherePacking(int nspheres,
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

void ComputeLIBBqDistances_sp(char * LIBB_id,
                              double * qDistances, double * signedDistance, double * signedDistance2, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm) {
    
    
    
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
    double volume = double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
    
    RayTrace rt;
    double * sphere_center;  sphere_center = new double[3];
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    double sphere_radius;
    
    size_t cross_count = 0;
    double *ccx,*ccy,*ccz,*r;
    ccx = new double[nspheres];
    ccy = new double[nspheres];
    ccz = new double[nspheres];
    r = new double[nspheres];
    ReadSpherePacking(nspheres,
                      ccx,ccy,ccz,r);
    double cx,cy,cz,rad;
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
    double x,y,z;
    
    double dat = 0;
    double a = 0;
    
    //  cout << "Lx=" << Lx << endl;
    double hx = Lx/double((Nx-2)*nprocx);
    
    double sr = r[0];
    // if (rank == 0) cout << "All sphere radii are " << sr << " lattice units" << endl;
    
    double analytical_porosity = (1.0-double(2)* (4./3.*3.14159265359 * sr*sr*sr))/1.0;
    // if (rank == 0) printf("Analytical porosity=%.15f (valid for equally sized spheres)\n",analytical_porosity);
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    /* GENERATE CORRECT ID */
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5;
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;
        
        /* Check for range of a particular sphere */
        imin = int ((cx-rad)/hx)-2;
        imax = int ((cx+rad)/hx)+3;
        jmin = int ((cy-rad)/hx)-2;
        jmax = int ((cy+rad)/hx)+3;
        kmin = int ((cz-rad)/hx)-2;
        kmax = int ((cz+rad)/hx)+3;
        if (imin<0)     imin = 0;
        if (imin>Nx)    imin = Nx;
        if (imax<0)     imax = 0;
        if (imax>Nx)    imax = Nx;
        if (jmin<0)     jmin = 0;
        if (jmin>Ny)    jmin = Ny;
        if (jmax<0)     jmax = 0;
        if (jmax>Ny)    jmax = Ny;
        if (kmin<0)     kmin = 0;
        if (kmin>Nz)    kmin = Nz;
        if (kmax<0)     kmax = 0;
        if (kmax>Nz)    kmax = Nz;
        
        rt.setSphereCenter(sphere_center);
        rt.setSphereRadius(sphere_radius);
        
        // Loop over the domain for this sphere (may be null)
        
        for (int i=imin;i<imax;i++){
            for (int j=jmin;j<jmax;j++){
                for (int k=kmin;k<kmax;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    /* ray origin: LBM units */
                    ray_origin[0] = double(i);
                    ray_origin[1] = double(j);
                    ray_origin[2] = double(k);
                    
                    rt.setRayOrigin(ray_origin);
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat > 0) LIBB_id[m] = 0;
                    
                }
            }
        }
    }
    
    
    
    
    
    int rank_size = nprocx*nprocy*nprocz;
    int * rankArray; rankArray = new int[rank_size];
    for (int i = 0; i < rank_size; i++) rankArray[i] = 0;
    
    vector<int> whichSphere;
    
    // my rank
    int myRank = iproc + jproc*nprocx + kproc*nprocx*nprocy;
    
    int rankspheres = 0;
    int negCount = 0;
    for (int p=0;p<nspheres;p++) {
        rankspheres = 0;
        
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5+0.5;
        sphere_center[2] = cz/hx+0.5+0.5;
        sphere_radius = rad/hx;
        
        negCount = 0;
        for (int i=1;i<Nx-1;i++) {
            for (int j=1;j<Ny-1;j++) {
                for (int k=1;k<Nz-1;k++) {
                    m = k*Nx*Ny+j*Nx+i;
                    
                    double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0])
                                + (double(j)-sphere_center[1])*(double(j)-sphere_center[1])
                                + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) )
                                - sphere_radius;
                    
                    if (dist < signedDistance[m]) signedDistance[m] = dist;
                    if (dist < 0) negCount++;
                }
            }
        }
        
        if (negCount > 0) {
            rankArray[myRank] = rankArray[myRank] + 1;
            whichSphere.push_back(p);
        }
    } // n spheres loop
    
    

    MPI_Barrier(Dm.Comm);
    MPI_Allreduce(MPI_IN_PLACE, rankArray, rank_size, MPI_INT, MPI_SUM, Dm.Comm);
    MPI_Barrier(Dm.Comm);
    
  int sumSpheres = 0;
    for (int a = 0; a < rank_size; a++)  {
        
        sumSpheres += rankArray[a];
        
    }
    
    if (sumSpheres == 0) {
        for (int a = 0; a < rank_size; a++)  {
            
            rankArray[a] = 1;
        }
    }
    
    
    // Allocate SignedDistance field with multiple levels given the number of spheres in each rank
    double * RankRayTraceDistance; RankRayTraceDistance = new double[rankArray[myRank]*Nx*Ny*Nz];
    for (int i = 0; i < rankArray[myRank]*Nx*Ny*Nz; i++) RankRayTraceDistance[i] = 1000;
    
    // should not be triggered if no spheres present
    for (int s=0; s < whichSphere.size(); s++) {
        int p = whichSphere.at(s);
        
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;
        
        for (int i=0;i<Nx;i++) {
            for (int j=0;j<Ny;j++) {
                for (int k=0;k<Nz;k++) {
                    m = k*Nx*Ny+j*Nx+i;
                    double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0])
                                + (double(j)-sphere_center[1])*(double(j)-sphere_center[1])
                                + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) )
                                - sphere_radius;
                    
                    if (dist < RankRayTraceDistance[m+s*(Nx*Ny*Nz)]) {
                        RankRayTraceDistance[m+s*(Nx*Ny*Nz)] = dist;
                    }
                }
            }
        }
        
        
        
    }
    
    
    if (sumSpheres == 0) whichSphere.push_back(0); // resizes the whichSphere to 1
    
       
    char LocalRankString[8]; sprintf(LocalRankString,"%05d",rank);
    char LocalRankFilename1[40];
    
    sprintf(LocalRankFilename1,"%s%s","RankRTD.",LocalRankString);
    FILE *RTFILE = fopen(LocalRankFilename1,"wb");
    if (RTFILE==NULL) ERROR("Error opening file: rankRT.xxxxx");
    fwrite(RankRayTraceDistance,8,(Nx*Ny*Nz)*whichSphere.size(),RTFILE);
    fclose(RTFILE);
    
    
    

    
    if (rank == 0) {
        char LocalRankString[8]; sprintf(LocalRankString,"%05d",rank);
        char LocalRankFilename1[40];
        
        // how many spheres on a rank
        sprintf(LocalRankFilename1,"%s%s","NSpheresPerRankArray.",LocalRankString);
        FILE *NSPRA = fopen(LocalRankFilename1,"wb");
        if (NSPRA==NULL) ERROR("Error opening file: NSPRA.xxxxx");
        fwrite(rankArray,4,rank_size,NSPRA);
        fclose(NSPRA);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // Output for mesh
    for (int p=0;p<nspheres;p++) {
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){
                for (int k=1;k<Nz-1;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
                    if (dist < signedDistance2[m]) signedDistance2[m] = dist;
                }
            }
        }
    } // n spheres loop
    
    
    
    
    
    
    // .........Loop over the spheres.............
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5;
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;
        
        /* Check for range of a particular sphere */
        imin = int ((cx-rad)/hx)-2;
        imax = int ((cx+rad)/hx)+3;
        jmin = int ((cy-rad)/hx)-2;
        jmax = int ((cy+rad)/hx)+3;
        kmin = int ((cz-rad)/hx)-2;
        kmax = int ((cz+rad)/hx)+3;
        if (imin<0)     imin = 0;
        if (imin>Nx)    imin = Nx;
        if (imax<0)     imax = 0;
        if (imax>Nx)    imax = Nx;
        if (jmin<0)     jmin = 0;
        if (jmin>Ny)    jmin = Ny;
        if (jmax<0)     jmax = 0;
        if (jmax>Ny)    jmax = Ny;
        if (kmin<0)     kmin = 0;
        if (kmin>Nz)    kmin = Nz;
        if (kmax<0)     kmax = 0;
        if (kmax>Nz)    kmax = Nz;
        
        rt.setSphereCenter(sphere_center);
        rt.setSphereRadius(sphere_radius);
        
        
        double D3Q19[18][3] = {{1,0,0},{-1,0,0},
               {0,1,0},{0,-1,0},
               {0,0,1},{0,0,-1},
               {1,1,0},{-1,-1,0},
               {1,-1,0},{-1,1,0},
               {1,0,1},{-1,0,-1},
               {1,0,-1},{-1,0,1},
               {0,1,1},{0,-1,-1},
               {0,1,-1},{0,-1,1}};
        
        double thresh = 0;
        // Loop over the domain for this sphere (may be null)
        for (int i=imin;i<imax;i++){
            for (int j=jmin;j<jmax;j++){
                for (int k=kmin;k<kmax;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    
                    /* ray origin: LBM units */
                    ray_origin[0] = double(i);
                    ray_origin[1] = double(j);
                    ray_origin[2] = double(k);
                    rt.setRayOrigin(ray_origin);
                    
                    /* Generate ray trace distance map */
                    for (int q = 0; q < 18; q++) {
                        if (q < 6) thresh = 1; else thresh = sqrt(2);
                        
                        ray_direction[0] = D3Q19[q][0];
                        ray_direction[1] = D3Q19[q][1];
                        ray_direction[2] = D3Q19[q][2];
                        rt.setRayDirection(ray_direction);
                        dat = rt.getRayTraceDistance_sp();
                        if (dat <= thresh && dat > 0) qDistances[m+q*Nx*Ny*Nz] = dat;
                    }
                }
            }
        }
    }
}
