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

#include "common/ScaLBL.h"

using namespace std;

class RayTrace {
public:
    RayTrace() {
        ray_origin = new double[3];
        ray_direction = new double[3];
        geometry = 0;
        
        /* s: geometry = 2 */
        plane_s = new double[3];
        /* */
        plane_pp1 = new double[3];
        plane_pp2 = new double[3];
        
        
        plane_st1 = new double[3];
        plane_st2 = new double[3];
        plane_st3 = new double[3];
        plane_st4 = new double[3];
        
        
        
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
    double * getRayDirection_sp() { return ray_direction; }
    double getRayTraceDistance_sp() {
        return ray_trace_sp(ray_origin, sphere_center, ray_direction, sphere_radius);
    }
    
    void setPlane_s(double A, double B, double C, double D) {
        plane_s[0] = A;
        plane_s[1] = B;
        plane_s[2] = C;
        D_s = D;
    }
    
    
    void setPlane_pp1(double A, double B, double C, double D) {
        plane_pp1[0] = A;
        plane_pp1[1] = B;
        plane_pp1[2] = C;
        D_pp1 = D;
    }
    
    void setPlane_pp2(double A, double B, double C, double D) {
        plane_pp2[0] = A;
        plane_pp2[1] = B;
        plane_pp2[2] = C;
        D_pp2 = D;
    }
    
    
    
    void setPlane_st1(double A, double B, double C, double D) {
        plane_st1[0] = A;
        plane_st1[1] = B;
        plane_st1[2] = C;
        D_st1 = D;
    }
    void setPlane_st2(double A, double B, double C, double D) {
        plane_st2[0] = A;
        plane_st2[1] = B;
        plane_st2[2] = C;
        D_st2 = D;
    }
    void setPlane_st3(double A, double B, double C, double D) {
        plane_st3[0] = A;
        plane_st3[1] = B;
        plane_st3[2] = C;
        D_st3 = D;
    }
    void setPlane_st4(double A, double B, double C, double D) {
        plane_st4[0] = A;
        plane_st4[1] = B;
        plane_st4[2] = C;
        D_st4 = D;
    }
    
    
    
    double * getRayDirection_s() { return ray_direction; }
    double getRayTraceDistance_s() {
        return ray_trace_s(ray_origin, plane_s, D_s);
    }
    
    double getRayTraceDistance_s_q() {
        return ray_trace_s_q(ray_origin, plane_s, D_s);
    }
    
    double * getRayDirection_pp() { return ray_direction; }
    double getRayTraceDistance_pp() {
        return ray_trace_pp(ray_origin, plane_pp1, D_pp1, plane_pp2, D_pp2);
    }
    
    double getRayTraceDistance_pp_q() {
        return ray_trace_pp_q(ray_origin, ray_direction, plane_pp1, D_pp1, plane_pp2, D_pp2);
    }
    
    
    double * getRayDirection_st() { return ray_direction; }
    double getRayTraceDistance_st() {
        return ray_trace_st(ray_origin,
                            plane_st1, D_st1,
                            plane_st2, D_st2,
                            plane_st3, D_st3,
                            plane_st4, D_st4);
    }
    
    double getRayTraceDistance_st_q() {
        return ray_trace_st_q(ray_origin,
                              plane_st1, D_st1,
                              plane_st2, D_st2,
                              plane_st3, D_st3,
                              plane_st4, D_st4);
    }
    
    double * getRayDirection_ct() { return ray_direction; }
    double getRayTraceDistance_ct() {
        return ray_trace_ct(ray_origin, sphere_center, ray_direction, sphere_radius);
    }
    
private:
    // 1 periodic sphere pack, 2 slab, 3 parallel plates, 4 square tube, 5 capillary tube
    int geometry;
    
    /* Operators needed to compute the ray_trace function */
    double dot(double *  vec1, double *  vec2);
    double mag(double *  vec1 /* Get || . ||  of this vector */);
    double magSq(double *  vec1 /* Get || . ||^2 of this vector */);
    void subtract_vectors(double *  vec1 /* ray_origin */, double *  vec2 /* cell_center */,double *  vec3 /* Return vector */);
    
    
    
    /* Signatures */
    double dist_sp(double & term1 /* u_dot_o_minus_c */, double & term2 /* under_sqrt */);
    
    
    
    double dist_ct(double & term1 /* "u.(o-c)" */,
                   double & term2 /* discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )" */);
    
    
    double ray_trace_sp(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */);
    
    double ray_trace_s(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* plane translation */);
    double ray_trace_s_q(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* plane translation */);
    
    double ray_trace_pp(double *  vec1 /* ray_origin */, double *  vec2 /* plane_pp1 */, double D_pp1, double *  vec3 /* plane_pp2 */, double D_pp2);
    double ray_trace_pp_q(double *  lnaught /* ray_origin */, double * l /* ray direction */,
                          double *  plane_pp1 /* plane_pp1 */, double D_pp1 /* distance from origin- need to generalize this for phi,theta!=0 */,
                          double *  plane_pp2 /* plane_pp2 */, double D_pp2);
    
    double ray_trace_st(double *  vec1 /* ray_origin */,
                        double *  vec2 /* plane_pp1 */, double D_st1,
                        double *  vec3 /* plane_pp2 */, double D_st2,
                        double *  vec4 /* plane_st3 */, double D_st3,
                        double *  vec5 /* plane_st4 */, double D_st4);
    
    double ray_trace_st_q(double *  vec1 /* ray_origin */,
                          double *  vec2 /* plane_pp1 */, double D_st1,
                          double *  vec3 /* plane_pp2 */, double D_st2,
                          double *  vec4 /* plane_st3 */, double D_st3,
                          double *  vec5 /* plane_st4 */, double D_st4);
    
    double ray_trace_ct(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */);
    
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
    
    /* parallel plates */
    double * plane_pp1, * plane_pp2;
    double D_pp1,D_pp2;
    
    /* square tube */
    double * plane_st1, * plane_st2, * plane_st3, * plane_st4;
    double D_st1,D_st2, D_st3, D_st4;
    
    /* capillary tube */
    // NYI
    
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



double RayTrace::dist_ct(double & term1 /* "u.(o-c)" */,
                         double & term2 /* discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )" */) { return 0.0; }





double RayTrace::ray_trace_sp(double *  vec1 /* ray_origin */,
                              double *  vec2 /*sphere_center */,
                              double *  vec3 /* ray direction -> will become unit vector */,
                              double rad /* analytical sphere radius */) {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // Note: using "u" for unit vector instead of "l"
    
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

double RayTrace::ray_trace_s(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* displacement */) {
    
    double tmp = mag(vec2); // Denominator of distance
    if (tmp == 0) tmp = 1.0;
    
    
    dp = dot(vec2, vec1);
    // cout << "dp=" << dp << endl;
    
    return (-(dp - D_s)/tmp);
}


double RayTrace::ray_trace_s_q(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* displacement */)
{
    double tmp = mag(vec2); if (tmp == 0) tmp = 1.0;
    subtract_vectors(vec2 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    o_minus_c[0] += 0; o_minus_c[1] += 0; o_minus_c[2] += 0;
    
    double * normal; normal = new double[3];
    normal[0] = vec2[0]/=tmp;
    normal[1] = vec2[1]/=tmp;
    normal[2] = vec2[2]/=tmp;
    
    tmp = dot(ray_direction,normal); if (tmp == 0) tmp = 1.0;
    dp = dot(o_minus_c, normal) / tmp;
    
    delete[] normal;
    
    // if (dp < 1.5) cout << "dp=" << dp << endl;
    
    return (-(dp-(D_s)));
}





//double RayTrace::ray_trace_st(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */) { return 0.0; }
//
//double RayTrace::ray_trace_ct(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */) { return 0.0; }


static inline void fgetl( char * str,
                         int num, FILE * stream ) {
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
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

void SetIdAndCount(char * id,
                   size_t & count, double & hx, int & Nx, int & Ny, int & Nz, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz,int rank) {} // Not used

void ComputeLIBBqDistances_sp(char * LIBB_id,
                              double * qDistances, double * signedDistance, double * signedDistance2, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank) {
    
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
    if (rank == 0) cout << "All sphere radii are " << sr << " lattice units" << endl;
    
    double analytical_porosity = (1.0-double(2)* (4./3.*3.14159265359 * sr*sr*sr))/1.0;
    if (rank == 0) printf("Analytical porosity=%.15f (valid for equally sized spheres)\n",analytical_porosity);
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    //define for later
    int dp = nspheres;
    int di;
    int dj;
    int dk;
    int idx_max;
    
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
        
        di = imax-imin;
        dj = jmax-jmin;
        dk = kmax-kmin;
        idx_max = di*dj*dk;
        for (int idx=0;idx<idx_max;idx++)
        {
            int i = idx / (dj*dk);
            int j = (idx - i*dj*dk) / dk;
            int k = idx - j*dk - i*dj*dk;

            i += imin;
            j += jmin;
            k += kmin;

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
        }//end of idx loop
    }
    
    di = Nx-2;
    dj = Ny-2;
    dk = Nz-2;
    idx_max = di*dj*dk;
    for (int p=0;p<nspheres;p++)
    {
        for (int idx=0;idx<idx_max;idx++)
        {
            int i = idx / (dj*dk);
            int j = (idx - i*dj*dk) / dk;
            int k = idx - j*dk - i*dj*dk;

            i += 1; //imin=1
            j += 1; //jmin=1
            k += 1; //kmin=1
            
            cx = ccx[p] - min_x;
            cy = ccy[p] - min_y;
            cz = ccz[p] - min_z;
            rad = r[p];

            /* converting sphere center to LBM units */
            sphere_center[0] = cx/hx+0.5+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
            sphere_center[1] = cy/hx+0.5+0.5;
            sphere_center[2] = cz/hx+0.5+0.5;
            sphere_radius = rad/hx;

            m = k*Nx*Ny+j*Nx+i;
            double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
            if (dist < signedDistance[m]) signedDistance[m] = dist;
        }//end of idx loop
    }
    
    //idx_max, dp, di, dj, dk same as before
    for (int p=0;p<nspheres;p++)
    {
        for (int idx=0;idx<idx_max;idx++)
        {
            int i = idx / (dj*dk);
            int j = (idx - i*dj*dk) / dk;
            int k = idx - j*dk - i*dj*dk;

            i += 1; //imin=1
            j += 1; //jmin=1
            k += 1; //kmin=1

            cx = ccx[p] - min_x;
            cy = ccy[p] - min_y;
            cz = ccz[p] - min_z;
            rad = r[p];

            /* converting sphere center to LBM units */
            sphere_center[0] = cx/hx+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
            sphere_center[1] = cy/hx+0.5;
            sphere_center[2] = cz/hx+0.5;
            sphere_radius = rad/hx;

            m = k*Nx*Ny+j*Nx+i;
            double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
            if (dist < signedDistance2[m]) signedDistance2[m] = dist;
        }
    }
    
    
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
        
        di = imax-imin;
        dj = jmax-jmin;
        dk = kmax-kmin;
        idx_max = di*dj*dk;
        for (int idx=0;idx<idx_max;idx++){
            int i = idx / (dj*dk);
            int j = (idx - i*dj*dk) / dk;
            int k = idx - j*dk - i*dj*dk;

            i += imin;
            j += jmin;
            k += kmin;
            
            m = k*Nx*Ny+j*Nx+i;

            /* ray origin: LBM units */
            ray_origin[0] = double(i);
            ray_origin[1] = double(j);
            ray_origin[2] = double(k);
            rt.setRayOrigin(ray_origin);

            /* Generate ray trace distance map */
            // 1
            ray_direction[0] = 1.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= 1 && dat > 0) qDistances[m] = dat;

            // 2
            ray_direction[0] = -1.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= 1 && dat > 0) qDistances[m+Nx*Ny*Nz] = dat;

            // 3
            ray_direction[0] = 0.0;
            ray_direction[1] = 1.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= 1 && dat > 0) qDistances[m+2*Nx*Ny*Nz] = dat;

            // 4
            ray_direction[0] = 0.0;
            ray_direction[1] = -1.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= 1 && dat > 0) qDistances[m+3*Nx*Ny*Nz] = dat;

            // 5
            ray_direction[0] = 0.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = 1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= 1 && dat > 0) qDistances[m+4*Nx*Ny*Nz] = dat;

            // 6
            ray_direction[0] = 0.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = -1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= 1 && dat > 0) qDistances[m+5*Nx*Ny*Nz] = dat;

            // 7
            ray_direction[0] = 1.0;
            ray_direction[1] = 1.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+6*Nx*Ny*Nz] = dat;

            // 8
            ray_direction[0] = -1.0;
            ray_direction[1] = -1.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+7*Nx*Ny*Nz] = dat;

            // 9
            ray_direction[0] = 1.0;
            ray_direction[1] = -1.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+8*Nx*Ny*Nz] = dat;

            // 10
            ray_direction[0] = -1.0;
            ray_direction[1] = 1.0;
            ray_direction[2] = 0.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+9*Nx*Ny*Nz] = dat;

            // 11
            ray_direction[0] = 1.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = 1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+10*Nx*Ny*Nz] = dat;

            // 12
            ray_direction[0] = -1.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = -1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+11*Nx*Ny*Nz] = dat;

            // 13
            ray_direction[0] = 1.0;
            ray_direction[1] = 0.0;
            ray_direction[2] = -1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+12*Nx*Ny*Nz] = dat;

            // 14
            ray_direction[0] = -1;
            ray_direction[1] = 0.0;
            ray_direction[2] = 1;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+13*Nx*Ny*Nz] = dat;

            // 15
            ray_direction[0] = 0.0;
            ray_direction[1] = 1;
            ray_direction[2] = 1;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+14*Nx*Ny*Nz] = dat;

            // 16
            ray_direction[0] = 0.0;
            ray_direction[1] = -1;
            ray_direction[2] = -1;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            dat = abs(dat);
            if (dat <= sqrt(2) && dat > 0) qDistances[m+15*Nx*Ny*Nz] = dat;

            // 17
            ray_direction[0] = 0.0;
            ray_direction[1] = 1.0;
            ray_direction[2] = -1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+16*Nx*Ny*Nz] = dat;

            // 18
            ray_direction[0] = 0.0;
            ray_direction[1] = -1.0;
            ray_direction[2] = 1.0;
            rt.setRayDirection(ray_direction);
            dat = rt.getRayTraceDistance_sp();
            if (dat <= sqrt(2) && dat > 0) qDistances[m+17*Nx*Ny*Nz] = dat;
        }//idx loop
    }//nspheres loop
}

void run_sp_raytrace(/*arguments for ComputeLIBBqDistances_sp*/
                     char * LIBB_id,
                     double * qDistances,
                     double * signedDistance,
                     double * signedDistance2,
                     int & Nx, int & Ny, int & Nz,
                     double & Lx,
                     int & nspheres,
                     int iproc, int jproc, int kproc,
                     int nprocx, int nprocy, int nprocz,
                     int rank,
                     /*extra arguments*/
                     size_t N, size_t NqN)
{
    char* d_LIBB_id;
    double* d_qDistances;
    double* d_signedDistance;
    double* d_signedDistance2;

    //allocate memory on device
    ScaLBL_AllocateDeviceMemory((void**)&d_LIBB_id, sizeof(char)*N);
    ScaLBL_AllocateDeviceMemory((void**)&d_qDistances, sizeof(double)*NqN);
    ScaLBL_AllocateDeviceMemory((void**)&d_signedDistance, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void**)&d_signedDistance2, sizeof(double)*N);

    //copy memory to device
    ScaLBL_CopyToDevice(d_LIBB_id,LIBB_id, sizeof(char)*N);
    ScaLBL_CopyToDevice(d_qDistances,qDistances, sizeof(double)*NqN);
    ScaLBL_CopyToDevice(d_signedDistance,signedDistance, sizeof(double)*N);
    ScaLBL_CopyToDevice(d_signedDistance2,signedDistance2, sizeof(double)*N);

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
    
    //allocate and copy some more arrays
    double *d_sphere_center, *d_ccx, *d_ccy, *d_ccz, *d_r;
//    ScaLBL_AllocateDeviceMemory((void**)&d_sphere_center, sizeof(double)*3);
    ScaLBL_AllocateDeviceMemory((void**)&d_ccx, sizeof(double)*nspheres);
    ScaLBL_AllocateDeviceMemory((void**)&d_ccy, sizeof(double)*nspheres);
    ScaLBL_AllocateDeviceMemory((void**)&d_ccz, sizeof(double)*nspheres);
    ScaLBL_AllocateDeviceMemory((void**)&d_r, sizeof(double)*nspheres);
    
//    ScaLBL_CopyToDevice(d_sphere_center, sphere_center, sizeof(double)*3);
    ScaLBL_CopyToDevice(d_ccx, ccx, sizeof(double)*nspheres);
    ScaLBL_CopyToDevice(d_ccy, ccy, sizeof(double)*nspheres);
    ScaLBL_CopyToDevice(d_ccz, ccz, sizeof(double)*nspheres);
    ScaLBL_CopyToDevice(d_r, r, sizeof(double)*nspheres);
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
    double x,y,z;
    
    double dat = 0;
    double a = 0;
    
    //  cout << "Lx=" << Lx << endl;
    double hx = Lx/double((Nx-2)*nprocx);
    
    double sr = r[0];
    if (rank == 0) cout << "All sphere radii are " << sr << " lattice units" << endl;
    
    double analytical_porosity = (1.0-double(2)* (4./3.*3.14159265359 * sr*sr*sr))/1.0;
    if (rank == 0) printf("Analytical porosity=%.15f (valid for equally sized spheres)\n",analytical_porosity);
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    double sphere_center_x, sphere_center_y, sphere_center_z;
    
    /* GENERATE CORRECT ID */
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
//        sphere_center[0] = cx/hx+0.5;
//        sphere_center[1] = cy/hx+0.5;
//        sphere_center[2] = cz/hx+0.5;
        sphere_center_x = cx/hx+0.5;
        sphere_center_y = cy/hx+0.5;
        sphere_center_z = cz/hx+0.5;
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
        
//        rt.setSphereCenter(sphere_center);
//        rt.setSphereRadius(sphere_radius);
        
        IDFunction(imax, imin, jmax, jmin, kmax, kmin, Nx, Ny, Nz, sphere_radius, sphere_center_x, sphere_center_y, sphere_center_z, d_LIBB_id);
    }
    
    for (int p=0;p<nspheres;p++)
    {
        SignedDistanceFunction(p, hx, Nx, Ny, Nz, min_x, min_y, min_z, d_ccx, d_ccy, d_ccz, d_r, d_signedDistance);
    }
    
    for (int p=0;p<nspheres;p++)
    {
        SignedDistanceFunction2(p, hx, Nx, Ny, Nz, min_x, min_y, min_z, d_ccx, d_ccy, d_ccz, d_r, d_signedDistance2);
    }
     
    // .........Loop over the spheres.............
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
//        sphere_center[0] = cx/hx+0.5;
//        sphere_center[1] = cy/hx+0.5;
//        sphere_center[2] = cz/hx+0.5;
        sphere_center_x = cx/hx+0.5;
        sphere_center_y = cy/hx+0.5;
        sphere_center_z = cz/hx+0.5;
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
        
//        rt.setSphereCenter(sphere_center);
//        rt.setSphereRadius(sphere_radius);
        
        QDistancesFunction(imax, imin, jmax, jmin, kmax, kmin, Nx, Ny, Nz, sphere_radius, sphere_center_x, sphere_center_y, sphere_center_z, d_qDistances);
    }//nspheres loop

    //copy memory back to host
    ScaLBL_CopyToHost(LIBB_id,d_LIBB_id, sizeof(char)*N);
    ScaLBL_CopyToHost(qDistances,d_qDistances, sizeof(double)*NqN);
    ScaLBL_CopyToHost(signedDistance,d_signedDistance, sizeof(double)*N);
    ScaLBL_CopyToHost(signedDistance2,d_signedDistance2, sizeof(double)*N);
    
//    ScaLBL_CopyToHost(sphere_center, d_sphere_center, sizeof(double)*3);
    ScaLBL_CopyToHost(ccx, d_ccx, sizeof(double)*nspheres);
    ScaLBL_CopyToHost(ccy, d_ccy, sizeof(double)*nspheres);
    ScaLBL_CopyToHost(ccz, d_ccz, sizeof(double)*nspheres);
    ScaLBL_CopyToHost(r, d_r, sizeof(double)*nspheres);
}

void ComputeLIBBqDistances_s(char * LIBB_id,
                             double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm) {
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
    double volume = double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
    
    RayTrace rt;
    
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    
    /* Should read this in I guess */
    double A,B,C,D;
    A = 0;
    B = 0;
    C = 1;
    D = 16.75;  // intercept
    cout << "set plane_s\n" << endl;
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
    
    cout << " calculating distance for ID field..." << endl;
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
                dat = rt.getRayTraceDistance_s();
                if (dat > 0) LIBB_id[m] = 0;
                if (dat < signedDistance[m]) signedDistance[m] = dat;
            }
        }
    }
    
    
    /* STEP 3: q distances */
    /* Check for range of a particular sphere */
    //        imin = int ((cx-rad)/hx)-2;
    //        imax = int ((cx+rad)/hx)+3;
    //        jmin = int ((cy-rad)/hx)-2;
    //        jmax = int ((cy+rad)/hx)+3;
    //        kmin = int ((cz-rad)/hx)-2;
    //        kmax = int ((cz+rad)/hx)+3;
    //        if (imin<0)     imin = 0;
    //        if (imin>Nx)    imin = Nx;
    //        if (imax<0)     imax = 0;
    //        if (imax>Nx)    imax = Nx;
    //        if (jmin<0)     jmin = 0;
    //        if (jmin>Ny)    jmin = Ny;
    //        if (jmax<0)     jmax = 0;
    //        if (jmax>Ny)    jmax = Ny;
    //        if (kmin<0)     kmin = 0;
    //        if (kmin>Nz)    kmin = Nz;
    //        if (kmax<0)     kmax = 0;
    //        if (kmax>Nz)    kmax = Nz;
    
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
                // 1
                ray_direction[0] = 1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= 1 && dat > 0) qDistances[m] = dat;
                
                // 2
                ray_direction[0] = -1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= 1 && dat > 0) qDistances[m+Nx*Ny*Nz] = dat;
                
                // 3
                ray_direction[0] = 0.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= 1 && dat > 0) qDistances[m+2*Nx*Ny*Nz] = dat;
                
                // 4
                ray_direction[0] = 0.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= 1 && dat > 0) qDistances[m+3*Nx*Ny*Nz] = dat;
                
                // 5
                ray_direction[0] = 0.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= 1 && dat > 0) qDistances[m+4*Nx*Ny*Nz] = dat;
                
                // 6
                ray_direction[0] = 0.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= 1 && dat > 0) qDistances[m+5*Nx*Ny*Nz] = dat;
                
                // 7
                ray_direction[0] = 1.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+6*Nx*Ny*Nz] = dat;
                
                // 8
                ray_direction[0] = -1.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+7*Nx*Ny*Nz] = dat;
                
                // 9
                ray_direction[0] = 1.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+8*Nx*Ny*Nz] = dat;
                
                // 10
                ray_direction[0] = -1.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+9*Nx*Ny*Nz] = dat;
                
                // 11
                ray_direction[0] = 1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+10*Nx*Ny*Nz] = dat;
                
                // 12
                ray_direction[0] = -1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+11*Nx*Ny*Nz] = dat;
                
                // 13
                ray_direction[0] = 1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+12*Nx*Ny*Nz] = dat;
                
                // 14
                ray_direction[0] = -1;
                ray_direction[1] = 0.0;
                ray_direction[2] = 1;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+13*Nx*Ny*Nz] = dat;
                
                // 15
                ray_direction[0] = 0.0;
                ray_direction[1] = 1;
                ray_direction[2] = 1;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+14*Nx*Ny*Nz] = dat;
                
                // 16
                ray_direction[0] = 0.0;
                ray_direction[1] = -1;
                ray_direction[2] = -1;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                dat = abs(dat);
                if (dat <= sqrt(2) && dat > 0) qDistances[m+15*Nx*Ny*Nz] = dat;
                
                // 17
                ray_direction[0] = 0.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+16*Nx*Ny*Nz] = dat;
                
                // 18
                ray_direction[0] = 0.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_s_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+17*Nx*Ny*Nz] = dat;
            }
        }
    }
    
}




double RayTrace::ray_trace_pp(double *  vec1 /* ray_origin */,
                              double *  vec2 /* plane_pp1 */,
                              double D_pp1,
                              double *  vec3 /* plane_pp2 */,
                              double D_pp2) {
    
    double tmp1 = mag(vec2); if (tmp1 == 0) tmp1 = 1.0;
    double tmp2 = mag(vec3); if (tmp2 == 0) tmp2 = 1.0;
    
    dp1 = dot(vec1, vec2); dp2 = dot(vec1, vec3);
    
    double distance1 = -(dp1 - D_pp1)/tmp1; double distance2 = -(dp2 + D_pp2)/tmp2;
    
    return std::max(distance1,distance2);
}


/*  Ray-Plane Intersection  https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
 */
double RayTrace::ray_trace_pp_q(double *  lnaught /* ray_origin */, double * l /* ray direction */,
                                double *  plane_pp1 /* plane_pp1 */, double D_pp1 /* distance from origin- need to generalize this for phi,theta!=0 */,
                                double *  plane_pp2 /* plane_pp2 */, double D_pp2) {
    
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
    
    tmp = mag(plane_pp1); if (tmp == 0) tmp = 1.0; // used if plane is rotated
    /* normal vector to the plane pp1 */
    normal[0] = plane_pp1[0]/tmp;
    normal[1] = plane_pp1[1]/tmp;
    normal[2] = plane_pp1[2]/tmp;
    
    
    /* distance point from world origin (needs to be parallelized!) */
    p0x = 0;  p0y = 0;  p0z = D_pp1;
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
    dp1 = numerator / denominator;  // Result for plane 1
    if (dp1 < 0 ) dp1 = 9;
    
    
    tmp = mag(plane_pp2); if (tmp == 0) tmp = 1.0; // used if plane is rotated
    /* normal vector to the plane pp1 */
    normal[0] = plane_pp2[0]/tmp;
    normal[1] = plane_pp2[1]/tmp;
    normal[2] = plane_pp2[2]/tmp;
    
    
    /* distance point from world origin (needs to be parallelized!) */
    p0x = 0;  p0y = 0;  p0z = D_pp2;
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
    dp2 = numerator / denominator;  // Result for plane 1
    if (dp2 < 0 ) dp2 = 9;
    
    
    
    // std::cout << "dp1=" << dp1 << " dp2=" << dp2 << std::endl;
    double lg = 1000;
    delete[] pnaught;
    delete[] normal;
    return std::min(dp2,dp1);
}


void ComputeLIBBqDistances_pp(char * LIBB_id,
                              double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm, double offset)
{
    size_t m;
    size_t count = 0;
    size_t c = 0;
    
    double volume = double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
    
    RayTrace rt;
    
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    
    bool symmetric_walls = false;
    
    /* Should read this in I guess */
    double A1,B1,C1,D1;
    A1 = 0;
    B1 = 0;
    C1 = 1;
    D1 = 8+0.5+offset;  // intercept
    rt.setPlane_pp1(A1,B1,C1,D1);
    
    
    
    double A2,B2,C2,D2;
    A2 = 0;
    B2 = 0;
    C2 = -1;
    if (symmetric_walls == true) D2 = 16+0.5+(1-offset);  // intercept:   symmetric walls use 1-offset
    else D2 = 16+0.5+offset;
    rt.setPlane_pp2(A2,B2,C2,D2);
    
    std::cout << "Bottom Plate=" << D1 << " Top Plate=" << D2 << std::endl;
    
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
    
    imin = 1;
    imax = (Nx-1);
    jmin = 1;
    jmax = (Ny-1);
    kmin = 1;
    kmax = (Nz-1);
    
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
                dat = rt.getRayTraceDistance_pp();
                if (dat > 0) LIBB_id[m] = 0;
            }
        }
    }
    
    for (int i=imin;i<imax;i++){
        for (int j=jmin;j<jmax;j++){
            for (int k=kmin;k<kmax;k++){
                m = k*Nx*Ny+j*Nx+i;
                /* ray origin: LBM units */
                ray_origin[0] = double(i);
                ray_origin[1] = double(j);
                ray_origin[2] = double(k-1); // Fixes the meshing issue at the solid-fluid boundary
                
                rt.setRayOrigin(ray_origin);
                dat = rt.getRayTraceDistance_pp();
                if (dat < signedDistance[m]) signedDistance[m] = -dat-0.5;
            }
        }
    }
    
    
    
    
    for (int i=imin;i<imax;i++){
        for (int j=jmin;j<jmax;j++){
            for (int k=kmin;k<kmax;k++){
                m = k*Nx*Ny+j*Nx+i;
                
                /* ray origin: LBM units */
                ray_origin[0] = double(i);
                ray_origin[1] = double(j);
                ray_origin[2] = double(k);
                rt.setRayOrigin(ray_origin);
                
                dat = 0;
                /* Generate ray trace distance map */
                // 1
                ray_direction[0] = 1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= 1 && dat > 0) qDistances[m] = dat;
                dat = 0;
                
                // 2
                ray_direction[0] = -1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= 1 && dat > 0) qDistances[m+Nx*Ny*Nz] = dat;
                dat = 0;
                // 3
                ray_direction[0] = 0.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= 1 && dat > 0) qDistances[m+2*Nx*Ny*Nz] = dat;
                dat = 0;
                // 4
                ray_direction[0] = 0.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= 1 && dat > 0) qDistances[m+3*Nx*Ny*Nz] = dat;
                dat = 0;
                // 5
                ray_direction[0] = 0.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= 1 && dat > 0) qDistances[m+4*Nx*Ny*Nz] = dat;
                dat = 0;
                // 6
                ray_direction[0] = 0.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= 1 && dat > 0) qDistances[m+5*Nx*Ny*Nz] = dat;
                dat = 0;
                // 7
                ray_direction[0] = 1.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+6*Nx*Ny*Nz] = dat;
                dat = 0;
                // 8
                ray_direction[0] = -1.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+7*Nx*Ny*Nz] = dat;
                dat = 0;
                // 9
                ray_direction[0] = 1.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+8*Nx*Ny*Nz] = dat;
                dat = 0;
                // 10
                ray_direction[0] = -1.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = 0.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+9*Nx*Ny*Nz] = dat;
                dat = 0;
                // 11
                ray_direction[0] = 1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = 1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+10*Nx*Ny*Nz] = dat;
                dat = 0;
                // 12
                ray_direction[0] = -1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+11*Nx*Ny*Nz] = dat;
                dat = 0;
                // 13
                ray_direction[0] = 1.0;
                ray_direction[1] = 0.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+12*Nx*Ny*Nz] = dat;
                dat = 0;
                // 14
                ray_direction[0] = -1;
                ray_direction[1] = 0.0;
                ray_direction[2] = 1;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+13*Nx*Ny*Nz] = dat;
                dat = 0;
                // 15
                ray_direction[0] = 0.0;
                ray_direction[1] = 1;
                ray_direction[2] = 1;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+14*Nx*Ny*Nz] = dat;
                dat = 0;
                // 16
                ray_direction[0] = 0.0;
                ray_direction[1] = -1;
                ray_direction[2] = -1;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                dat = abs(dat);
                if (dat <= sqrt(2) && dat > 0) qDistances[m+15*Nx*Ny*Nz] = dat;
                dat = 0;
                // 17
                ray_direction[0] = 0.0;
                ray_direction[1] = 1.0;
                ray_direction[2] = -1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+16*Nx*Ny*Nz] = dat;
                dat = 0;
                // 18
                ray_direction[0] = 0.0;
                ray_direction[1] = -1.0;
                ray_direction[2] = 1.0;
                rt.setRayDirection(ray_direction);
                dat = rt.getRayTraceDistance_pp_q();
                if (dat <= sqrt(2) && dat > 0) qDistances[m+17*Nx*Ny*Nz] = dat;
                dat = 0;
            }
        }
    }
    
}




double RayTrace::ray_trace_st(double *  vec1 /* ray_origin */,
                              double *  vec2 /* plane_pp1 */, double D_st1,
                              double *  vec3 /* plane_pp2 */, double D_st2,
                              double *  vec4 /* plane_st3 */, double D_st3,
                              double *  vec5 /* plane_st4 */, double D_st4) {
    
    double tmp1 = mag(vec2); if (tmp1 == 0) tmp1 = 1.0;
    double tmp2 = mag(vec3); if (tmp2 == 0) tmp2 = 1.0;
    double tmp3 = mag(vec4); if (tmp3 == 0) tmp3 = 1.0;
    double tmp4 = mag(vec5); if (tmp4 == 0) tmp4 = 1.0;
    
    dp1 = dot(vec1, vec2);
    dp2 = dot(vec1, vec3);
    dp3 = dot(vec1, vec4);
    dp4 = dot(vec1, vec5);
    
    
    
    double distance1 = (dp1 - D_st1)/tmp1;
    double distance2 = (dp2 - D_st2)/tmp2;
    double distance3 = (dp3 - D_st3)/tmp3;
    double distance4 = (dp4 - D_st4)/tmp4;
    
    double max12 = std::max(distance1,distance3);
    double max34 = std::max(distance2,distance4);
    
    return std::max(max12,max34);
}

double RayTrace::ray_trace_st_q(double *  vec1 /* ray_origin */,
                                double *  vec2 /* plane_st1 */, double D_st1,
                                double *  vec3 /* plane_st2 */, double D_st2,
                                double *  vec4 /* plane_st3 */, double D_st3,
                                double *  vec5 /* plane_st4 */, double D_st4) {
    
    double tmp1 = mag(vec2); if (tmp1 == 0) tmp1 = 1.0;
    subtract_vectors(vec2 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    o_minus_c[0] += 0; o_minus_c[1] -= (D_st1-1); ; o_minus_c[2] += 0;   // NEGATIVE PLANE
    double * normal; normal = new double[3]; normal[0] = vec2[0]/=tmp1; normal[1] = vec2[1]/=tmp1; normal[2] = vec2[2]/=tmp1;
    tmp1 = dot(ray_direction,normal); if (tmp1 == 0) tmp1 = 1.0;
    dp1 = abs(dot(o_minus_c, normal) / tmp1);
    
    
    
    
    o_minus_c[0] = 0; o_minus_c[1] = 0; o_minus_c[2] = 0;
    normal[0]    = 0; normal[1]    = 0; normal[2]    = 0;
    
    
    
    
    double tmp2 = mag(vec3); if (tmp2 == 0) tmp2 = 1.0;
    subtract_vectors(vec3 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    o_minus_c[0] += 0; o_minus_c[1] += (D_st2-1); o_minus_c[2] += 0;   // POSITIVE PLANE
    normal[0] = vec3[0]/=tmp2; normal[1] = vec3[1]/=tmp2; normal[2] = vec3[2]/=tmp2;
    tmp2 = dot(ray_direction,normal); if (tmp2 == 0) tmp2 = 1.0;
    dp2 = abs(dot(o_minus_c, normal) / tmp2);
    
    
    
    double minA = std::min(dp1,dp2);
    
    
    
    o_minus_c[0] = 0; o_minus_c[1] = 0; o_minus_c[2] = 0;
    normal[0]    = 0; normal[1]    = 0; normal[2]    = 0;
    
    
    
    
    double tmp3 = mag(vec4); if (tmp3 == 0) tmp3 = 1.0;
    subtract_vectors(vec4 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    o_minus_c[0] -= (D_st3-1); o_minus_c[1] += 0; o_minus_c[2] +=0;   // POSITIVE PLANE
    normal[0] = vec4[0]/=tmp3; normal[1] = vec4[1]/=tmp3; normal[2] = vec4[2]/=tmp3;
    tmp3 = dot(ray_direction,normal); if (tmp3 == 0) tmp3 = 1.0;
    dp3 = abs(dot(o_minus_c, normal) / tmp3);
    
    
    
    
    o_minus_c[0] = 0; o_minus_c[1] = 0; o_minus_c[2] = 0;
    normal[0]    = 0; normal[1]    = 0; normal[2]    = 0;
    
    
    
    
    double tmp4 = mag(vec5); if (tmp4 == 0) tmp4 = 1.0;
    subtract_vectors(vec5 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    o_minus_c[0] += (D_st4-1); o_minus_c[1] += 0; o_minus_c[2] += 0;   // POSITIVE PLANE
    normal[0] = vec5[0]/=tmp4; normal[1] = vec5[1]/=tmp4; normal[2] = vec5[2]/=tmp4;
    tmp4 = dot(ray_direction,normal); if (tmp4 == 0) tmp4 = 1.0;
    dp4 = abs(dot(o_minus_c, normal) / tmp4);
    
    
    
    
    double minB = std::min(dp3,dp4);
    
    return std::min(minA,minB);
}


void ComputeLIBBqDistances_st(char * LIBB_id,
                              double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm) {}

void ComputeLIBBqDistances_ct(char * LIBB_id,
                              double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm) {
    // NYI
}

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


void ComputeHWBBqDistances(char * sid, double * qDistances, int & Nx, int & Ny, int & Nz, size_t & count) {
    size_t n,nn,nn2,N;
    int i,j,k;
    N = Nx*Ny*Nz;
    
    size_t c = 0; // For compressed qDistance array
    /*
     You can pre-compute some q's based upon the neighbors, so
     go get your periodic BC code and run through all pairs and populate qDistance
     if both left and right neighbors are solid
     */
    for (n=0; n<N; n++){  // Maybe this is not correct - not sure yet.
        if (!(sid[n] > 0)) {
            k = n/(Nx*Ny);
            j = (n-Nx*Ny*k)/Nx;
            i = n-Nx*Ny*k-Nx*j;
            
            
            nn = n-1;                               // neighbor index (get convention)
            if (i-1<0)        nn += Nx;            // periodic BC along the x-boundary
            //........................................................................
            nn2 = n+1;                              // neighbor index (get convention)
            if (!(i+1<Nx))    nn2 -= Nx;            // periodic BC along the x-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+0*count] = 1.0;  qDistances[c+1*count] = 1.0;  }
            
            
            nn = n-Nx;                            // neighbor index (get convention)
            if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
            //........................................................................
            nn2 = n+Nx;                            // neighbor index (get convention)
            if (!(j+1<Ny))    nn2 -= Nx*Ny;        // Perioidic BC along the y-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+2*count] = 1.0;  qDistances[c+3*count] = 1.0;  }
            
            
            nn = n-Nx*Ny;                        // neighbor index (get convention)
            if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n+Nx*Ny;                        // neighbor index (get convention)
            if (!(k+1<Nz))    nn2 -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+4*count] = 1.0;  qDistances[c+5*count] = 1.0;  }
            
            
            nn = n-Nx-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (j-1<0)            nn += Nx*Ny;    // Perioidic BC along the y-boundary
            //........................................................................
            nn2 = n+Nx+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (!(j+1<Ny))        nn2 -= Nx*Ny;    // Perioidic BC along the y-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+6*count] = 1.0;  qDistances[c+7*count] = 1.0;  }
            
            
            nn = n+Nx-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (!(j+1<Ny))        nn -= Nx*Ny;    // Perioidic BC along the y-boundary
            //........................................................................
            nn2 = n-Nx+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (j-1<0)            nn2 += Nx*Ny;    // Perioidic BC along the y-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+8*count] = 1.0;  qDistances[c+9*count] = 1.0;  }
            
            
            nn = n-Nx*Ny-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (k-1<0)            nn += Nx*Ny*Nz;    // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n+Nx*Ny+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (!(k+1<Nz))        nn2 -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+10*count] = 1.0;  qDistances[c+11*count] = 1.0;  }
            
            
            nn = n+Nx*Ny-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (!(k+1<Nz))        nn -= Nx*Ny*Nz;    // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n-Nx*Ny+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (k-1<0)            nn2 += Nx*Ny*Nz;    // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+12*count] = 1.0;  qDistances[c+13*count] = 1.0;  }
            
            
            nn = n-Nx*Ny-Nx;                    // neighbor index (get convention)
            if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
            if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n+Nx*Ny+Nx;                    // neighbor index (get convention)
            if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
            if (!(k+1<Nz))    nn2 -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+14*count] = 1.0;  qDistances[c+15*count] = 1.0;  }
            
            
            nn = n+Nx*Ny-Nx;                    // neighbor index (get convention)
            if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
            if (!(k+1<Nz))    nn -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n-Nx*Ny+Nx;                    // neighbor index (get convention)
            if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
            if (k-1<0)        nn2 += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+16*count] = 1.0;  qDistances[c+17*count] = 1.0;  }
            
            
            
        }
        c++; // For compressed qDistance array
    }
    
    
}


inline void save_3Ddouble_parallel(int rank,
                                   double* field, int Nx, int Ny, int Nz, int depth, std::string s) {
    //assert(depth >= 0);
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    FILE *f2 = fopen(LocalRankFilename, "w");
    
    double var = 0;
    size_t c = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z;
                
                var = field[c]; c++;
                fprintf(f2, "%d %d %d %.2f\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}
