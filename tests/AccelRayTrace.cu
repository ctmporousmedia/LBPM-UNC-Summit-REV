/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <cuda_profiler_api.h>
#include <algorithm>
#include <cmath>

//declare pointers for device
//allocate
//copy
//copy back


using namespace std;

#define NBLOCKS 1024
#define NTHREADS 256


/*
RayTrace Code
*/

__global__ void dvc_IDFunction(int imax, int imin,
                                int jmax, int jmin,
                                int kmax, int kmin,
                                int Nx, int Ny, int Nz,
                                double sphere_radius,
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

__global__ void dvc_SignedDistanceFunction(int p,
                                            double hx,
                                            int Nx, int Ny, int Nz,
                                            double min_x, double min_y, double min_z,
                                            double* ccx, double* ccy, double* ccz,
                                            double* r,
                                            double* signedDistance)
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
        sphere_center[0] = cx/hx+0.5+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5+0.5;
        sphere_center[2] = cz/hx+0.5+0.5;
        double sphere_radius = rad/hx;
        
        int m = k*Nx*Ny+j*Nx+i;
        double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
        if (dist < signedDistance[m]) signedDistance[m] = dist;
    }//end of idx loop
}

__global__ void dvc_SignedDistanceFunction2(int p,
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


__global__ void dvc_QDistancesFunction(int imax, int imin,
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

extern "C" void IDFunction(int imax, int imin, int jmax, int jmin, int kmax, int kmin, int Nx, int Ny, int Nz, double sphere_radius, double sphere_center_x, double sphere_center_y, double sphere_center_z, char * LIBB_id)
{
    dvc_IDFunction<<<NBLOCKS,NTHREADS >>>(imax, imin, jmax, jmin, kmax, kmin, Nx, Ny, Nz, sphere_radius, sphere_center_x, sphere_center_y, sphere_center_z, LIBB_id);
}

extern "C" void SignedDistanceFunction(int p, double hx, int Nx, int Ny, int Nz, double min_x, double min_y, double min_z, double* ccx, double* ccy, double* ccz, double* r, double* signedDistance)
{
    dvc_SignedDistanceFunction<<<NBLOCKS,NTHREADS >>>(p, hx, Nx, Ny, Nz, min_x, min_y, min_z, ccx, ccy, ccz, r, signedDistance);
}

extern "C" void SignedDistanceFunction2(int p, double hx, int Nx, int Ny, int Nz, double min_x, double min_y, double min_z, double* ccx, double* ccy, double* ccz, double* r, double* signedDistance)
{
    dvc_SignedDistanceFunction2<<<NBLOCKS,NTHREADS >>>(p, hx, Nx, Ny, Nz, min_x, min_y, min_z, ccx, ccy, ccz, r, signedDistance);
}

extern "C" void QDistancesFunction(int imax, int imin, int jmax, int jmin, int kmax, int kmin, int Nx, int Ny, int Nz, double sphere_radius, double sphere_center_x, double sphere_center_y, double sphere_center_z, double* qDistances)
{
    dvc_QDistancesFunction<<<NBLOCKS,NTHREADS >>>(imax, imin, jmax, jmin, kmax, kmin, Nx, Ny, Nz, sphere_radius, sphere_center_x, sphere_center_y, sphere_center_z, qDistances);
}
